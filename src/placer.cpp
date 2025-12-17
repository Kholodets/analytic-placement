#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "suraj_parser.h"
#include "conj_grad.h"

#define DELTA 1.5f

//TODO pass the stuff in as arguments instead of using the globals
int make_Qd(float **Q_h, float **dx_h, float **dy_h, int *n_h)
{
	//first, count the number of hyperedges that will be turned into stars
	//this tells us how big q needs to be
	
	int stars = 0;
	for (int i = 0; i < numhyper; i++) {
		if (hEdge_idxToFirstEntryInPinArray[i+1] - hEdge_idxToFirstEntryInPinArray[i] > 2 ) {
			stars++;
		}
	}

	int n = stars + numCells_noPads;
	float *Q = (float *) malloc(sizeof(float) * n * n);
	float *dx = (float *) malloc(sizeof(float) * n);
	float *dy = (float *) malloc(sizeof(float) * n);
	
	memset(Q, 0, n*n * sizeof(float));
	memset(dx, 0, n * sizeof(float));
	memset(dy, 0, n * sizeof(float));

	//going to count stars back down for indicies into Q
	stars--;

	for (int i = 0; i < numhyper; i++) {
		int first = hEdge_idxToFirstEntryInPinArray[i];
		int last = hEdge_idxToFirstEntryInPinArray[i+1] - 1;
		float w = hyperwts[i];
		int k = last + 1 - first;
		float gamma = 1.0f / (((float) k) - 1.0f);

		if (last - first  == 1 ) {
			//normal edge
			int cell1 = cellPinArray[first];
			int cell2 = cellPinArray[last];
			//printf("normal edge, nodes %d and %d\n", cell1, cell2);
			if (cell1 >= numCells_noPads) {
				dx[cell2] += pinLocations[cell1 - numCells_noPads].x * w;
				dy[cell2] += pinLocations[cell1 - numCells_noPads].y * w;
				Q[n * cell2 + cell2] += w;
			} else if (cell2 >= numCells_noPads) {
				dx[cell1] += pinLocations[cell2 - numCells_noPads].x * w;
				dy[cell1] += pinLocations[cell2 - numCells_noPads].y * w;
				Q[n * cell1 + cell1] += w;
			} else {
				Q[n * cell1 + cell2] -= w;
				Q[n * cell2 + cell1] -= w;
			
				Q[n * cell1 + cell1] += w;
				Q[n * cell2 + cell2] += w;
			}

		} else if (last - first == 2 ) {
			w *= gamma;
			//degree 3 hyperedge (clique)
			//printf("clique, nodes ");
			int cells[3];
			int pads[2];
			int ncells = 0;
			int npads = 0;

			for (int i = first; i <= last; i++) {
				int cell = cellPinArray[i];
				//printf("%d, ", cell);
				if (cell >= numCells_noPads) {
					pads[npads] = cell;
					npads++;
				} else {
					cells[ncells] = cell;
					ncells++;
				}
			}
			//printf("\n");

			for (int j = 0; j < ncells; j++) {
				for (int k = j+1; k < ncells; k++) {
					Q[n * cells[j] + cells[k]] -= w;
					Q[n * cells[k] + cells[j]] -= w;

					Q[n * cells[j] + cells[j]] += w;
					Q[n * cells[k] + cells[k]] += w;
				}

				for (int k = 0; k < npads; k++) {
					dx[cells[j]] += pinLocations[pads[k] - numCells_noPads].x * w;
					dy[cells[j]] += pinLocations[pads[k] - numCells_noPads].y * w;
					Q[n * cells[j] + cells[j]] += w;
				}
			}

		} else if (last - first > 2 ) {
			w *= gamma * k;
			//degree >3 hyperedge (star)
			//Q[numCells_noPads + stars] is our star cell
			//printf("star %d, nodes ", stars);
			for (int i = first; i <= last; i++) {
				int cell = cellPinArray[i];
				//printf("%d, ", cell);
				if (cell >= numCells_noPads) {
					dx[stars + numCells_noPads] += pinLocations[cell - numCells_noPads].x * w;
					dy[stars + numCells_noPads] += pinLocations[cell - numCells_noPads].y * w;
					Q[n * (stars + numCells_noPads) + (stars + numCells_noPads)] += w;
				} else {
					Q[n * cell + (stars + numCells_noPads)] -= w;
					Q[n * (stars + numCells_noPads) + cell] -= w;
					
					Q[n * (stars + numCells_noPads) + (stars + numCells_noPads)] += w;
					Q[n * cell + cell] += w;
				}
			}
			//printf("\n");
			stars--;
		}
	}

	*Q_h = Q;
	*dx_h = dx;
	*dy_h = dy;
	*n_h = n;
	return 0;
}

int free_Qd(float *Q, float *dx, float *dy)
{
	free(Q);
	free(dx);
	free(dy);
	return 0;
}

float calc_wirelength(float *Q, float *x, float *y, int n, int stride)
{
	float sum = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			float dx = x[i] - x[j];
			float dy = y[i] - y[j];
			sum -= Q[i * stride + j] * (dx * dx + dy * dy);
		}
	}
	return sqrt(sum);
}

int spread(float *ox, float *oy, float **nx_h, float **ny_h, int n, int all, float a, float max_x, float max_y)
{
	float *nx = (float *) malloc(sizeof(float) * all);
	float *ny = (float *) malloc(sizeof(float) * all);

	float U[25];
	memset(U, 0.0f, sizeof(int) * 25);

	float ar = (max_x / 5.0f) * (max_y / 5.0f);


	
	for (int i = 0; i < n; i++) {
		int xid = (int) (5.0f * ox[i] / max_x);
		int yid = (int) (5.0f * oy[i] / max_y);

		U[yid * 5 + xid] += 1.0f; /* / ar  DO NOT NORMALIZE UTILIZATION*/;
	}

	
	/*
	for (int i = 0; i < 25; i++) {
		printf("%f, ", U[i]);
		if (i % 5 == 4) {
			printf("\n");
		}
	}
	*/	
	

	//perform x spreading
	float xob[6];
	float xnb[6];

	xnb[0] = 0;
	xnb[5] = max_x;



	//printf("xob: ");
	for (int i = 0; i < 6; i++) {
		xob[i] = i * (max_x / 5.0f);
		//printf("%f, ", xob[i]);
	}
	//printf("\n");
	

	for (int i = 0; i < 5; i++) {
		//printf("xnb: 0, ");
		for (int j = 1; j < 5; j++) {
			//compute new x boundaries for current row
			xnb[j] = ( xob[j-1]*(U[i * 5 + j] + DELTA) + xob[j+1]*(U[i * 5 + j - 1] + DELTA) ) 
				/ ( U[i * 5 + j - 1] + U[i * 5 + j] + 2 * DELTA );
			//printf("%f, ", xnb[j]);
		}
		//printf(", %f\n", xnb[5]);

		for (int j = 0; j < n; j++) {
			//printf("i = %d, oy[j] = %f\n", i, oy[j]);
			//move this cell if its in this row
			if ( ((int)(5.0f * oy[j] / max_y)) == i ) {
				float xj = ox[j];
				int ui = ((int) (5.0f * xj / max_x)) + 1;
				float xjp = (xnb[ui] * (xj - xob[ui - 1]) + xnb[ui - 1] * (xob[ui] - xj)) 
					/ (xob[ui] - xob[ui - 1]);
				nx[j] = xj + a * (xjp - xj);
				//printf("xj %f, xjp %f, nx %f\n", xj, xjp, nx[j]);
			}
		}

	}

	//printf("spread x\n");
	
	//perform y spreading
	float yob[6];
	float ynb[6];

	ynb[0] = 0;
	ynb[5] = max_y;

	for (int i = 0; i < 6; i++)
		yob[i] = i * (max_y / 5.0f);

	for (int i = 0; i < 5; i++) {
		//printf("ynb: %f, ", ynb[0]);
		for (int j = 1; j < 5; j++) {
			//compute new y boundaries for current columb
			ynb[j] = ( yob[j-1]*(U[j * 5 + i] + DELTA) + yob[j+1]*(U[(j-1) * 5 + i ] + DELTA) ) 
				/ ( U[(j-1) * 5 + i] + U[j * 5 + i] + 2 * DELTA );
			//printf("%f, ", ynb[j]);
		}
		//printf("%f\n", ynb[5]);

		for (int j = 0; j < n; j++) {
			//move this cell if its in this column
			float yj = oy[j];
			int ui = ((int) (5.0f * yj / max_y));
			if ( ((int)(5.0f * ox[j] / max_x)) == i /*nx[j] >= xnb_save[ui * 6 + i] && nx[j] <= xnb_save[ui * 6 + i + 1] */) {
				ui++;
				float yjp = (ynb[ui] * (yj - yob[ui - 1]) + ynb[ui - 1] * (yob[ui] - yj)) 
					/ (yob[ui] - yob[ui - 1]);
				ny[j] = yj + a * (yjp - yj);
				//printf("yj %f, yjp %f, ny %f\n", yj, yjp, ny[j]);
			}
		}

	}

	//printf("spread y\n");
	*nx_h = nx;
	*ny_h = ny;
	return 0;
}

void print_out(float *x, float *y, int n, FILE *out)
{
	for (int i = 0; i < n; i++) {
		fprintf(out, "c %6.3f %6.3f\n", x[i], y[i]);
	}
}

void print_pads(FILE *out)
{
	for (int i = 0; i < numCellsAndPads - numCells_noPads; i++) {
		SPinLocation spl = pinLocations[i];
		fprintf(out, "p %6.3f %6.3f\n", (float) spl.x, (float) spl.y);
	}
}

int main(int argc, char **argv)
{
	char inareFileName[64], innetFileName[64], inPadLocationFileName[64];

	if (argc < 2) {
		printf("incorrect usage\n");
		return 1;
	}

	//discard first argument, its the name of the binary
	argc--, argv++;

	sprintf(inareFileName, "%s.are", *argv);
	sprintf(innetFileName, "%s.net", *argv);
	sprintf(inPadLocationFileName, "%s.kiaPad", *argv);
	argc--, argv++;

	//printf("%s\n", inareFileName);
	//printf("%s\n", innetFileName);
	//printf("%s\n", inPadLocationFileName);

	parseIbmFile(inareFileName, innetFileName, inPadLocationFileName);

	float *Q, *dx, *dy;
	int n;

	make_Qd(&Q, &dx, &dy, &n);

	
	printf("calculated Q, n = %d\n", n);

	/*
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%6.2f, ", Q[n * i + j]);
		}
		printf("\n");
	}
	
	printf("  dx     dy\n");
	for (int i = 0; i < n; i++) {
		printf("%6.2f   %6.2f\n", dx[i], dy[i]);
	}
	*/
	
	

	float max_x = 0;
	float max_y = 0;

	for (int i = 0; i < numCellsAndPads - numCells_noPads; i++) {
		SPinLocation loc = pinLocations[i];
		if (loc.x > max_x) {
			max_x = loc.x;
		}

		if (loc.y > max_y) {
			max_y = loc.y;
		}
	}

	//printf("max x: %f, max y: %f\n", max_x, max_y);

	float *x0 = (float *) malloc(sizeof(float) * n);
	float *y0 = (float *) malloc(sizeof(float) * n);


	for (int i = 0; i < n; i++) {
		x0[i] = max_x / 2.0;
		y0[i] = max_y / 2.0;
	}

	printf("performing gradient descent\n");
	float *x = conj_grad(Q, dx, x0, n);
	float *y = conj_grad(Q, dy, y0, n);

	free(x0);
	free(y0);

	
/*	
	printf("   x      y\n");
	for (int i = 0; i < n; i++) {
		printf("%6.2f   %6.2f\n", x[i], y[i]);
	}
*/


	FILE *wirelength = fopen("wire_length.txt", "w");
	FILE *pre = fopen("pre_spreading.txt", "w");
	FILE *post = fopen("post_spreading.txt", "w");
	float wl = calc_wirelength(Q, x, y, n, n);


	print_out(x, y, n, pre);
	print_pads(pre);

	fprintf(wirelength, "before spreading: %f\n", wl);

	float *xs;
	float *ys;


	//int spread(float *ox, float *oy, float **nx_h, float **ny_h, int n, int all, float a, float max_x, float max_y)
	spread(x, y, &xs, &ys, numCells_noPads, n, 0.8, max_x, max_y);
	//spread(x, y, &xs, &ys, numCells_noPads, 0.8, max_x, max_y);
	//my intuition says that moving the star nodes would be useful, at least for this exercise
	//however, i can see how you might want to ignore star models
	//and for the sake of the correctness of the numbers of this output, I am not considering star nodes in the
	//spread operation, thus passint numCells_noPads instead of n
	
	//passing a new value all to the spread function to tell it how long to alloc them
	//this copies over the extra
	for (int i = numCells_noPads; i < n; i++) {
		xs[i] = x[i];
		ys[i] = y[i];
	}

	float wls = calc_wirelength(Q, xs, ys, n, n);
	//float wls = calc_wirelength(Q, xs, ys, numCells_noPads, n);

	fprintf(wirelength, "after spreading: %f\n", wls);
	print_out(xs, ys, n, post);
	print_pads(post);
	fclose(wirelength);
	fclose(pre);
	fclose(post);

	free_Qd(Q, dx, dy);
	free(x);
	free(y);
	free(xs);
	free(ys);

	return 0;
}
