#include <math.h>
#include <stdio.h>

#define TOL 0.011

struct CSR
{
	int n;
	float *nonzero;
	int *col;
	int *row_idx;
};

CSR to_CSR(float *A, int n)
{

	int nz = 0;
	for (int i = 0; i < n * n; i++) {
		if(A[i] != 0.0f) {
			nz++;
		}
	}


	float *B = (float *) malloc(sizeof(float) * nz);
	int *cols = (int *) malloc(sizeof(int) * nz);

	int *row_idx = (int *) malloc(sizeof(int) * n);

	nz = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (A[i*n + j] != 0.0f) {
				B[nz] = A[i*n + j];
				cols[nz] = j;
				if (A[i *n + j] > 0.001f) {
					if (i != j) {
						printf("Q is broken, positive not on identity\n");
					}
				}

				if (A[i * n + j] < -0.001f) {
					if (A[i * n + j] != A[j * n + i]) {
						printf("Q is broken, not symmetrical, i = %d, j = %d\n", i, j);
					}
				}
				nz++;
			}
		}

		row_idx[i] = nz;
	}

	CSR csr;
	csr.n = n;
	csr.nonzero = B;
	csr.col = cols;
	csr.row_idx = row_idx;

/*	
	for (int i = 0; i < 12; i++) {
		printf("%5.2f, ", B[i]);
	}
	printf("\n");

	for (int i = 0; i < 12; i++) {
		printf("%5d, ", cols[i]);
	}
	printf("\n");

	for (int i = 0; i < 5; i++) {
		printf("%d, ", row_idx[i]);
	}
	printf("\n");
	printf("made csr\n");
*/

	return csr;
}

float *matcol_csr(CSR csr, float *b, int n)
{
	float *prod = (float *) malloc(sizeof(float) * n);

	for (int i = 0; i < n; i++) {
		float sum = 0;
		for (int j = (i == 0 ? 0 : csr.row_idx[i-1]); j < csr.row_idx[i]; j++) {
			sum += csr.nonzero[j] * b[csr.col[j]];
		}
		prod[i] = sum;
	}

	return prod;
}

int free_csr(CSR csr)
{
	free(csr.nonzero);
	free(csr.col);
	free(csr.row_idx);
	return 0;
}


float *matcol(float *A, float *b, int n)
{
	float *prod = (float *) malloc(sizeof(float) * n);

	for (int i = 0; i < n; i++) {
		float sum = 0;
		for (int j = 0; j < n; j++) {
			sum += b[j] * A[n * i + j];
			//printf("multiplying %4.2f and %4.2f\n", b[j], A[n * i + j]); 
		}
		prod[i] = sum;
	}

	return prod;
}

float *a_psb(float *a, float *b, int n, float s)
{
	float *sum = (float *) malloc(sizeof(float) * n);

	for (int i = 0; i < n; i++) {
		sum[i] = a[i] + s * b[i];
	}

	return sum;
}

float dot(float *a, float *b, int n)
{
	float sum = 0;
	for (int i = 0; i < n; i++) {
		sum += a[i] * b[i];
	}
	return sum;
}

float *conj_grad(float *A, float *b, float *x0, int n)
{
	//returns solution to A * x = b

	CSR A_csr = to_CSR(A, n);

	//float *Ap = matcol(A, x0, n);
	float *Ap = matcol_csr(A_csr, x0, n);
	//exit(1);
	float *r = a_psb(b, Ap, n, -1);
	free(Ap);

	float *p = a_psb(r, x0, n, 0);

	float rsold = dot(r, r, n);


	float *x = a_psb(x0, r, n, 0);

	int run = 0;
	while (rsold > TOL * TOL) {
		//Ap = matcol(A, p, n);
		Ap = matcol_csr(A_csr, p, n);

		float alpha = rsold / dot(p, Ap, n);

		float *xtemp = x;
		x = a_psb(x, p, n, alpha);
		free(xtemp);

		float *rtemp = r;
		r = a_psb(r, Ap, n, -alpha);
		//float *ax_temp = matcol_csr(A_csr, x, n);
		//r = a_psb(b, ax_temp, n, -1);
		//free(ax_temp);
		free(rtemp);
		free(Ap);

		float rsnew = dot(r, r, n);

		float *ptemp = p;
		p = a_psb(r, p, n, rsnew / rsold);
		free(ptemp);
		rsold = rsnew;

		//printf("run = %d, n = %d, rsold = %10.3f\n", run++, n, sqrt(rsold));
	}
	//printf("\n");

	free(r);
	free(p);
	free_csr(A_csr);

	return x;
}
