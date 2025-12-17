# Analytic Placement
## Alexis MacLean

This is a project for EE 5301 involving implementing part of the algorithm described in the paper "FastPlace: Efficient Analytical Placement using Cell Shifting, Iterative Local Refinement and a Hybrid Net Model"

I may add more features to this project after the semester, but as of now what is implemented is as required by the assigment.

**ASSUMPTIONS:** I used the canvas announcements list of assumptions to design the spreading part of the algorithm.
That is to say, I consider cells as points and calculate utilization as such.

## Usage

`$ make` or `$ make placer` will compile the binary `placer`.

`$ ./placer benchmarks/toy01_02/toy01/toy01` will run the placer on the first toy example problem, producing the files `pre_spreading.txt`, `post_spreading.txt`, and `wire_length.txt` in the main directory.
These files will not be labeled according to their example, and will be overwritten with each run.

`$ python visualizer.py` will generate an image `visualizer_plot.png` based on the pre and post spreading locations of the cells.

`$ make clean` will delete the binary file and all text and image files in the main directory.

`

## File descriptions

`src/conj_grad.c/h` describe an implementation of the conjugate gradient method to optimize a sparse quadratic problem.
I implemented this using a CSR format for performance.

`src/placer.cpp` is the main body of the placement, it has functions for generating the $Q$ matrix, spreading, and writing out to the files.

`visualizer.py` is a short python script which uses matplotlib to generate a visualization of the spreading step of the FastPlace algorithm.

`src/suraj_parser.cpp/h` is a parser for the provided filetypes and was provided as a part of the assignment.

`benchmarks/*` folder of provided examples to test the placer on
