placer: src/placer.cpp src/conj_grad.c
	g++ -o placer src/placer.cpp src/suraj_parser.cpp src/conj_grad.c

zip:
	tar -cvf macle119_analytic_placement.tar.gz src Makefile visualizer.py benchmarks README.md

clean:
	rm placer *.txt *.png
