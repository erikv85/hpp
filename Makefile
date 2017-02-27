FLAGS=-lm -fopenmp -O0

galsim: galsim.c
	gcc galsim.c -o galsim $(FLAGS)

omp_gal: omp_gal.c
	gcc omp_gal.c -o omp_gal $(FLAGS)
