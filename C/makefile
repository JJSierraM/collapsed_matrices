CC=gcc
LIBCODE=vector2D.c body.c
DEBUG=-g

build: 
	$(CC) -O3 -march=native $(LIBCODE) -lm collapsed_matrices.c -o collapsed_matrices -fopenmp

exec:
	./collapsed_matrices

debug-build:
	$(CC) $(DEBUG) -O0 -march=native $(LIBCODE) -lm main.c -o main -fopenmp

valgrind:
	make debug-build
	valgrind ./collapsed_matrices --leak-check=full