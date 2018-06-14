all: use_vectors use_matrices use_gmres solve_transport

#For debugging
# OPT=-g -Wall
#For optimistaion
OPT=-O2

#All objects (except main) come from cpp and hpp 
%.o:	%.cpp %.hpp
	g++ ${OPT} -c -o $@ $<
#use_vectors relies on objects which rely on headers
use_vectors:	use_vectors.cpp Vector.o Exception.o
		g++ ${OPT} -o use_vectors use_vectors.cpp Vector.o Exception.o

use_matrices: use_matrices.cpp Matrix.o Exception.o
		g++ ${OPT} -o use_matrices use_matrices.cpp Vector.o Matrix.o Exception.o

use_gmres: use_gmres.cpp Exception.o
	g++ ${OPT} -o use_gmres use_gmres.cpp Vector.o Matrix.o Exception.o

solve_transport: solve_transport.cpp transport.o Exception.o
	g++ ${OPT} -o solve_transport solve_transport.cpp transport.o Vector.o Matrix.o Exception.o


clean:
	rm -f *.o *~ use_vectors, use_matrices, use_gmres