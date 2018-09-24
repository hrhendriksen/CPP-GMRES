all: use_vectors use_matrices use_gmres solve_transport
#For debugging
# FLAGS=-g -Wall
#For optimistaion
FLAGS=-O2 -std=c++11

#All objects (except main) come from cpp and hpp 
%.o:	%.cpp %.hpp
	g++ ${FLAGS} -c -o $@ $<
#use_vectors relies on objects which rely on headers
use_vectors:	use_vectors.cpp Vector.o Exception.o
		g++ ${FLAGS} -o use_vectors use_vectors.cpp Vector.o Exception.o

use_matrices: use_matrices.cpp Matrix.o Exception.o
		g++ ${FLAGS} -o use_matrices use_matrices.cpp Vector.o Matrix.o Exception.o

use_gmres: use_gmres.cpp Givens.o Exception.o
	g++ ${FLAGS} -o use_gmres use_gmres.cpp Givens.o Vector.o Matrix.o Exception.o

solve_transport: solve_transport.cpp transport.o Givens.o Exception.o
	g++ ${FLAGS} -o solve_transport solve_transport.cpp transport.o Givens.o Vector.o Matrix.o Exception.o


clean:
	rm -f *.o *~ use_vectors,use_matrices,use_gmres,solve_transport