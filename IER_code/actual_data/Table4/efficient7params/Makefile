FORTRAN = ifort
TARGET  = output.out
FLAGS  = -O3
#FLAGS  = -O3 -ip -xN
#FLAGS   = -O3 -fpstkchk -ftrapuv -check all
#DFLAGS = ${FLAGS} -debug extended
DFLAGS  = ${FLAGS} -ggdb
SOURCES =  mparams.f90 model.f90 model_sol.f90 mempirical.f90 mainMCMC.f90

all: ${SOURCES} sfunctions.o lapack_d.o linpack_z.o blas1_z.o
	${FORTRAN} ${FLAGS} -o ${TARGET} sfunctions.o lapack_d.o linpack_z.o blas1_z.o ${SOURCES}

clean:
	'rm' -f *.o *.13.cfg *.life1 *.life2 *.mod *.out*

debug: ${SOURCES} sfunctions.o lapack_d.o blas1_z.o linpack_z.o
	${FORTRAN} ${DFLAGS} -o ${TARGET} sfunctions.o lapack_d.o linpack_z.o blas1_z.o ${SOURCES}

mainMCMC.o: mainMCMC.f90
	${FORTRAN} ${FLAGS} -c mainMCMC.f90

sfunctions.o: sfunctions.f90
	${FORTRAN} ${FLAGS} -c sfunctions.f90

mempirical.o: mempirical.f90
	${FORTRAN} ${FLAGS} -c mempirical.f90

model.o: model.f90
	${FORTRAN} ${FLAGS} -c model.f90

model_sol.o: model_sol.f90
	${FORTRAN} ${FLAGS} -c model_sol.f90

lapack_d.o: lapack_d.f90
	${FORTRAN} ${FLAGS} -c lapack_d.f90

linpack_z.o: linpack_z.f90
	${FORTRAN} ${FLAGS} -c linpack_z.f90

lapack_s.o: lapack_s.f90
	${FORTRAN} ${FLAGS} -c lapack_s.f90

blas1_z.o: blas1_z.f90
	${FORTRAN} ${FLAGS} -c blas1_z.f90

