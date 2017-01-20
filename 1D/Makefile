PROG =	run

SRCS =	es.f90 main.f90 numerique.f90 donnees.f90

OBJS =  main.o es.o numerique.o donnees.o

LIBS =	

F90 = gfortran
#F90FLAGS = -O3 -march=native
F90FLAGS = -O0 -pedantic -Wall -ffpe-trap=invalid,zero,overflow,underflow -g -fbounds-check -fbacktrace -fdump-core 
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

es.o: donnees.o numerique.o
main.o: es.o numerique.o donnees.o
numerique.o: donnees.o
donnees.o: 
