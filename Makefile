CC = mpiCC -std=c++0x -O3
CFLAGS= -Wall -fopenmp -c
LFLAGS= -Wall -fopenmp

INC=-I/usr/local/include -I/usr/local/Healpix_3.31/src/cxx/generic_gcc/include
LIB=-L/usr/local/lib -L/usr/local/Healpix_3.31/src/cxx/generic_gcc/lib -lhealpix_cxx -lcxxsupport -lsharp -lfftpack -lc_utils -lcfitsio -lgsl -lgslcblas -lm

BINDIR = bin
MAINDIR = src/main
HEADDIR = src/include
SSDIR = src/mods

SRC=$(wildcard $(SSDIR)/*.cpp $(MAINDIR)/*.cpp)
OBJ=$(patsubst $(SSDIR)/%.cpp, $(BINDIR)/%.o, $(SRC))

ALL = Assassi

all: $(ALL)

### Compute angular power spectrum and bispectrum predictions
Assassi: $(OBJ)
	$(CC) $(LFLAGS) $(INC) $(OBJ) $(LIB) -o exec.out

$(BINDIR)/%.o: $(SSDIR)/%.cpp $(HEADDIR)/%.h
	$(CC) $(CFLAGS) $(INC) $< -o $@ $(LIB)

clean:
	rm -f $(BINDIR)/*
