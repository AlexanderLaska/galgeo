CC = g++
OPT = -fPIC -O3 -march=native
CXXFLAGS = $(OPT) -std=c++11 -pedantic -Wall

HDIR = ../include/galgeo/

ODIR = ../lib/obj/
OBJ = GaloisScalar.o GaloisOneZeroTensor.o GaloisZeroOneTensor.o GaloisOneOneTensor.o GaloisZeroTwoTensor.o GaloisTwoZeroTensor.o FiniteProjectiveGeometry.o
OBJS = $(patsubst %,$(ODIR)%,$(OBJ))

LIBNAME = libgalgeo
OUT = $(PWD)/../lib/


all: library

library: $(OBJS)
	$(CC) -shared $(CXXFLAGS) $(OBJS) -o $(OUT)$(LIBNAME).so

$(ODIR)%.o: %.cpp $(HDIR)%.hpp
	mkdir -p $(ODIR)
	$(CC) $(CXXFLAGS) -c $< -o $@


clean: clear

clear:
	rm -f $(OUT)*.so
	rm -f $(ODIR)*.o
