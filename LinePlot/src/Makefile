CC = g++
OPT = -fPIC -O3 -march=native
CXXFLAGS = $(OPT) -std=c++11 -pedantic -Wall -Wshadow
LDFLAGS = $(OPT) -std=c++11 -pedantic -Wall -Wshadow

INCPATH=../../libgalois/include/galois
SRCPATH=../src
OBJPATH=../obj
LIBPATH=../../libgalois/lib
BINPATH=../bin

INC=
SRC=$(SRCPATH)/main.cpp
OBJ=$(OBJPATH)/main.o
OUT=$(BINPATH)/LinePlot

INCLUDES= -I $(INCPATH)
LIBRARIES= -L $(LIBPATH)
LIB= -lgalois 


default: $(OUT)

$(OUT): $(OBJ)
		$(CC) $(LDFLAGS) $(LIBRARIES) -o $@ $^ $(LIB)

$(OBJPATH)/%.o: $(SRCPATH)/%.cpp
		$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


clean:
		rm -f $(OBJPATH)/*.o
