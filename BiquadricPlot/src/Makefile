CC = g++
OPT = -fPIC -O3 -march=native
CXXFLAGS = $(OPT) -std=c++11 -pedantic -Wall -Wshadow
LDFLAGS = $(OPT) -std=c++11 -pedantic -Wall -Wshadow

INCPATH=../../libgalgeo/include/galgeo
SRCPATH=../src
OBJPATH=../obj
LIBPATH=../../libgalgeo/lib
BINPATH=../bin

INC=
SRC=$(SRCPATH)/main.cpp
OBJ=$(OBJPATH)/main.o
OUT=$(BINPATH)/BiquadricPlot

INCLUDES= -I $(INCPATH)
LIBRARIES= -L $(LIBPATH)
LIB= -lgalgeo


default: $(OUT)

$(OUT): $(OBJ)
		$(CC) $(LDFLAGS) $(LIBRARIES) -o $@ $^ $(LIB)

$(OBJPATH)/%.o: $(SRCPATH)/%.cpp
		$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


clean:
		rm -f $(OBJPATH)/*.o
