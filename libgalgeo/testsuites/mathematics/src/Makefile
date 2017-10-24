CC = g++
OPT =
CXXFLAGS = $(OPT) -std=c++11 -pedantic -Wall
LDFLAGS = $(OPT) -std=c++11 -pedantic -Wall

INCPATH=../../include/galgeo
SRCPATH=../src
OBJPATH=../obj
LIBPATH=$(PWD)/../../../lib/
BINPATH=../bin

INC=../inc
SRC=$(SRCPATH)/main.cpp
OBJ=$(OBJPATH)/main.o
OUT=$(BINPATH)/test

INCLUDES= -I $(INCPATH)
LIBRARIES= -L $(LIBPATH)
LIB=-lgalgeo


default: $(OUT)

$(OUT): $(OBJ)
		$(CC) $(LDFLAGS) $(LIBRARIES) -o $@ $^ $(LIB)

$(OBJPATH)/%.o: $(SRCPATH)/%.cpp
		$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


clean:
		rm -f $(OBJPATH)/*.o
