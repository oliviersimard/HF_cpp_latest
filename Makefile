SRC=$(PWD)/src
OBJ=$(PWD)/obj

CXX=g++ -std=c++11 -Wall -DONED
INC=${HOME}/include
LIB=${HOME}/lib
CXXFLAGS= -L$(LIB) -larmadillo 

PROG=PiTon

all: $(PROG)

$(PROG): $(OBJ)/main.o $(OBJ)/green_utils.o $(OBJ)/susceptibility_utils.o
	$(CXX) $(OBJ)/main.o $(OBJ)/green_utils.o $(OBJ)/susceptibility_utils.o $(CXXFLAGS) -o main.out

$(OBJ)/main.o: $(PWD)/main.cpp
	$(CXX) -c -I$(INC) $(PWD)/main.cpp -o $(OBJ)/main.o

$(OBJ)/green_utils.o: $(SRC)/green_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/green_utils.cpp -o $(OBJ)/green_utils.o

$(OBJ)/susceptibility_utils.o: $(SRC)/susceptibility_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/susceptibility_utils.cpp -o $(OBJ)/susceptibility_utils.o

clean:
	rm $(OBJ)/* $(PWD)/main.out
