SRC=$(PWD)/src
OBJ=$(PWD)/obj

CXX=g++ -std=c++14 -g -O0 -Wall #-Xpreprocessor -fopenmp #-DONED
INC=${HOME}/include
LIB=${HOME}/lib
CXXFLAGS= -L$(LIB) -larmadillo -lfftw3 -ljson_spirit -lpthread #-lboost_system -lboost_thread#-lomp

PROG=PiTon

all: $(PROG)

$(PROG): $(OBJ)/main_plot_ladder_parts.o $(OBJ)/green_utils.o $(OBJ)/susceptibility_utils.o $(OBJ)/fft.o $(OBJ)/json_utils.o $(OBJ)/integral_utils.o $(OBJ)/thread_utils.o
	$(CXX) $(OBJ)/main_plot_ladder_parts.o $(OBJ)/green_utils.o $(OBJ)/susceptibility_utils.o $(OBJ)/fft.o $(OBJ)/json_utils.o $(OBJ)/integral_utils.o $(OBJ)/thread_utils.o $(CXXFLAGS) -o main_plot_ladder_parts.out

$(OBJ)/main_plot_ladder_parts.o: $(PWD)/main_plot_ladder_parts.cpp
	$(CXX) -c -I$(INC) $(PWD)/main_plot_ladder_parts.cpp -o $(OBJ)/main_plot_ladder_parts.o

$(OBJ)/green_utils.o: $(SRC)/green_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/green_utils.cpp -o $(OBJ)/green_utils.o

$(OBJ)/susceptibility_utils.o: $(SRC)/susceptibility_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/susceptibility_utils.cpp -o $(OBJ)/susceptibility_utils.o

$(OBJ)/fft.o: $(SRC)/fft.cpp
	$(CXX) -c -I$(INC) $(SRC)/fft.cpp -o $(OBJ)/fft.o

$(OBJ)/json_utils.o: $(SRC)/json_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/json_utils.cpp -o $(OBJ)/json_utils.o

$(OBJ)/integral_utils.o: $(SRC)/integral_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/integral_utils.cpp -o $(OBJ)/integral_utils.o

$(OBJ)/thread_utils.o: $(SRC)/thread_utils.cpp
	$(CXX) -c -I$(INC) $(SRC)/thread_utils.cpp -o $(OBJ)/thread_utils.o

clean:
	rm $(OBJ)/* $(PWD)/main_plot_ladder_parts.out
