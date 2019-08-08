#!/bin/bash

dir=${PWD}/obj

SRC=${PWD}/src

CPPFLAGS=-std=c++11

CXX=g++

echo $dir

if [ ! -d $dir ]; then
        mkdir ${dir}
else
        echo -e "Directory $dir already created\n"
fi

make -j 2

# ${CXX} ${CPPFLAGS} -c ${SRC}/green_utils.cpp -o ${dir}/green_utils.o

# ${CXX} ${CPPFLAGS} -c ${SRC}/susceptibility_utils.cpp -o ${dir}/susceptibility_utils.o

# ${CXX} ${CPPFLAGS} -c main.cpp -o ${dir}/main.o

# ${CXX} ${CPPFLAGS} ${dir}/main.o ${dir}/green_utils.o ${dir}/susceptibility_utils.o -larmadillo -o main.out
