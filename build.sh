#!/bin/bash

dir="${PWD}/obj"

echo $dir

if [ ! -d $dir ]; then
        mkdir ${PWD}/obj
else
        echo -e "Directory $dir already created\n"
fi

make
