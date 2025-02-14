#!/bin/bash

OMP_NUM_THREADS=4 # <- change this to available cores
#output="./out" # <- change this to your desired output folder
#output="/mnt/c/Users/nico/Desktop/PolymorphOutput"
output="/mnt/c/Users/muell/Desktop/PolymorphOutput"

move_files() {
    echo "moving output files to folder ..."
    local timestamp=$(date +%Y-%m-%d_%H-%M)
    mkdir $output/$timestamp
    mv *.vtp *.vts *.cfg log.txt "$output/$timestamp/"
}

cleanup_files() {
    echo "cleaning up leftover output files ..."
    rm *.vtp *.vts *.cfg log.txt #*.off
    make clean
}

run() {
    echo "compiling src/$1.cpp ..."
    make clean
    make
    export OMP_NUM_THREADS
    echo "running $1 with $OMP_NUM_THREADS threads... "
    ./polymorph
}

if [ "$1" = "c" ]; then
    cleanup_files
    
elif [ "$1" = "m" ]; then
    move_files
else
    run
    move_files
fi

echo "all done."