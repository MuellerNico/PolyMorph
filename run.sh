#!/bin/bash

OMP_NUM_THREADS=4 # <- change this to available cores
output="/mnt/c/Users/muell/Desktop/PolymorphOutput" # <- change this to your desired output folder

move_files() {
    echo "moving output files ..."
    local timestamp=$(date +%Y-%m-%d_%H-%M)
    mkdir $output/$timestamp
    mv /out/*.vtp /out/*.vts /out/*.cfg /out/log.txt /out/*.off "$output/$timestamp/"
}

cleanup_files() {
    echo "cleaning up leftover output files ..."
    rm /out/*.vtp /out/*.vts /out/*.cfg /out/log.txt /out/*.off
    make clean
}

run() {
    make clean
    make
    export OMP_NUM_THREADS
    echo "running PolyMorph with $OMP_NUM_THREADS threads... "
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