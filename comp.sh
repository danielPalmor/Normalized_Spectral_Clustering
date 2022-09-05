#!/bin/bash
# Script to compile and execute a c program
gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c vectorOperations.c kmeans.c matrixOperations.c -lm -o spkmeans
