CC = clang

CFLAGS = -Wall -Wextra -std=gnu2x -O3 -march=native -g -fopenmp
LDLIBS = -lm


all: smooth-ref smooth-opt
