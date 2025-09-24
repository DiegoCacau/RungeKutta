#!/bin/bash

set -e

CC=g++
CFLAGS=-O3


$CC -o exec.out $CFLAGS teste.cpp

./exec.out

gnuplot -p grafico.p