#!/bin/bash


make clean
./configure --with-problem=convection --with-gas=hydro --enable-conduction --enable-ghost
make all
