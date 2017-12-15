#!/bin/bash


make clean
./configure --with-problem=convection --with-flux=roe --with-gas=hydro --enable-conduction --enable-ghost $1 $2
make all
