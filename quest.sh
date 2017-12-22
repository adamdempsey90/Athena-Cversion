#!/bin/bash


make clean MACHINE=quest 
./configure --with-problem=convection --with-flux=roe --with-gas=hydro --enable-conduction --enable-mpi --enable-viscosity
make all MACHINE=quest 
