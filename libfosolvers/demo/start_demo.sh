#!/bin/sh

mpirun -np 2 demo --grid grid.msh --result rst.msh --data data.tab --condition conditions.cod
