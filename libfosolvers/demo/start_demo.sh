#!/bin/sh

mpirun -np 2 demo --grid grid.msh --rst rst.msh --data data.tab
