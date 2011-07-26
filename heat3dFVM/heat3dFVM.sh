#!/bin/bash

# Starter of heat3dFVM solver
workspace=$(pwd)
if [ ! -r conditions.cod ]; then
  echo 'ERROR: can not find (or read) "conditions.cod".'
elif [ ! -r grid.msh ]; then
  echo 'ERROR: can not find (or read) "grid.msh".'
else
  if [ ! -r heat3dFVM ]; then
    cp $HOME/.FOSolverS/heat3dFVM/heat3dFVM ${workspace}
  fi
  if [ ! -r materials.mtl ]; then
    cp $HOME/.FOSolverS/heat3dFVM/materials.mtl ${workspace}
  fi
  ./heat3dFVM
fi
