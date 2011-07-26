#!/bin/bash

# Starter of se3dFEM solver
workspace=$(pwd)
if [ ! -r conditions.cod ]; then
  echo 'ERROR: can not find (or read) "conditions.cod".'
elif [ ! -r grid.msh ]; then
  echo 'ERROR: can not find (or read) "grid.msh".'
else
  if [ ! -r se3dFEM ]; then
    cp $HOME/.FOSolverS/se3dFEM/se3dFEM ${workspace}
  fi
  if [ ! -r materials.mtl ]; then
    cp $HOME/.FOSolverS/se3dFEM/materials.mtl ${workspace}
  fi
  ./se3dFEM
fi
