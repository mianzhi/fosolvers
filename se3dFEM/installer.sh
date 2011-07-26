#!/bin/bash

# installer of se3dFEM solver
sudo cp -v se3dFEM.sh /usr/bin/se3dFEM
sudo chmod +x /usr/bin/se3dFEM
if [ ! -r $HOME/.FOSolverS ]; then
  mkdir -v $HOME/.FOSolverS
fi
if [ ! -r $HOME/.FOSolverS/se3dFEM ]; then
  mkdir -v $HOME/.FOSolverS/se3dFEM
fi
cp -v se3dFEM $HOME/.FOSolverS/se3dFEM/
cp -v materials.mtl $HOME/.FOSolverS/se3dFEM/
chmod +x $HOME/.FOSolverS/se3dFEM/se3dFEM
echo "done!"
