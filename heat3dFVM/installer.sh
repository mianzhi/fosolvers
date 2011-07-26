#!/bin/bash

# installer of heat3dFVM solver
sudo cp -v heat3dFVM.sh /usr/bin/heat3dFVM
sudo chmod +x /usr/bin/heat3dFVM
if [ ! -r $HOME/.FOSolverS ]; then
  mkdir -v $HOME/.FOSolverS
fi
if [ ! -r $HOME/.FOSolverS/heat3dFVM ]; then
  mkdir -v $HOME/.FOSolverS/heat3dFVM
fi
cp -v heat3dFVM $HOME/.FOSolverS/heat3dFVM/
cp -v materials.mtl $HOME/.FOSolverS/heat3dFVM/
chmod +x $HOME/.FOSolverS/heat3dFVM/heat3dFVM
echo "done!"
