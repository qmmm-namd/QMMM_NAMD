#!/usr/bin/env bash 

cpptraj << EOF 
parm prmtop 
trajin coord.nc mdvel vel.nc 
trajout test.ncrst
run 
exit 
EOF


