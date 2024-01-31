#!/bin/bash

#PBS -l nodes=node78:ppn=1
#PBS -q eth
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR
export CUDA_VISIBLE_DEVICES=1
#mpirun -np 128  \
cd min
    pmemd.cuda -O -i min.in -o min.out -p prmtop -c \
    min.rst7 -r min.ncrst -inf min.mdinfo -ref min.rst7

cp min.ncrst ../nvt 
cd ../nvt 
    pmemd.cuda -O -i nvt.in -o nvt.out -p prmtop -c \
    min.ncrst -r nvt.ncrst -x \
        nvt.nc -inf nvt.mdinfo -ref min.ncrst

cp nvt.ncrst ../npt
cd ../npt 
    pmemd.cuda -O -i npt.in -o npt.out -p prmtop -c \
    nvt.ncrst -r npt.ncrst -x \
        npt.nc -inf npt.mdinfo -ref nvt.ncrst

if [ -d ../product ] 
then 
cp npt.ncrst ../product
cd ../product 
    pmemd.cuda -O -i pro.in -o pro.out -p prmtop -c \
    npt.ncrst -r pro.ncrst -x \
        pro.nc -inf pro.mdinfo -ref npt.ncrst \
        -x coord.nc -v vel.nc
fi 


