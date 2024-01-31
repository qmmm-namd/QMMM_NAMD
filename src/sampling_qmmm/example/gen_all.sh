#!/usr/bin/env bash 
qmmm_path=$JADE_SAMPLE_QMMM
top_file=top.inp
group_file=group.inp
ap_file=ap.inp
qi_file=qi.inp
ini_file=ini.inp

if [ -f $top_file ] ; then top_file=$PWD/$top_file; else top_file=$qmmm_path/example/top.inp ;fi
if [ -f $group_file ] ; then group_file=$PWD/$group_file; else group_file=$qmmm_path/example/group.inp ; fi
if [ -f $ap_file ] ; then ap_file=$PWD/$ap_file; else ap_file=$qmmm_path/example/ap.inp ; fi
if [ -f $qi_file ] ; then qi_file=$PWD/$qi_file; else qi_file=$qmmm_path/example/qi.inp ; fi
if [ -f $ini_file ] ; then ini_file=$PWD/$ini_file; else ini_file=$qmmm_path/example/ini.inp ; fi


for i in *
do 
if [[ $i =~ [0-9] ]] && [ -d $i ] ; then
cd $i
# the first step should be the generation of new top files, the order of other processes can be random
$qmmm_path/gen_top_file.py $top_file

# generate group file 
$qmmm_path/gen_group_file.py $group_file
# generate atom_pair file
$qmmm_path/gen_atom_pair_file.py $ap_file
# generate qmmm_index file
$qmmm_path/gen_qmmm_index.py $qi_file
# generate stru_xyz.in and vel_xyz.in
$qmmm_path/gen_init_cond.py $ini_file

# move top files to the store directory 
mkdir QMMM_EXAM; mv *.top group.json QMMM_EXAM
cd ..
fi
done 

echo 'top input file :' $top_file
echo 'group input file :' $group_file
echo 'atom_pair input file :' $ap_file 
echo 'qmmm_index input file :' $qi_file
echo 'init condision input file :' $ini_file

