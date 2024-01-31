for i in *
do 
if [[ $i =~ [0-9] ]] && [ -d $i ] ; then
        cd $i
        gen_atom_pair_file.py ../ap.inp 
        cd ..
fi 
done 


