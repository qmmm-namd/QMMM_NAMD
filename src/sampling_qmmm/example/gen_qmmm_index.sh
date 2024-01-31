for i in *
do 
if [[ $i =~ [0-9] ]] && [ -d $i ] ; then
        cd $i
        gen_qmmm_index.py ../qi.inp
        cd ..
fi 
done 


