for i in *
do 
if [[ $i =~ [0-9] ]] && [ -d $i ] ; then
        cd $i
        gen_init_cond.py ../ini.inp
        cd ..
fi 
done 


