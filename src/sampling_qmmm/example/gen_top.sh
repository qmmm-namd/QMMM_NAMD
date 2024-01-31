for i in *
do 
if [[ $i =~ [0-9] ]] && [ -d $i ] ; then
        cd $i
        gen_top_file.py ../top.inp 
        cd ..
fi 
done 


