#!/bin/bash
num=0
>log
find -name "prod.info"|while read line
do
    info=`tail -n 2 $line|head -n 1`
    if [[ $info =~ "seconds" ]]
    then
        echo "\n" >> log
    fi
done
lines=`wc -l log|awk '{print $1}'`
if [[ $lines == '2' ]]
    then 
    echo "all finish"
    bash /pubhome/yjcheng02/pdb_rings/mmpbsa_test/xcf2z55.sh
else 
    echo "not finish"
fi
