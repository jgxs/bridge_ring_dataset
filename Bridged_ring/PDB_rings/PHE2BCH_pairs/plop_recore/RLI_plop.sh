#!/bin/bash

# RLI is the short for receptor and ligand interaction.
# Here, we use the plop to compute the interaction energy of complex system.
export plop_path=/pubhome/yjcheng02/plop25.1   # need to be changed
export plop_data=/pubhome/yjcheng02/plop25.1/data
export LD_LIBRARY_PATH=$plop_path:$LD_LIBRARY_PATH
PATH=/pubhome/soft/amber/amber20_21/bin:$PATH
schrodinger_path="/pubhome/yjcheng02/schrodinger2016-2"

# the input file contains 3 file. protein structure ends with *sin.pdb; ligand 1 ends with phe.pdb; ligand 2 end with bch.pdb
pro_file=`ls|grep sin`
phe_file=`ls|grep phe`
bch_file=`ls|grep bch`

lig_id=${phe_file:0:3}
sed -i "s/UNL/${lig_id}/g" $bch_file
line=`head -n 1 $phe_file`
sed -i "s/${line:20:3}/   /g" $phe_file

sed "s/${lig_id}/001/g" $phe_file | grep ^HETATM > 001.het
echo TER >> 001.het
echo 001 > list

sed "s/${lig_id}/002/g" $bch_file |grep ^HETATM > 002.het
echo TER >> 002.het
echo 002 >> list

charge=`grep REMARK $phe_file|awk '{print $2}'`
for lig in 001 002
    do
        {
            antechamber -i ${lig}.het -fi pdb -o ${lig}.mol2 -fo mol2 -c bcc -pf y -at sybyl -nc $charge
        } ||{
            antechamber -i ${lig}.het -fi pdb -o ${lig}.mol2 -fo mol2 -c bcc -pf y -at sybyl -nc $charge -m 2
        }
        cat ${lig}.het |grep HETATM > ${lig}.bak
        # python /pubhome/yjcheng02/interact/renumber.py --mol ${lig}.mol2 --pdb ${lig}.bak 
        $schrodinger_path/utilities/mol2convert -imol2 ${lig}.mol2 -omae ${lig}.mae
        $schrodinger_path/utilities/hetgrp_ffgen 2005 ${lig}.mae
    done 

echo """rand seed 1111111111
logfile rec_h_opt.log
datadir $plop_data
load pdb ${pro_file} het no opt yes
write pdb rec_h_opt.pdb
"""> rec.con
echo -e rec.con| $plop_path/plop

/pubhome/yjcheng02/interact/plop_tool/nwdock_pdbconstruct rec_h_opt.pdb list
/pubhome/yjcheng02/interact/plop_tool/nwdock_confile_2009 false . rigid
sed -i "s?plop/6.0?plop/25.1?" *.con
mkdir log_files
mkdir out_files

for lig in 001 002
do
echo ${lig}.cmxmin.con | /pubhome/yjcheng02/plop25.1/plop32-static
echo ${lig}.ligmin.con | /pubhome/yjcheng02/plop25.1/plop32-static
/pubhome/yjcheng02/interact/plop_tool/extract_lig.pl out_files/001_cmxmin.out out_files/${lig}_cmxmin.lig
echo ${lig}.ligfix.con | /pubhome/yjcheng02/plop25.1/plop32-static
done 
/pubhome/yjcheng02/interact/plop_tool/nrg_set_2009 res.out rigid

