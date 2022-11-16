#!/bin/bash

# RLI is the short for receptor and ligand interaction.
# Here, we use the plop to compute the interaction energy of complex system.
export plop_path=/pubhome/yjcheng02/plop25.1   # need to be changed
export plop_data=/pubhome/yjcheng02/plop25.1/data
export LD_LIBRARY_PATH=$plop_path:$LD_LIBRARY_PATH
PATH=/pubhome/soft/amber/amber20_21/bin:$PATH
schrodinger_path="/pubhome/yjcheng02/schrodinger2016-2"

interpret="/pubhome/yjcheng02/iso_database/tool/I-interpret/bin/I-interpret"
interpret_para="/pubhome/yjcheng02/iso_database/tool/I-interpret/bin/parameter.txt"

# the input file contains 3 file. protein structure ends with *sin.pdb; ligand 1 ends with phe.pdb; ligand 2 end with bch.pdb
pro_file=`ls|grep sin`
phe_file=`ls|grep phe`
bch_file=`ls|grep bch`

lig_id=${phe_file:0:3}
cp "/pubhome/yjcheng02/iso_database/tool/I-interpret/bin/parameter.txt" .
${interpret} ${phe_file} ${lig_id}_phe.pdb

grep "ATOM  " ${lig_id}_phe.pdb>001.het
echo "TER" >> 001.het
sed -i "s/${lig_id}/001/g" 001.het
echo 001 > list

antechamber -i 001.het -fi pdb -o 001.mol2 -fo mol2 -c bcc -pf y -at sybyl -nc
cp 001.het 001.bak
python /pubhome/yjcheng02/interact/renumber.py --mol 001.mol2 --pdb 001.bak > 001.het

$schrodinger_path/utilities/mol2convert -imol2 001.mol2 -omae 001.mae
$schrodinger_path/utilities/hetgrp_ffgen 2005 001.mae

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
echo 001.cmxmin.con | /pubhome/yjcheng02/plop25.1/plop
echo 001.ligmin.con | /pubhome/yjcheng02/plop25.1/plop
/pubhome/yjcheng02/interact/plop_tool/extract_lig.pl out_files/001_cmxmin.out out_files/001_cmxmin.lig
echo 001.ligfix.con | /pubhome/yjcheng02/plop25.1/plop

/pubhome/yjcheng02/interact/plop_tool/nrg_set_2009 001.out rigid

