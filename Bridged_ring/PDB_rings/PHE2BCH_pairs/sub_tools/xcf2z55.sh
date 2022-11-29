#!/bin/bash

env=`pwd`
parts=(${env//\// })
sys=${parts[-1]}
if [ ! -d ${sys} ]; then
    mkdir ${sys}

    for lig in P01 B02
        do 
            cd $lig
            pwd
            mkdir ${lig}_res 
            prmtop=`ls|grep bonded.parm7`
            complex=`ls|grep bonded.rst7`
            pdbfile=`ls|grep bonded.pdb`
            absoluate=`pwd`
            ln -sf ${absoluate}/$prmtop ${lig}_res
            ln -sf ${absoluate}/$pdbfile ${lig}_res
            ln -sf ${absoluate}/$complex ${lig}_res 
            ln -sf ${absoluate}/press.nc ${lig}_res
            ln -sf ${absoluate}/prod.nc ${lig}_res
            mv ${lig}_res ../${sys}
            cd ..
        done
    pwd
    scp -r $sys z55:/pubhome/yjcheng02/bridge_ring_dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/mmpbsa_analysis/
    scp /pubhome/yjcheng02/pdb_rings/mmpbsa_test/finished z55:/pubhome/yjcheng02/bridge_ring_dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/mmpbsa_analysis/$sys
else
    echo "it has done already"
fi


