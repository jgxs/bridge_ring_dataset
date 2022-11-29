#!/bin/bash

cd /pubhome/yjcheng02/bridge_ring_dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/lig_bch_est

ls | while read lig
    do
        cd $lig
        ls | while read pdb
            do
                cd $pdb 
                qsub /pubhome/yjcheng02/bridge_ring_dataset/Bridged_ring/PDB_rings/PHE2BCH_pairs/plop_recore/sub.qsub
                #echo test
                cd ..
            done
        cd ..
    done