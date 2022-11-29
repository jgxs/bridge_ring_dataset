#!/bin/bash
while read line
do
	cd $line
	bash /pubhome/yjcheng02/pdb_rings/mmpbsa_test/find_finished.sh
	cd ..
done < menu_all

