#!/bin/bash

structypes=('300' 'q400' '500')
nstrucs=('216' '263' '300')
motypes=('lo' 'virtual_w_HOMO' 'hi')


for k in {0..2} ; do
	st="${structypes[k]}"
	nn="${nstrucs[k]}"
	dn="sAMC-$st"
	mkdir -p "$dn/hamiltonians/hvals"
	mkdir -p "$dn/hamiltonians/inds"
	for mt in ${motypes[@]}; do
		mkdir -p "$dn/MOs/$mt"
		mkdir -p "$dn/energies/$mt"
	done
	cd $dn 
	sed -e "s/##NNN##/${nn}/g" -e "s/##TTT##/${st}/g" ../submit.template > submit_array.sh
	ln -s ../*.jl .
	ln -s ../requirements.txt .
	sbatch submit_array.sh
	cd -
done


