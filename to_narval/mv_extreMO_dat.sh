#!/bin/bash

structypes=('40x40' 'tempdot6' 'tempdot5')
motypes=('hi' 'lo')
outdirs=('MOs' 'energies')

for st in ${structypes[@]}; do
    for d in ${outdirs}; do
        for mt in ${motypes[@]}; do
            mv "${st}/${d}/${mt}" "${st}/${d}/${mt}100"
        done
    done
done