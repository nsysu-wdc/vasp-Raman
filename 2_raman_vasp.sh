#!/bin/bash

dir="disp_POSCAR"
alist=$(cat alist)
#alist=$(sed -n '1,$p' alist)
arrlist=(${alist})

#num=${#arrlist[@]}

## step: submit
#cd ${dir}
#for i in ${arrlist[@]} ; do
#    #echo $i
#	mkdir str.${i}
#	cd str.${i}
#	ln -s ../POSCAR.${i} POSCAR
#	ln -s ../../POTCAR .
#	ln -s ../../INCAR .
#	ln -s ../../KPOINTS .
#	ln -s ../../vasp.huang.default.sh .
#	sbatch -J $i vasp.huang.default.sh
#	cd ..
#done
#cd ..

## step: summary OUTCAR
for i in ${arrlist[@]} ; do
    #echo $i
	ln -s ${dir}/str.${i}/OUTCAR OUTCAR.${i}.out
done

