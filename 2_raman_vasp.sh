#!/bin/bash

dir="disp_POSCAR"
alist=$(cat alist)
#alist=$(sed -n '1,$p' alist)
arrlist=(${alist})

#num=${#arrlist[@]}

## step: submit
cd ${dir}
for i in ${arrlist[@]} ; do
    #echo $i
	mkdir str.${i}
	cd str.${i}

	ln -s ../POSCAR.${i} POSCAR
	## If necessary, please use phonopy standardized POSCAR structure for calculation. ##
	#module purge
	#module load phonopy
	#cp ../POSCAR.${i} POSCAR
	#phonopy --symmetry
	#rm POSCAR BPOSCAR ; mv PPOSCAR POSCAR

	ln -s ../../POTCAR .
	ln -s ../../INCAR .
	ln -s ../../KPOINTS .

	echo " ========== " ${i}" ========== "
	## submit file ##
	ln -s ../../vasp.huang.default.sh .
	sbatch -J $i vasp.huang.default.sh

	cd ..
done
cd ..

## step: summary OUTCAR
for i in ${arrlist[@]} ; do
    #echo $i
	ln -s ${dir}/str.${i}/OUTCAR OUTCAR.${i}.out
done

