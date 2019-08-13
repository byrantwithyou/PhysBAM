#!/bin/bash

NAME=mumps-cmp
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
tests=(simple grid20 rgrid0 voronoi-s4)
res2=(16 4 6 8)
res3=(4 2 3 3)
ARGS="../fem -q -mu 8.9e-4"
MUMPS="../mumps/mumps-par"
UMFPACK="../umfpack/ldivide"

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME

    dump="-dump_sysbin"
    for th in 1 2 4 8 16 ; do
        for p in `seq 0 3` ; do
            c=${tests[$p]}
            r=${res2[$p]}
            $ARGS -o $NAME/$c-t$th -refine $r -threads $th $dump ../$c.txt;
        done

        for p in `seq 0 3` ; do
            c=${tests[$p]}
            r=${res3[$p]}
            $ARGS -3d -o $NAME/$c-3d-t$th -refine $r -threads $th $dump ../$c.txt;
        done
        dump=""
    done

    for th in 1 2 4 8 16 ; do
        for p in `seq 0 3` ; do
            c=${tests[$p]}
            mkdir -p $NAME/$c-mumps-t$th/common
            mpirun -n $th $MUMPS $NAME/$c-t1 $th > $NAME/$c-mumps-t$th/common/log.txt
        done

        for p in `seq 0 3` ; do
            c=${tests[$p]}
            mkdir -p $NAME/$c-3d-mumps-t$th/common
            mpirun -n $th $MUMPS $NAME/$c-3d-t1 $th > $NAME/$c-3d-mumps-t$th/common/log.txt
        done
    done

    for th in 1 2 4 8 16 ; do
        for p in `seq 0 3` ; do
            c=${tests[$p]}
            mkdir -p $NAME/$c-umfpack-t$th/common
            $UMFPACK $NAME/$c-t1 $th > $NAME/$c-umfpack-t$th/common/log.txt
        done

        for p in `seq 0 3` ; do
            c=${tests[$p]}
            mkdir -p $NAME/$c-3d-umfpack-t$th/common
            $UMFPACK $NAME/$c-3d-t1 $th > $NAME/$c-3d-umfpack-t$th/common/log.txt
        done
    done
fi


cat <<EOF > $NAME/result.txt
EOF

for th in 1 2 4 8 16 ; do
    for c in ${tests[@]} ; do
        for dim in "" "-3d"; do
            mumps_time=`grep "solve" $NAME/$c$dim-mumps-t$th/common/log.txt | sed 's/.*solve *\(.*\) ms.*/\1/g'`
            umfpack_time=`grep "solve" $NAME/$c$dim-umfpack-t$th/common/log.txt | sed 's/.*solve *\(.*\) ms.*/\1/g'`

            comp_matrix=`grep "compute matrix" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*matrix *\(.*\) ms.*/\1/g'`
            elim_irreg=`grep "elim irreg" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*irreg *\(.*\) ms.*/\1/g'`
            elim_non_sep=`grep "elim non sep" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*sep *\(.*\) ms.*/\1/g'`
            elim_3=`grep "elim 3" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*elim 3 *\(.*\) ms.*/\1/g'`
            back=`grep "back solve" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*solve *\(.*\) ms.*/\1/g'`
            exec_jobs=`grep "exec jobs" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*jobs *\(.*\) ms.*/\1/g'`
            sol=`perl -e "{printf('%f',$comp_matrix+$elim_irreg+$elim_non_sep+$elim_3+$back+$exec_jobs);}"`

            echo $c$dim $th $sol $mumps_time $umfpack_time >> $NAME/result.txt
        done
    done
done
