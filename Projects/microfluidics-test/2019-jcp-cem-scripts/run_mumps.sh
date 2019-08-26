#!/bin/bash

NAME=mumps-cmp
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
tests=(simple grid20 rgrid0 voronoi-s4)
names=(wide grid20 rgrid0 voronoi-s4)
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
            r=${res2[$p]}
            $ARGS -o $NAME/$c-nc-t$th -refine $r -threads $th -force_blk_ref ../$c.txt;
        done

        for p in `seq 0 3` ; do
            c=${tests[$p]}
            r=${res3[$p]}
            $ARGS -3d -o $NAME/$c-3d-nc-t$th -refine $r -threads $th -force_blk_ref ../$c.txt;
        done
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

cat <<EOF > $NAME/stats.csv
Name,Resolution,Total dofs,Blocks,Tasks,Avg dofs/block
EOF
for p in `seq 0 3` ; do
    c=${tests[$p]}
    name=${names[$p]}
    r=${res2[$p]}
    tot=`grep "dofs 2d" $NAME/$c-t1/common/log.txt | sed 's/.*total: \([^<]*\).*/\1/g'`
    blks=`grep ">blocks:" $NAME/$c-t1/common/log.txt | sed 's/.*blocks: \([^<]*\).*/\1/g'`
    jbs=`grep "jobs:" $NAME/$c-t1/common/log.txt | sed 's/.*jobs: \([^<]*\).*/\1/g'`
    dofspb=`perl -e "{printf('%.1f',$tot/$blks);}"`
    rndtot=`perl -e "{printf('%.1f M',$tot/1e6);}"`
    rndjbs=`perl -e "{printf('%.1f K',$jbs/1e3);}"`
    res=`perl -e "{print 2*$r;}"`
    echo "$name,$res,$rndtot,$blks,$rndjbs,$dofspb" >> $NAME/stats.csv
done
for p in `seq 0 3` ; do
    c=${tests[$p]}
    name=${names[$p]}
    r=${res3[$p]}
    tot=`grep "dofs 3d" $NAME/$c-3d-t1/common/log.txt | sed 's/.*total: \([^<]*\).*/\1/g'`
    blks=`grep ">blocks:" $NAME/$c-3d-t1/common/log.txt | sed 's/.*blocks: \([^<]*\).*/\1/g'`
    jbs=`grep "jobs:" $NAME/$c-3d-t1/common/log.txt | sed 's/.*jobs: \([^<]*\).*/\1/g'`
    dofspb=`perl -e "{printf('%.1f',$tot/$blks);}"`
    rndtot=`perl -e "{printf('%.1f M',$tot/1e6);}"`
    rndjbs=`perl -e "{printf('%.1f K',$jbs/1e3);}"`
    res=`perl -e "{print 2*$r;}"`
    echo "$name-3d,$res,$rndtot,$blks,$rndjbs,$dofspb" >> $NAME/stats.csv
done


cat <<EOF > $NAME/result.txt
EOF

for th in 1 2 4 8 16 ; do
    for i in `seq 0 $((${#tests[@]}-1))` ; do
        c=${tests[$i]}
        name=${names[$i]}
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

            nc_comp_matrix=`grep "compute matrix" $NAME/$c$dim-nc-t$th/common/log.txt | sed 's/.*matrix *\(.*\) ms.*/\1/g'`
            nc_elim_irreg=`grep "elim irreg" $NAME/$c$dim-nc-t$th/common/log.txt | sed 's/.*irreg *\(.*\) ms.*/\1/g'`
            nc_elim_non_sep=`grep "elim non sep" $NAME/$c$dim-nc-t$th/common/log.txt | sed 's/.*sep *\(.*\) ms.*/\1/g'`
            nc_elim_3=`grep "elim 3" $NAME/$c$dim-nc-t$th/common/log.txt | sed 's/.*elim 3 *\(.*\) ms.*/\1/g'`
            nc_back=`grep "back solve" $NAME/$c$dim-nc-t$th/common/log.txt | sed 's/.*solve *\(.*\) ms.*/\1/g'`
            nc_exec_jobs=`grep "exec jobs" $NAME/$c$dim-nc-t$th/common/log.txt | sed 's/.*jobs *\(.*\) ms.*/\1/g'`
            nc_sol=`perl -e "{printf('%f',$nc_comp_matrix+$nc_elim_irreg+$nc_elim_non_sep+$nc_elim_3+$nc_back+$nc_exec_jobs);}"`

            echo $name$dim $th $sol $mumps_time $umfpack_time $nc_sol >> $NAME/result.txt
        done
    done
done


for i in `seq 0 $((${#tests[@]}-1))` ; do
    c=${tests[$i]}
    name=${names[$i]}
    for dim in "" "-3d"; do
        grep "$name$dim " $NAME/result.txt | awk '{print $2,$3}' > $NAME/$c$dim.txt
        grep "$name$dim " $NAME/result.txt | awk '{print $2,$4}' > $NAME/$c$dim-mumps.txt
        grep "$name$dim " $NAME/result.txt | awk '{print $2,$5}' > $NAME/$c$dim-umfpack.txt
        grep "$name$dim " $NAME/result.txt | awk '{print $2,$6}' > $NAME/$c$dim-nc.txt
        sed -e "s/CCCC/$c$dim/g" -e "s/MMMM/$c$dim-mumps/g" -e "s/FFFF/$c$dim-umfpack/g" -e "s/NNNN/$c$dim-nc/g" \
            -e "s/TTTT/$name$dim/g" \
            comparison_timing_plot.tex > $NAME/plot-$c$dim.tex
    done
done

cat <<EOF > $NAME/SConstruct
import os
import re
env=Environment(ENV = os.environ)
env['PSSUFFIX']=".eps"
r=re.compile(".*\.tex$")
for f in [x for x in os.listdir(".") if r.match(x)]:
    t=env.DVI(f)
    env.PostScript(t)
    env.PDF(t)
EOF
(
    cd $NAME
    scons
)

