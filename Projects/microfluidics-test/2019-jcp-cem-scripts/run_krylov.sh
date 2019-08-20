#!/bin/bash

NAME=krylov-cmp
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
tests=(simple grid20 rgrid0 voronoi-s4)
res2=(16 4 6 8)
res3=(4 2 3 3)
ARGS="../fem -q -mu 8.9e-4"

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME

    for th in 1 2 4 8 16 ; do
        for p in `seq 0 3` ; do
            c=${tests[$p]}
            r=${res2[$p]}
            $ARGS -o $NAME/$c-krylov-t$th -k -kn 10 -refine $r -threads $th ../$c.txt;
            $ARGS -o $NAME/$c-t$th -refine $r -threads $th ../$c.txt;
        done

        for p in `seq 0 3` ; do
            c=${tests[$p]}
            r=${res3[$p]}
            $ARGS -3d -o $NAME/$c-3d-krylov-t$th -k -kn 10 -refine $r -threads $th ../$c.txt;
            $ARGS -3d -o $NAME/$c-3d-t$th -refine $r -threads $th ../$c.txt;
        done
    done
fi


cat <<EOF > $NAME/result.txt
EOF

for th in 1 2 4 8 16 ; do
    for c in ${tests[@]} ; do
        for dim in "" "-3d"; do
            krylov_time=`grep "krylov solve" $NAME/$c$dim-krylov-t$th/common/log.txt | sed 's/.*solve *\(.*\) ms.*/\1/g'`
            iter_time=`perl -e "{printf('%f',$krylov_time/10);}"`

            comp_matrix=`grep "compute matrix" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*matrix *\(.*\) ms.*/\1/g'`
            elim_irreg=`grep "elim irreg" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*irreg *\(.*\) ms.*/\1/g'`
            elim_non_sep=`grep "elim non sep" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*sep *\(.*\) ms.*/\1/g'`
            elim_3=`grep "elim 3" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*elim 3 *\(.*\) ms.*/\1/g'`
            back=`grep "back solve" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*solve *\(.*\) ms.*/\1/g'`
            exec_jobs=`grep "exec jobs" $NAME/$c$dim-t$th/common/log.txt | sed 's/.*jobs *\(.*\) ms.*/\1/g'`
            sol=`perl -e "{printf('%f',$comp_matrix+$elim_irreg+$elim_non_sep+$elim_3+$back+$exec_jobs);}"`

            eq_iters=`perl -e "{printf('%f',$sol/$iter_time);}"`
            echo $c$dim $th $sol $iter_time $eq_iters >> $NAME/result.txt
        done
    done
done


for c in ${tests[@]} ; do
    for dim in "" "-3d"; do
        grep "$c$dim " $NAME/result.txt | awk '{print $2,$4}' > $NAME/$c$dim-iter-time.txt
        grep "$c$dim " $NAME/result.txt | awk '{print $2,$5}' > $NAME/$c$dim-num-iters.txt
    done
done

cp krylov_iter_time.tex $NAME/plot-iter-time.tex
cp krylov_iter_time.tex $NAME/plot-iter-time-3d.tex
cp num_krylov_iter.tex $NAME/plot-num-iters.tex
cp num_krylov_iter.tex $NAME/plot-num-iters-3d.tex

for dim in "" "-3d" ; do
    for p in `seq 0 3` ; do
        c=${tests[$p]}
        sed -i -e "s/XXXX$p/$c$dim-iter-time/g" -e "s/TTTT$p/$c$dim/g" $NAME/plot-iter-time$dim.tex
        sed -i -e "s/XXXX$p/$c$dim-num-iters/g" -e "s/TTTT$p/$c$dim/g" $NAME/plot-num-iters$dim.tex
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

