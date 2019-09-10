#!/bin/bash

NAME=scale-refine-mumps

export OPENBLAS_NUM_THREADS=1
ARGS="../fem_profile -q -dump_sysbin -mu 8.9e-4 -threads 16"
MUMPS="../mumps/mumps-par"
UMFPACK="../umfpack/ldivide"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(grid20 rgrid0 voronoi-s4)
R2_a=9
R2_b=18
R3_a=2
R3_b=4

# run_umfpack <test-name> <begin refine> <end refine> <-3d or none>
function run_umfpack {
    c=$1
    low=$2
    hi=$3
    dim=$4
    for r in `seq $hi -1 $low`; do
        out="$NAME/$c$dim-umfpack-r$r"
        cmd="$UMFPACK $NAME/$c$dim-r$r 16"
        mkdir -p $out/common
        $cmd > $out/common/log.txt
    done
}

# run_mumps <test-name> <begin refine> <end refine> <-3d or none>
function run_mumps {
    c=$1
    low=$2
    hi=$3
    dim=$4
    for r in `seq $hi -1 $low`; do
        out="$NAME/$c$dim-mumps-r$r"
        cmd="mpirun -n 16 $MUMPS $NAME/$c$dim-r$r 16"
        mkdir -p $out/common
        $cmd > $out/common/log.txt
    done
}

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        for r in `seq $R2_b -1 $R2_a`; do
            echo $c $r
            $ARGS -refine $r ../$c.txt -o $NAME/$c-r$r > /dev/null;
        done
    done

    for c in ${tests[@]} ; do
        for r in `seq $R3_b -1 $R3_a`; do
            $ARGS -refine $r ../$c.txt -3d -o $NAME/$c-3d-r$r > /dev/null;
        done
    done

    run_umfpack grid20 9 15 ""
    run_umfpack rgrid0 9 18 ""
    run_umfpack voronoi-s4 9 18 ""

    run_mumps grid20 9 18 ""
    run_mumps rgrid0 9 18 ""
    run_mumps voronoi-s4 9 18 ""

    run_umfpack grid20 2 3 "-3d"
    run_umfpack rgrid0 2 4 "-3d"
    run_umfpack voronoi-s4 2 4 "-3d"

    run_mumps grid20 2 3 "-3d"
    run_mumps rgrid0 2 4 "-3d"
    run_mumps voronoi-s4 2 4 "-3d"
fi


# collect <umfpack or mumps> <test-name> <begin refine> <end refine> <-3d or none>
function collect {
    solver=$1
    c=$2
    low=$3
    hi=$4
    dim=$5
    cat <<EOF > $NAME/$c$dim-$solver.txt
res solving
EOF
    for r in `seq $low $hi`; do
        res=`perl -e "{print 2*$r;}"`
        t=`grep "solve" $NAME/$c$dim-$solver-r$r/common/log.txt | sed 's/.*solve *\(.*\) ms.*/\1/g'`
        sec=`perl -e "{print $t/1000;}"`
        echo "$res $sec" >> $NAME/$c$dim-$solver.txt
    done
}

collect umfpack grid20 9 15 ""
collect umfpack rgrid0 9 18 ""
collect umfpack voronoi-s4 9 18 ""

collect mumps grid20 9 18 ""
collect mumps rgrid0 9 18 ""
collect mumps voronoi-s4 9 18 ""

collect umfpack grid20 2 3 "-3d"
collect umfpack rgrid0 2 4 "-3d"
collect umfpack voronoi-s4 2 4 "-3d"

collect mumps grid20 2 3 "-3d"
collect mumps rgrid0 2 4 "-3d"
collect mumps voronoi-s4 2 4 "-3d"

for m in mumps umfpack ; do
    sed -e "s/DDDD/18:36/g" -e "s/TTTT/$m/g" timing_vs_res_plot_other.tex > $NAME/plot-$m.tex
    sed -e "s/DDDD/4:10/g" -e "s/TTTT/$m-3d/g" timing_vs_res_plot_other.tex > $NAME/plot-$m-3d.tex
    for i in `seq 0 $((${#tests[@]}-1))` ; do
         for dim in "" "-3d" ; do
            c=${tests[$i]}$dim-$m
            rm fit.log
            gnuplot -e "fit a*x+b \"$NAME/$c.txt\" u (log10(\$1)):(log10(\$2)) via a,b" 2>/dev/null
            a=`grep Final -A 4 fit.log| grep "^a" | sed 's/a *= \([^ ]*\).*/\1/g'`
            b=`grep Final -A 4 fit.log| grep "^b" | sed 's/b *= \([^ ]*\).*/\1/g'`
            order=`perl -e "{printf('%.2f',$a);}"`
            b10=`perl -e "{print 10**$b}"`
            sed -i -e "s/LLLL$i/${tests[$i]}/g"\
                -e "s/XXXX$i/$c/g" -e "s/EEEE$i/$a/g" -e "s/CCCC$i/$b10/g" -e "s/OOOO$i/$order/g" $NAME/plot-$m$dim.tex
        done
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
