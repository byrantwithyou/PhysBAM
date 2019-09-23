#!/bin/bash

# run ./run_scale_refinement.sh and ./run_scale_refinement_mumps.sh first.
NAME=scale-refine-all
CELIM_NAME=scale-refine
MUMPS_NAME=scale-refine-mumps

tests=(grid20 rgrid0 voronoi-s4)

rm -rf $NAME
mkdir -p $NAME

sed -e "s/DDDD/18:36/g" timing_vs_res_plot_all.tex > $NAME/plot.tex
sed -e "s/DDDD/4:10/g" timing_vs_res_plot_all.tex > $NAME/plot-3d.tex
cp scale-refinement-all-legends.tex $NAME/legends.tex

for dim in "" "-3d" ; do
    for i in `seq 0 $((${#tests[@]}-1))` ; do
        c=${tests[$i]}$dim
        dim_digit="2"
        if [ "X$dim" = "X-3d" ] ; then
            dim_digit="3"
        fi

        cp $CELIM_NAME/$c.txt $NAME
        rm fit.log
        gnuplot -e "fit a*x+b \"$NAME/$c.txt\" u (log10(\$1)):(log10(\$3)) via a,b" 2>/dev/null
        a=`grep Final -A 4 fit.log| grep "^a" | sed 's/a *= \([^ ]*\).*/\1/g'`
        b=`grep Final -A 4 fit.log| grep "^b" | sed 's/b *= \([^ ]*\).*/\1/g'`
        order=`perl -e "{printf('%.2f',$a);}"`
        b10=`perl -e "{print 10**$b}"`
        sed -i -e "s/XXXX$i/$c/g" -e "s/XXX$i$dim_digit/\$$order\$/g" $NAME/legends.tex
        sed -i -e "s/XXXX$i/$c/g" -e "s/EEEE$i/$a/g" -e "s/CCCC$i/$b10/g" -e "s/OOOO$i/$order/g" $NAME/plot$dim.tex

        cp $MUMPS_NAME/$c-mumps.txt $NAME
        rm fit.log
        gnuplot -e "fit a*x+b \"$NAME/$c-mumps.txt\" u (log10(\$1)):(log10(\$2)) via a,b" 2>/dev/null
        a=`grep Final -A 4 fit.log| grep "^a" | sed 's/a *= \([^ ]*\).*/\1/g'`
        b=`grep Final -A 4 fit.log| grep "^b" | sed 's/b *= \([^ ]*\).*/\1/g'`
        order=`perl -e "{printf('%.2f',$a);}"`
        b10=`perl -e "{print 10**$b}"`
        sed -i -e "s/YYYY$i/$c-mumps/g" -e "s/YYY$i$dim_digit/\$$order\$/g" $NAME/legends.tex
        sed -i -e "s/YYYY$i/$c-mumps/g" -e "s/FFFF$i/$a/g" -e "s/KKKK$i/$b10/g" -e "s/PPPP$i/$order/g" $NAME/plot$dim.tex

        cp $MUMPS_NAME/$c-umfpack.txt $NAME
        rm fit.log
        gnuplot -e "fit a*x+b \"$NAME/$c-umfpack.txt\" u (log10(\$1)):(log10(\$2)) via a,b" 2>/dev/null
        a=`grep Final -A 4 fit.log| grep "^a" | sed 's/a *= \([^ ]*\).*/\1/g'`
        b=`grep Final -A 4 fit.log| grep "^b" | sed 's/b *= \([^ ]*\).*/\1/g'`
        order=`perl -e "{printf('%.2f',$a);}"`
        b10=`perl -e "{print 10**$b}"`
        sed -i -e "s/ZZZZ$i/$c-umfpack/g" -e "s/ZZZ$i$dim_digit/\$$order\$/g" $NAME/legends.tex
        sed -i -e "s/ZZZZ$i/$c-umfpack/g" -e "s/GGGG$i/$a/g" -e "s/JJJJ$i/$b10/g" -e "s/QQQQ$i/$order/g" $NAME/plot$dim.tex
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
