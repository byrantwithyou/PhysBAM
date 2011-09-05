#!/bin/bash

DIR=$1
FRAME=$2
CUT=$3

CWD=`pwd`
TMPDIR=`mktemp -d`
pushd $TMPDIR
    SOLID_POSITION=`rigid_body_simplicial_particles_dump_nocona -d 2 -frame ${FRAME} ${CWD}/${DIR} | tail -n 1 | sed s,^.*\ ,,g | sed s,\],,g | sed s,\-,,g`
    ${PHYSBAM}/Projects/compressible_thin_shells_2d/to_gnuplot_nocona -in ${CWD}/${DIR} -frame ${FRAME} -xcut ${CUT} -data analytic_rho | tail -n +5 > analytic_rho.data
    ${PHYSBAM}/Projects/compressible_thin_shells_2d/to_gnuplot_nocona -in ${CWD}/${DIR} -frame ${FRAME} -xcut ${CUT} -data rho_extrapolated | tail -n +5 > extrapolation_rho.data
    ${PHYSBAM}/Projects/compressible_thin_shells_2d/to_gnuplot_nocona -in ${CWD}/${DIR} -frame ${FRAME} -xcut ${CUT} -data rho_fixed | tail -n +5 > new_rho.data

cat <<EOF > file
${SOLID_POSITION}   -.2
${SOLID_POSITION}   2
EOF

#set yrange [-.5:2]
/usr/bin/gnuplot <<EOF
set yrange [-.2:2]
set key top left
set term png enhanced size 1600,800 font "Vera,12"
set output '${FRAME}.png'
plot 'analytic_rho.data' using (-\$1):2 with l t 'analytic', \
     'new_rho.data' using (-\$1):2 with l lt 3 t 'new approach', \
     'extrapolation_rho.data' using (-\$1):2 with l lt 4 t 'extrapolation', \
     'file' using 1:2 with filledcurve lt 5 t 'solid'
EOF
    mv ${FRAME}.png ${CWD}
popd
