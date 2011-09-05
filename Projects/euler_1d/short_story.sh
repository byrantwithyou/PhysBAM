#!/bin/bash

force=false

# Generate solution data, if it doesn't exist
# if $force || [[ ! -e Sod_ST/Test_1_512_semiimplicit_0_3.00 ]];   then euler_1d_nocona -sod -timesplit -last_frame 20 -resolution 512  -eno_scheme 0 -cfl 3;  fi
# if $force || [[ ! -e Sod_ST/Test_1_512_semiimplicit_0_0.90 ]];   then euler_1d_nocona -sod -timesplit -last_frame 20 -resolution 512  -eno_scheme 0 -cfl .9;  fi
# if $force || [[ ! -e Sod_ST/Test_1_512_semiimplicit_1_0.90 ]];   then euler_1d_nocona -sod -timesplit -last_frame 20 -resolution 512 -eno_scheme 1 -cfl .9;  fi
# if $force || [[ ! -e Sod_ST/Test_1__Resolution_8192_explicit ]]; then euler_1d_nocona -sod            -last_frame 20 -resolution 8192 -eno_scheme 1 -cfl .9;  fi

function render_data () {
        echo $2
	gnuplot  <<-EOF
		set term png enhanced font "Vera,12" $4
		set output '$1.png'
		$3
		plot $2
	EOF
}

function render_data_zoomed () {
        echo $2
	gnuplot  <<-EOF
		set term png enhanced font "Vera,12" $5
		set output '$1.png'
                set multiplot
                set size 1,1
                set origin 0,0
		$3
		plot $2
                set size .4,.4
                set origin 0.08,0.08
                unset key
                $4
		plot $2
	EOF
}

SEMIIMPLICIT_HIGHRES=`mktemp -q /tmp/generate_figures-XXXXX`
parsing_advection_nocona -double -start 0 -end 1 -frame 15 Sod_ST/Test_1_12800_semiimplicit_0_0.50 > ${SEMIIMPLICIT_HIGHRES}

EXPLICIT_HIGHRES=`mktemp -q /tmp/generate_figures-XXXXX`
parsing_advection_nocona -double -start 0 -end 1 -frame 15 Sod_ST/Test_1__Resolution_12800_explicit > ${EXPLICIT_HIGHRES}

parsing_advection_nocona -double -start 0 -end 1 -frame 20 Sod_ST/Test_1_12800_semiimplicit_0_0.50 > ${SEMIIMPLICIT_HIGHRES}-2
parsing_advection_nocona -double -start 0 -end 1 -frame 20 Sod_ST/Test_1__Resolution_12800_explicit > ${EXPLICIT_HIGHRES}-2

# IMPLICIT_STD=`mktemp -q /tmp/generate_figures-XXXXX`
# ${PHYSBAM}/Tools/parsing_advection/parsing_advection_nocona -start 0 -end 1 -frame 15 Sod_ST/Test_1_8192_semiimplicit_1_0.90 > ${IMPLICIT_STD}
# 
# IMPLICIT_NEW=`mktemp -q /tmp/generate_figures-XXXXX`
# ${PHYSBAM}/Tools/parsing_advection/parsing_advection_nocona -start 0 -end 1 -frame 15 Sod_ST/Test_1_512_semiimplicit_0_3.00 > ${IMPLICIT_NEW}
# 
# IMPLICIT_NEW_LOW=`mktemp -q /tmp/generate_figures-XXXXX`
# ${PHYSBAM}/Tools/parsing_advection/parsing_advection_nocona -start 0 -end 1 -frame 15 Sod_ST/Test_1_512_semiimplicit_0_0.90 > ${IMPLICIT_NEW_LOW}

CONVERGENCE_PREFIX=`mktemp -q /tmp/generate_figures-XXXXX`
for res in 100 200 400 800 1600 3200 6400; do
    parsing_advection_nocona -double -start 0 -end 1 -frame 15 Sod_ST/Test_1_${res}_semiimplicit_0_0.50 > ${CONVERGENCE_PREFIX}-${res}-0.50
    parsing_advection_nocona -double -start 0 -end 1 -frame 15 Sod_ST/Test_1_${res}_semiimplicit_0_3.00 > ${CONVERGENCE_PREFIX}-${res}-3.00
    parsing_advection_nocona -double -start 0 -end 1 -frame 15 Sod_ST/Test_1_${res}_semiimplicit_1_0.50 > ${CONVERGENCE_PREFIX}-${res}-meno
    parsing_advection_nocona -double -start 0 -end 1 -frame 15 Sod_ST/Test_1_${res}_semiimplicit_1_0.50-eno-1 > ${CONVERGENCE_PREFIX}-${res}-meno-1

    parsing_advection_nocona -double -start 0 -end 1 -frame 20 Sod_ST/Test_1_${res}_semiimplicit_0_0.50 > ${CONVERGENCE_PREFIX}-${res}-0.50-2
    parsing_advection_nocona -double -start 0 -end 1 -frame 20 Sod_ST/Test_1_${res}_semiimplicit_0_3.00 > ${CONVERGENCE_PREFIX}-${res}-3.00-2
    parsing_advection_nocona -double -start 0 -end 1 -frame 20 Sod_ST/Test_1_${res}_semiimplicit_1_0.50 > ${CONVERGENCE_PREFIX}-${res}-meno-2
    parsing_advection_nocona -double -start 0 -end 1 -frame 20 Sod_ST/Test_1_${res}_semiimplicit_1_0.50-eno-1 > ${CONVERGENCE_PREFIX}-${res}-meno-1-2
done

render_data_zoomed convergence_cons_s-L-CFL-0_50 \
         "'${CONVERGENCE_PREFIX}-100-0.50' using 1:2 with l ti 'cons. s-L, CFL .5, Res 100', \
          '${CONVERGENCE_PREFIX}-200-0.50' using 1:2 with l ti 'cons. s-L, CFL .5, Res 200', \
          '${CONVERGENCE_PREFIX}-400-0.50' using 1:2 with l ti 'cons. s-L, CFL .5, Res 400', \
          '${CONVERGENCE_PREFIX}-800-0.50' using 1:2 with l ti 'cons. s-L, CFL .5, Res 800', \
          '${CONVERGENCE_PREFIX}-1600-0.50' using 1:2 with l ti 'cons. s-L, CFL .5, Res 1600', \
          '${CONVERGENCE_PREFIX}-3200-0.50' using 1:2 with l ti 'cons. s-L, CFL .5, Res 3200', \
          '${CONVERGENCE_PREFIX}-6400-0.50' using 1:2 with l ti 'cons. s-L, CFL .5, Res 6400', \
          '${SEMIIMPLICIT_HIGHRES}' using 1:2 with l ti 'cons. s-L, CFL .5, Res 12800', \
          '${EXPLICIT_HIGHRES}' using 1:2 with l ti 'Explicit, Res 12800' lw 2" \
          "set xrange [0:1];set yrange [0:1.1];set xtics .25;set ytics .25;set obj 1 rect from 0.725,.1 to 0.775,.3 fillstyle empty;" \
          "set xrange [0.725:0.775];set yrange [.1:.3];unset xtics;unset ytics;unset obj 1" \
          "size 1200,800"

render_data_zoomed convergence_cons_s-L \
         "'${CONVERGENCE_PREFIX}-100-3.00' using 1:2 with l ti 'cons. s-L, CFL 3, Res 100', \
          '${CONVERGENCE_PREFIX}-200-3.00' using 1:2 with l ti 'cons. s-L, CFL 3, Res 200', \
          '${CONVERGENCE_PREFIX}-400-3.00' using 1:2 with l ti 'cons. s-L, CFL 3, Res 400', \
          '${CONVERGENCE_PREFIX}-800-3.00' using 1:2 with l ti 'cons. s-L, CFL 3, Res 800', \
          '${CONVERGENCE_PREFIX}-1600-3.00' using 1:2 with l ti 'cons. s-L, CFL 3, Res 1600', \
          '${CONVERGENCE_PREFIX}-3200-3.00' using 1:2 with l ti 'cons. s-L, CFL 3, Res 3200', \
          '${CONVERGENCE_PREFIX}-6400-3.00' using 1:2 with l ti 'cons. s-L, CFL 3, Res 6400', \
          '${SEMIIMPLICIT_HIGHRES}' using 1:2 with l ti 'cons. s-L, CFL .5, Res 12800', \
          '${EXPLICIT_HIGHRES}' using 1:2 with l ti 'Explicit, Res 12800' lw 2" \
          "set xrange [0:1];set yrange [0:1.1];set xtics .25;set ytics .25;set obj 1 rect from 0.725,.1 to 0.775,.3 fillstyle empty;" \
          "set xrange [0.725:0.775];set yrange [.1:.3];unset xtics;unset ytics;unset obj 1" \
          "size 1200,800"

render_data_zoomed convergence_meno-1 \
         "'${CONVERGENCE_PREFIX}-100-meno-1' using 1:2 with l ti 'MENO 1st order, Res 100', \
          '${CONVERGENCE_PREFIX}-200-meno-1' using 1:2 with l ti 'MENO 1st order, Res 200', \
          '${CONVERGENCE_PREFIX}-400-meno-1' using 1:2 with l ti 'MENO 1st order, Res 400', \
          '${CONVERGENCE_PREFIX}-800-meno-1' using 1:2 with l ti 'MENO 1st order, Res 800', \
          '${CONVERGENCE_PREFIX}-1600-meno-1' using 1:2 with l ti 'MENO 1st order, Res 1600', \
          '${CONVERGENCE_PREFIX}-3200-meno-1' using 1:2 with l ti 'MENO 1st order, Res 3200', \
          '${CONVERGENCE_PREFIX}-6400-meno-1' using 1:2 with l ti 'MENO 1st order, Res 6400', \
          '${SEMIIMPLICIT_HIGHRES}' using 1:2 with l ti 'cons. s-L, CFL .5, Res 12800', \
          '${EXPLICIT_HIGHRES}' using 1:2 with l ti 'Explicit, Res 12800' lw 2" \
          "set xrange [0:1];set yrange [0:1.1];set xtics .25;set ytics .25;set obj 1 rect from 0.725,.1 to 0.775,.3 fillstyle empty;" \
          "set xrange [0.725:0.775];set yrange [.1:.3];unset xtics;unset ytics;unset obj 1" \
          "size 1200,800"

render_data_zoomed convergence_meno \
         "'${CONVERGENCE_PREFIX}-100-meno' using 1:2 with l ti 'MENO, Res 100', \
          '${CONVERGENCE_PREFIX}-200-meno' using 1:2 with l ti 'MENO, Res 200', \
          '${CONVERGENCE_PREFIX}-400-meno' using 1:2 with l ti 'MENO, Res 400', \
          '${CONVERGENCE_PREFIX}-800-meno' using 1:2 with l ti 'MENO, Res 800', \
          '${CONVERGENCE_PREFIX}-1600-meno' using 1:2 with l ti 'MENO, Res 1600', \
          '${CONVERGENCE_PREFIX}-3200-meno' using 1:2 with l ti 'MENO, Res 3200', \
          '${CONVERGENCE_PREFIX}-6400-meno' using 1:2 with l ti 'MENO, Res 6400', \
          '${SEMIIMPLICIT_HIGHRES}' using 1:2 with l ti 'cons. s-L, CFL .5, Res 12800', \
          '${EXPLICIT_HIGHRES}' using 1:2 with l ti 'Explicit, Res 12800' lw 2" \
          "set xrange [0:1];set yrange [0:1.1];set xtics .25;set ytics .25;set obj 1 rect from 0.725,.1 to 0.775,.3 fillstyle empty;" \
          "set xrange [0.725:0.775];set yrange [.1:.3];unset xtics;unset ytics;unset obj 1" \
          "size 1200,800"

render_data_zoomed convergence_cons_s-L-CFL-0_50-2 \
         "'${CONVERGENCE_PREFIX}-100-0.50-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 100', \
          '${CONVERGENCE_PREFIX}-200-0.50-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 200', \
          '${CONVERGENCE_PREFIX}-400-0.50-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 400', \
          '${CONVERGENCE_PREFIX}-800-0.50-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 800', \
          '${CONVERGENCE_PREFIX}-1600-0.50-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 1600', \
          '${CONVERGENCE_PREFIX}-3200-0.50-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 3200', \
          '${CONVERGENCE_PREFIX}-6400-0.50-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 6400', \
          '${SEMIIMPLICIT_HIGHRES}-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 12800', \
          '${EXPLICIT_HIGHRES}-2' using 1:2 with l ti 'Explicit, Res 12800' lw 2" \
          "set xrange [0:1];set yrange [0:1.1];set xtics .25;set ytics .25;set obj 1 rect from 0.812,.1 to 0.862,.3 fillstyle empty;" \
          "set xrange [0.812:0.862];set yrange [.1:.3];unset xtics;unset ytics;unset obj 1" \
          "size 1200,800"

render_data_zoomed convergence_cons_s-L-2 \
         "'${CONVERGENCE_PREFIX}-100-3.00-2' using 1:2 with l ti 'cons. s-L, CFL 3, Res 100', \
          '${CONVERGENCE_PREFIX}-200-3.00-2' using 1:2 with l ti 'cons. s-L, CFL 3, Res 200', \
          '${CONVERGENCE_PREFIX}-400-3.00-2' using 1:2 with l ti 'cons. s-L, CFL 3, Res 400', \
          '${CONVERGENCE_PREFIX}-800-3.00-2' using 1:2 with l ti 'cons. s-L, CFL 3, Res 800', \
          '${CONVERGENCE_PREFIX}-1600-3.00-2' using 1:2 with l ti 'cons. s-L, CFL 3, Res 1600', \
          '${CONVERGENCE_PREFIX}-3200-3.00-2' using 1:2 with l ti 'cons. s-L, CFL 3, Res 3200', \
          '${CONVERGENCE_PREFIX}-6400-3.00-2' using 1:2 with l ti 'cons. s-L, CFL 3, Res 6400', \
          '${SEMIIMPLICIT_HIGHRES}-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 12800', \
          '${EXPLICIT_HIGHRES}-2' using 1:2 with l ti 'Explicit, Res 12800' lw 2" \
          "set xrange [0:1];set yrange [0:1.1];set xtics .25;set ytics .25;set obj 1 rect from 0.812,.1 to 0.862,.3 fillstyle empty;" \
          "set xrange [0.812:0.862];set yrange [.1:.3];unset xtics;unset ytics;unset obj 1" \
          "size 1200,800"

render_data_zoomed convergence_meno-1-2 \
         "'${CONVERGENCE_PREFIX}-100-meno-1-2' using 1:2 with l ti 'MENO 1st order, Res 100', \
          '${CONVERGENCE_PREFIX}-200-meno-1-2' using 1:2 with l ti 'MENO 1st order, Res 200', \
          '${CONVERGENCE_PREFIX}-400-meno-1-2' using 1:2 with l ti 'MENO 1st order, Res 400', \
          '${CONVERGENCE_PREFIX}-800-meno-1-2' using 1:2 with l ti 'MENO 1st order, Res 800', \
          '${CONVERGENCE_PREFIX}-1600-meno-1-2' using 1:2 with l ti 'MENO 1st order, Res 1600', \
          '${CONVERGENCE_PREFIX}-3200-meno-1-2' using 1:2 with l ti 'MENO 1st order, Res 3200', \
          '${CONVERGENCE_PREFIX}-6400-meno-1-2' using 1:2 with l ti 'MENO 1st order, Res 6400', \
          '${SEMIIMPLICIT_HIGHRES}-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 12800', \
          '${EXPLICIT_HIGHRES}-2' using 1:2 with l ti 'Explicit, Res 12800' lw 2" \
          "set xrange [0:1];set yrange [0:1.1];set xtics .25;set ytics .25;set obj 1 rect from 0.812,.1 to 0.862,.3 fillstyle empty;" \
          "set xrange [0.812:0.862];set yrange [.1:.3];unset xtics;unset ytics;unset obj 1" \
          "size 1200,800"

render_data_zoomed convergence_meno-2 \
         "'${CONVERGENCE_PREFIX}-100-meno-2' using 1:2 with l ti 'MENO, Res 100', \
          '${CONVERGENCE_PREFIX}-200-meno-2' using 1:2 with l ti 'MENO, Res 200', \
          '${CONVERGENCE_PREFIX}-400-meno-2' using 1:2 with l ti 'MENO, Res 400', \
          '${CONVERGENCE_PREFIX}-800-meno-2' using 1:2 with l ti 'MENO, Res 800', \
          '${CONVERGENCE_PREFIX}-1600-meno-2' using 1:2 with l ti 'MENO, Res 1600', \
          '${CONVERGENCE_PREFIX}-3200-meno-2' using 1:2 with l ti 'MENO, Res 3200', \
          '${CONVERGENCE_PREFIX}-6400-meno-2' using 1:2 with l ti 'MENO, Res 6400', \
          '${SEMIIMPLICIT_HIGHRES}-2' using 1:2 with l ti 'cons. s-L, CFL .5, Res 12800', \
          '${EXPLICIT_HIGHRES}-2' using 1:2 with l ti 'Explicit, Res 12800' lw 2" \
          "set xrange [0:1];set yrange [0:1.1];set xtics .25;set ytics .25;set obj 1 rect from 0.812,.1 to 0.862,.3 fillstyle empty;" \
          "set xrange [0.812:0.862];set yrange [.1:.3];unset xtics;unset ytics;unset obj 1" \
          "size 1200,800"



rm -rf /tmp/generate_figures-*
# mv *png ${PHYSBAM}/Papers/JCP/Fully_Conservative/figures/.
