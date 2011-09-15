#!/bin/bash

force=false

# Generate solution data, if it doesn't exist
if $force || [[ ! -e Sod_ST/Test_1_512_semiimplicit_0_3.00 ]];   then euler_1d_nocona -sod -timesplit -last_frame 20 -resolution 512  -eno_scheme 0 -cfl 3;  fi
if $force || [[ ! -e Sod_ST/Test_1_512_semiimplicit_0_0.90 ]];   then euler_1d_nocona -sod -timesplit -last_frame 20 -resolution 512  -eno_scheme 0 -cfl .9;  fi
if $force || [[ ! -e Sod_ST/Test_1_512_semiimplicit_1_0.90 ]];   then euler_1d_nocona -sod -timesplit -last_frame 20 -resolution 512 -eno_scheme 1 -cfl .9;  fi
if $force || [[ ! -e Sod_ST/Test_1__Resolution_8192_explicit ]]; then euler_1d_nocona -sod            -last_frame 20 -resolution 8192 -eno_scheme 1 -cfl .9;  fi

function render_data () {
        echo $2
	gnuplot <<-EOF
		set term png enhanced font "Vera,24" $4
		set output '$1.png'
		$3
		plot $2
	EOF
}

EXPLICIT_HIGHRES=`mktemp -q /tmp/generate_figures-XXXXX`
${PHYSBAM}/Tools/parsing_advection/parsing_advection_nocona -start 0 -end 1 -frame 15 Sod_ST/Test_1__Resolution_8192_explicit > ${EXPLICIT_HIGHRES}

IMPLICIT_STD=`mktemp -q /tmp/generate_figures-XXXXX`
${PHYSBAM}/Tools/parsing_advection/parsing_advection_nocona -start 0 -end 1 -frame 15 Sod_ST/Test_1_8192_semiimplicit_1_0.90 > ${IMPLICIT_STD}

IMPLICIT_NEW=`mktemp -q /tmp/generate_figures-XXXXX`
${PHYSBAM}/Tools/parsing_advection/parsing_advection_nocona -start 0 -end 1 -frame 15 Sod_ST/Test_1_512_semiimplicit_0_3.00 > ${IMPLICIT_NEW}

IMPLICIT_NEW_LOW=`mktemp -q /tmp/generate_figures-XXXXX`
${PHYSBAM}/Tools/parsing_advection/parsing_advection_nocona -start 0 -end 1 -frame 15 Sod_ST/Test_1_512_semiimplicit_0_0.90 > ${IMPLICIT_NEW_LOW}

#         "'${EXPLICIT_HIGHRES}' using 1:2 with l ti 'explicit',
render_data sod_density \
         "'${IMPLICIT_STD}' using 1:2 with l ti 'MENO advection, CFL .9', \
          '${IMPLICIT_NEW_LOW}' using 1:2 with l ti 'cons. s-L advection, CFL .9', \
          '${IMPLICIT_NEW}' using 1:2 with l ti 'cons. s-L advection, CFL 3' lt 3" "set yrange [0:1.1];set xtics .5;set xtics .5;" "size 1200,800"

rm -rf /tmp/generate_figures-*
mv *png ${PHYSBAM}/Papers/JCP/Fully_Conservative/figures/.
