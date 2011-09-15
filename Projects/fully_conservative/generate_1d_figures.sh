#!/bin/bash

# TODO: ENO sims

force=false
if $force; then for i in . ${PHYSBAM}/Tools/parsing_advection ${PHYSBAM}/Tools/convergence_order; do
    pushd $i;
    scons -Q --random -u -j 64 CXX="/opt/icecream/bin/g++"
    popd;
done; fi

# Generate solution data, if it doesn't exist
if [[ ! -e Advection_Tests/Test_1_128_1_0.90 ]] || $force;  then fully_conservative_${PLATFORM} -scale 128;  fi
if [[ ! -e Advection_Tests/Test_1_128_1_0.90 ]] || $force;  then fully_conservative_${PLATFORM} -scale 128;  fi
if [[ ! -e Advection_Tests/Test_1_256_1_0.90 ]] || $force;  then fully_conservative_${PLATFORM} -scale 256;  fi
if [[ ! -e Advection_Tests/Test_1_512_1_0.90 ]] || $force;  then fully_conservative_${PLATFORM} -scale 512;  fi
if [[ ! -e Advection_Tests/Test_1_1024_1_0.90 ]] || $force; then fully_conservative_${PLATFORM} -scale 1024; fi
if [[ ! -e Advection_Tests/Test_1_2048_1_0.90 ]] || $force; then fully_conservative_${PLATFORM} -scale 2048; fi
if [[ ! -e Advection_Tests/Test_1_4096_1_0.90 ]] || $force; then fully_conservative_${PLATFORM} -scale 4096; fi

if [[ ! -e Advection_Tests/Test_1_128_1_0.90_conservative ]] || $force;  then fully_conservative_${PLATFORM} -conservative -scale 128;  fi
if [[ ! -e Advection_Tests/Test_1_256_1_0.90_conservative ]] || $force;  then fully_conservative_${PLATFORM} -conservative -scale 256;  fi
if [[ ! -e Advection_Tests/Test_1_512_1_0.90_conservative ]] || $force;  then fully_conservative_${PLATFORM} -conservative -scale 512;  fi
if [[ ! -e Advection_Tests/Test_1_1024_1_0.90_conservative ]] || $force; then fully_conservative_${PLATFORM} -conservative -scale 1024; fi
if [[ ! -e Advection_Tests/Test_1_2048_1_0.90_conservative ]] || $force; then fully_conservative_${PLATFORM} -conservative -scale 2048; fi
if [[ ! -e Advection_Tests/Test_1_4096_1_0.90_conservative ]] || $force; then fully_conservative_${PLATFORM} -conservative -scale 4096; fi

if [[ ! -e Advection_Tests/Test_1_128_1_0.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -conservative -scale 128;  fi
if [[ ! -e Advection_Tests/Test_1_256_1_0.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -conservative -scale 256;  fi
if [[ ! -e Advection_Tests/Test_1_512_1_0.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -conservative -scale 512;  fi
if [[ ! -e Advection_Tests/Test_1_1024_1_0.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -conservative -scale 1024; fi
if [[ ! -e Advection_Tests/Test_1_2048_1_0.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -conservative -scale 2048; fi
if [[ ! -e Advection_Tests/Test_1_4096_1_0.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -conservative -scale 4096; fi

if [[ ! -e Advection_Tests/Test_1_128_1_2.90_conservative ]] || $force;  then fully_conservative_${PLATFORM} -cfl 2.9 -conservative -scale 128;  fi
if [[ ! -e Advection_Tests/Test_1_256_1_2.90_conservative ]] || $force;  then fully_conservative_${PLATFORM} -cfl 2.9 -conservative -scale 256;  fi
if [[ ! -e Advection_Tests/Test_1_512_1_2.90_conservative ]] || $force;  then fully_conservative_${PLATFORM} -cfl 2.9 -conservative -scale 512;  fi
if [[ ! -e Advection_Tests/Test_1_1024_1_2.90_conservative ]] || $force; then fully_conservative_${PLATFORM} -cfl 2.9 -conservative -scale 1024; fi
if [[ ! -e Advection_Tests/Test_1_2048_1_2.90_conservative ]] || $force; then fully_conservative_${PLATFORM} -cfl 2.9 -conservative -scale 2048; fi
if [[ ! -e Advection_Tests/Test_1_4096_1_2.90_conservative ]] || $force; then fully_conservative_${PLATFORM} -cfl 2.9 -conservative -scale 4096; fi

if [[ ! -e Advection_Tests/Test_1_128_1_2.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -cfl 2.9 -conservative -scale 128;  fi
if [[ ! -e Advection_Tests/Test_1_256_1_2.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -cfl 2.9 -conservative -scale 256;  fi
if [[ ! -e Advection_Tests/Test_1_512_1_2.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -cfl 2.9 -conservative -scale 512;  fi
if [[ ! -e Advection_Tests/Test_1_1024_1_2.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -cfl 2.9 -conservative -scale 1024; fi
if [[ ! -e Advection_Tests/Test_1_2048_1_2.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -cfl 2.9 -conservative -scale 2048; fi
if [[ ! -e Advection_Tests/Test_1_4096_1_2.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -cfl 2.9 -conservative -scale 4096; fi

if [[ ! -e Advection_Tests/Test_2_128_1_0.90 ]] || $force;  then fully_conservative_${PLATFORM} -scale 128  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_256_1_0.90 ]] || $force;  then fully_conservative_${PLATFORM} -scale 256  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_512_1_0.90 ]] || $force;  then fully_conservative_${PLATFORM} -scale 512  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_1024_1_0.90 ]] || $force; then fully_conservative_${PLATFORM} -scale 1024 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_2048_1_0.90 ]] || $force; then fully_conservative_${PLATFORM} -scale 2048 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_4096_1_0.90 ]] || $force; then fully_conservative_${PLATFORM} -scale 4096 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_8192_1_0.90 ]] || $force; then fully_conservative_${PLATFORM} -scale 8192 -test_number 2; fi

if [[ ! -e Advection_Tests/Test_2_128_1_0.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -conservative -scale 128  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_256_1_0.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -conservative -scale 256  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_512_1_0.90_conservative_high_order ]] || $force;  then fully_conservative_${PLATFORM} -order 2 -conservative -scale 512  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_1024_1_0.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -conservative -scale 1024 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_2048_1_0.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -conservative -scale 2048 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_4096_1_0.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -conservative -scale 4096 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_8192_1_0.90_conservative_high_order ]] || $force; then fully_conservative_${PLATFORM} -order 2 -conservative -scale 8192 -test_number 2; fi

if [[ ! -e Advection_Tests/Test_2_128_1_0.90_eno ]] || $force;  then fully_conservative_${PLATFORM} -eno -scale 128  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_256_1_0.90_eno ]] || $force;  then fully_conservative_${PLATFORM} -eno -scale 256  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_512_1_0.90_eno ]] || $force;  then fully_conservative_${PLATFORM} -eno -scale 512  -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_1024_1_0.90_eno ]] || $force; then fully_conservative_${PLATFORM} -eno -scale 1024 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_2048_1_0.90_eno ]] || $force; then fully_conservative_${PLATFORM} -eno -scale 2048 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_4096_1_0.90_eno ]] || $force; then fully_conservative_${PLATFORM} -eno -scale 4096 -test_number 2; fi
if [[ ! -e Advection_Tests/Test_2_8192_1_0.90_eno ]] || $force; then fully_conservative_${PLATFORM} -eno -scale 8192 -test_number 2; fi

if [[ ! -e Velocity_Tests/Test_2_128_2 ]] || $force;                     then fully_conservative_${PLATFORM} -scale 128 -velocity -2d -test_number 2; fi
if [[ ! -e Velocity_Tests/Test_2_128_2_conservative ]] || $force;        then fully_conservative_${PLATFORM} -scale 128 -velocity -2d -test_number 2 -conservative; fi
if [[ ! -e Velocity_Tests/Test_2_128_2_conservative_energy ]] || $force; then fully_conservative_${PLATFORM} -scale 128 -velocity -2d -test_number 2 -conservative -energy; fi

if [[ ! -e Advection_Tests/Test_3_64_2_conservative ]] || $force; then fully_conservative_${PLATFORM} -scale 64 -2d -conservative -test_number 3; fi
if [[ ! -e Advection_Tests/Test_3_128_2_conservative ]] || $force; then fully_conservative_${PLATFORM} -scale 128 -2d -conservative -test_number 3; fi
if [[ ! -e Advection_Tests/Test_3_256_2_conservative ]] || $force; then fully_conservative_${PLATFORM} -scale 256 -2d -conservative -test_number 3; fi
if [[ ! -e Advection_Tests/Test_3_512_2_conservative ]] || $force; then fully_conservative_${PLATFORM} -scale 512 -2d -conservative -test_number 3; fi
if [[ ! -e Advection_Tests/Test_3_1024_2_conservative ]] || $force; then fully_conservative_${PLATFORM} -scale 1024 -2d -conservative -test_number 3; fi
if [[ ! -e Advection_Tests/Test_3_2048_2_conservative ]] || $force; then fully_conservative_${PLATFORM} -scale 2048 -2d -conservative -test_number 3; fi

function render_data () {
        echo $1
	# cat <<-EOF
	gnuplot <<-EOF
		set style function lines
		set term png enhanced font "Vera" 18 $4
		set output '$1.png'
		$3
		plot $2
	EOF
}

# Generate the composite overlay figures a la Yabe
INITIAL_DATA=`mktemp -q /tmp/generate_figures-XXXXX`
parsing_advection_${PLATFORM} -double -start 0 -end 1 -frame 0 Advection_Tests/Test_1_2048_1_0.90 > ${INITIAL_DATA}

INITIAL_DATA_2=`mktemp -q /tmp/generate_figures-XXXXX`
parsing_advection_${PLATFORM} -double -start 0 -end 1 -frame 0 Advection_Tests/Test_2_2048_1_0.90 > ${INITIAL_DATA_2}

FINAL_DATA_2_ENO=`mktemp -q /tmp/generate_figures-XXXXX`
parsing_advection_${PLATFORM} -double -start 0 -end 15 -frame 120 Advection_Tests/Test_2_8192_1_0.90_eno > ${FINAL_DATA_2_ENO}

FINAL_DATA_1_NON=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048; do
	parsing_advection_${PLATFORM} -double -start 0 -end 5 -frame 72 Advection_Tests/Test_1_${res}_1_0.90 > ${FINAL_DATA_1_NON}_${res}
        RENDER_STRING="$RENDER_STRING, '${FINAL_DATA_1_NON}_${res}' using 1:2 with l lw 1.9 ti 'Res ${res}'"
done
RENDER_STRING="'$INITIAL_DATA' using 1:2 with l lw 1.9 ti 'Initial', '$INITIAL_DATA' using (\$1+3):2 with l lw 1.9 ti 'Analytic' ${RENDER_STRING}"
render_data advection_1_non_conservative "${RENDER_STRING}" "set xrange [0:5];set xtics 1;set yrange [0:1.5];set ytics .5;set size ratio .3" "size 2000,600"

FINAL_DATA_1=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048; do
	parsing_advection_${PLATFORM} -double -start 0 -end 5 -frame 72 Advection_Tests/Test_1_${res}_1_0.90 > ${FINAL_DATA_1}_${res}
        RENDER_STRING="$RENDER_STRING, '${FINAL_DATA_1}_${res}' using 1:2 with l lw 1.9 ti 'Res ${res}'"
done
RENDER_STRING="'$INITIAL_DATA' using 1:2 with l lw 1.9 ti 'Initial', '$INITIAL_DATA' using (\$1+3):2 with l lw 1.9 ti 'Analytic' ${RENDER_STRING}"
render_data advection_1_conservative "${RENDER_STRING}"  "set xrange [0:5];set xtics 1;set yrange [0:1.5];set ytics .5;set size ratio .3" "size 2000,600"

FINAL_DATA_2_NON=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048; do
	parsing_advection_${PLATFORM} -double -start 0 -end 5 -frame 120 Advection_Tests/Test_2_${res}_1_0.90 > ${FINAL_DATA_2_NON}_${res}
        RENDER_STRING="$RENDER_STRING,  '${FINAL_DATA_2_NON}_${res}' using 1:2 with l lw 1.9 ti 'Res ${res}'"
done
RENDER_STRING="'${INITIAL_DATA_2}' using 1:2 with l lw 1.9 ti 'Initial', '${FINAL_DATA_2_ENO}' using 1:2 with l lw 1.9 ti 'ENO' ${RENDER_STRING}, sin(pi*x/5) with l lw 1.9 ti ''"
render_data advection_2_non_conservative "${RENDER_STRING}" "set xrange [0:6];set xtics 1;set yrange [0:1.5];set ytics .5;set size ratio .3" "size 2000,600"

FINAL_DATA_2=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048; do
	parsing_advection_${PLATFORM} -double -start 0 -end 5 -frame 120 Advection_Tests/Test_2_${res}_1_0.90_conservative_high_order > ${FINAL_DATA_2}_${res}
        RENDER_STRING="$RENDER_STRING,  '${FINAL_DATA_2}_${res}' using 1:2 with l lw 1.9 ti 'Res ${res}'"
done
RENDER_STRING="'${INITIAL_DATA_2}' using 1:2 with l lw 1.9 ti 'Initial', '${FINAL_DATA_2_ENO}' using 1:2 with l lw 1.9 ti 'ENO' ${RENDER_STRING}, sin(pi*x/5) with l lw 1.9 ti ''"
render_data advection_2_conservative "${RENDER_STRING}" "set xrange [0:6];set xtics 1;set yrange [0:1.5];set ytics .5;set size ratio .3" "size 2000,600"

#The four below are error plots
FINAL_DATA_ERROR_1=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048 4096; do
	convergence_order_${PLATFORM} -double -start 0 -end 15 -frame 72 Advection_Tests/Test_1_${res}_1_0.90_conservative Advection_Tests/Test_1_${res}_1_0.90_conservative > ${FINAL_DATA_ERROR_1}_${res}
        RENDER_STRING="$RENDER_STRING '${FINAL_DATA_ERROR_1}_${res}' with l lw 1.9 ti 'Res ${res}',"
done
RENDER_STRING=${RENDER_STRING:0:${#RENDER_STRING}-1}
render_data error_1_conservative "${RENDER_STRING}" "set xrange [0:5];set xtics 1;set yrange [-.5:1];set ytics .5;set size ratio .4" "size 2000,800"

FINAL_DATA_ERROR_2=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048 4096; do
	convergence_order_${PLATFORM} -double -start 0 -end 15 -frame 72 Advection_Tests/Test_1_${res}_1_0.90_conservative_high_order Advection_Tests/Test_1_${res}_1_0.90_conservative_high_order > ${FINAL_DATA_ERROR_2}_${res}
        RENDER_STRING="$RENDER_STRING '${FINAL_DATA_ERROR_2}_${res}' with l lw 1.9 ti 'Res ${res}',"
done
RENDER_STRING=${RENDER_STRING:0:${#RENDER_STRING}-1}
render_data error_1_conservative_high_order "${RENDER_STRING}" "set xrange [0:5];set xtics 1;set yrange [-.5:1];set ytics .5;set size ratio .4" "size 2000,800"

FINAL_DATA_ERROR_1_HIGH=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048 4096; do
	convergence_order_${PLATFORM} -double -start 0 -end 15 -frame 72 Advection_Tests/Test_1_${res}_1_2.90_conservative Advection_Tests/Test_1_${res}_1_2.90_conservative > ${FINAL_DATA_ERROR_1_HIGH}_${res}
        RENDER_STRING="$RENDER_STRING '${FINAL_DATA_ERROR_1_HIGH}_${res}' with l lw 1.9 ti 'Res ${res}',"
done
RENDER_STRING=${RENDER_STRING:0:${#RENDER_STRING}-1}
render_data error_1_conservative_cfl3 "${RENDER_STRING}" "set xrange [0:5];set xtics 1;set yrange [-.5:1];set ytics .5;set size ratio .4" "size 2000,800"

FINAL_DATA_ERROR_2_HIGH=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048 4096; do
	convergence_order_${PLATFORM} -double -start 0 -end 15 -frame 72 Advection_Tests/Test_1_${res}_1_2.90_conservative_high_order Advection_Tests/Test_1_${res}_1_2.90_conservative_high_order > ${FINAL_DATA_ERROR_2_HIGH}_${res}
        RENDER_STRING="$RENDER_STRING '${FINAL_DATA_ERROR_2_HIGH}_${res}' with l lw 1.9 ti 'Res ${res}',"
done
RENDER_STRING=${RENDER_STRING:0:${#RENDER_STRING}-1}
render_data error_1_conservative_cfl3_high_order "${RENDER_STRING}" "set xrange [0:5];set xtics 1;set yrange [-.5:1];set ytics .5;set size ratio .4" "size 2000,800"

FINAL_DATA_ERROR_1_LATE=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048 4096; do
	convergence_order_${PLATFORM} -double -start 0 -end 15 -frame 216 Advection_Tests/Test_1_${res}_1_2.90_conservative Advection_Tests/Test_1_${res}_1_2.90_conservative > ${FINAL_DATA_ERROR_1_LATE}_${res}
        RENDER_STRING="$RENDER_STRING '${FINAL_DATA_ERROR_1_LATE}_${res}' with l lw 1.9 ti 'Res ${res}',"
done
RENDER_STRING=${RENDER_STRING:0:${#RENDER_STRING}-1}
render_data error_1_conservative_cfl3_late "${RENDER_STRING}" "set xrange [6:11];set xtics 1;set yrange [-.5:1];set ytics .5;set size ratio .4" "size 2000,800"

FINAL_DATA_ERROR_2_LATE=`mktemp -q /tmp/generate_figures-XXXXX`
RENDER_STRING=""
for res in 128 256 512 1024 2048 4096; do
	convergence_order_${PLATFORM} -double -start 0 -end 15 -frame 216 Advection_Tests/Test_1_${res}_1_2.90_conservative_high_order Advection_Tests/Test_1_${res}_1_2.90_conservative_high_order > ${FINAL_DATA_ERROR_2_LATE}_${res}
        RENDER_STRING="$RENDER_STRING '${FINAL_DATA_ERROR_2_LATE}_${res}' with l lw 1.9 ti 'Res ${res}',"
done
RENDER_STRING=${RENDER_STRING:0:${#RENDER_STRING}-1}
render_data error_1_conservative_cfl3_late_high_order "${RENDER_STRING}" "set xrange [6:11];set xtics 1;set yrange [-.5:1];set ytics .5;set size ratio .4" "size 2000,800"


# Movie-style figures, to show how bad std s-L can get
MOVIE=`mktemp -q /tmp/generate_figures-XXXXX`
for frame in 0 24 48 72 96 120; do
	parsing_advection_${PLATFORM} -double -start 0 -end 5 -frame $frame Advection_Tests/Test_2_2048_1_0.90 > ${MOVIE}_${frame}
  	parsing_advection_${PLATFORM} -double -start 0 -end 5 -frame $frame Advection_Tests/Test_2_8192_1_0.90_eno > ${MOVIE}_${frame}_eno
	parsing_advection_${PLATFORM} -double -start 0 -end 5 -frame $frame Advection_Tests/Test_2_8192_1_0.90_conservative_high_order > ${MOVIE}_${frame}_conservative_highres
# render_data advection_movie_${frame} "'${MOVIE}_${frame}' using 1:2 with l lw 1.9 ti 'standard semi-Lagrangian', '${MOVIE}_${frame}_conservative_highres' using 1:2 with l lw 1.9 ti 'conservative semi-Lagrangian', '${MOVIE}_${frame}_eno' using 1:2 with l lw 1.9 ti 'ENO',sin(pi*x/5) ti 'velocity'" "set xrange [0:5];set yrange [0:1.5]; set key top left;set ytics .5;set size ratio .3" "size 2000,600"
    render_data advection_movie_${frame} "'${MOVIE}_${frame}' using 1:2 with l lw 1.9 ti 'standard semi-Lagrangian', '${MOVIE}_${frame}_conservative_highres' using 1:2 with l lw 1.9 ti 'conservative semi-Lagrangian', sin(pi*x/5) lw 1.9 ti 'velocity'" "set xrange [0:5];set yrange [0:1.5]; set key top left;set ytics .5;set size ratio .3" "size 2000,600"
done


MASS_HISTORY=`mktemp -q /tmp/generate_figures-XXXXX`
parsing_advection_${PLATFORM} -double -m Advection_Tests/Test_2_8192_1_0.90 > ${MASS_HISTORY}_tmp; cat -n ${MASS_HISTORY}_tmp > ${MASS_HISTORY}
parsing_advection_${PLATFORM} -double -m Advection_Tests/Test_2_8192_1_0.90_conservative_high_order > ${MASS_HISTORY}_conservative_tmp; cat -n ${MASS_HISTORY}_conservative_tmp > ${MASS_HISTORY}_conservative
render_data mass_history "'${MASS_HISTORY}' using ((\$1-1)/24):2 with l lw 1.9 ti 'standard semi-Lagrangian','${MASS_HISTORY}_conservative' using ((\$1-1)/24):2 with l lw 1.9 ti 'conservative semi-Lagrangian'" "set xrange [0:10];set xtics 2;set format x '%.0f s';set ytics .5;set key bottom left" "size 1200,800"


MASS_HISTORY_DISK=`mktemp -q /tmp/generate_figures-XXXXX`
cat -n density_disk_sl.txt > ${MASS_HISTORY_DISK}
cat -n density_disk_cons.txt > ${MASS_HISTORY_DISK}_conservative
render_data mass_history_disk "'${MASS_HISTORY_DISK}' using ((\$1-1)/24):2 with l lw 1.9 ti 'standard semi-Lagrangian','${MASS_HISTORY_DISK}_conservative' using ((\$1-1)/24):2 with l lw 1.9 ti 'conservative semi-Lagrangian'" "set xrange [0:8];set xtics 2;set format x '%.0f s';set ytics 10;set key bottom left" "size 1200,800"


OPENGL_CMD_STRING="VVS/]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]6V<F12>"

echo "streamlines"
opengl_2d_nocona -double -keys "${OPENGL_CMD_STRING}" -w 1280 -h 1280 -start_frame 10 -stop_frame 10 -so streamlines.png Velocity_Tests/Test_2_128_2
echo "streamlines_conservative"
opengl_2d_nocona -double -keys "${OPENGL_CMD_STRING}" -w 1280 -h 1280 -start_frame 10 -stop_frame 10 -so streamlines_conservative.png Velocity_Tests/Test_2_128_2_conservative
echo "streamlines_conservative_energy"
opengl_2d_nocona -double -keys "${OPENGL_CMD_STRING}" -w 1280 -h 1280 -start_frame 10 -stop_frame 10 -so streamlines_conservative_energy.png Velocity_Tests/Test_2_128_2_conservative_energy

convert streamlines.png -negate -crop 1280x1152+0+128 streamlines.png
convert streamlines_conservative.png -negate -crop 1280x1152+0+128 streamlines_conservative.png
convert streamlines_conservative_energy.png -negate -crop 1280x1152+0+128 streamlines_conservative_energy.png


function render_contours () {
        echo $1
        # cat <<-EOF
        /usr/bin/gnuplot <<-EOF
                set term table
                set output '/tmp/table.dat'
                $3
                $2
                set term png enhanced font "Vera" 18 $4
                set output '$1.png'
                plot '/tmp/table.dat' using 1:2 with lines lt -1
EOF
}

opengl_2d_nocona -double -keys "M^q" -start_frame 0 -stop_frame 0 Advection_Tests/Test_3_2048_2_0.90_conservative/; mv _matlab_dump.0 initial.dat
for i in 64 128 256 512 1024 2048; do
    opengl_2d_nocona -double -keys "M^q" -start_frame 24 -stop_frame 24 Advection_Tests/Test_3_${i}_2_0.90_conservative/; mv _matlab_dump.0 /tmp/generate_figures-contours-cons_${i}.dat
done
render_contours composite_contours \
"splot 'initial.dat' using 1:2:3,'/tmp/generate_figures-contours-cons_128.dat' using 1:2:3,'/tmp/generate_figures-contours-cons_256.dat' using 1:2:3,'/tmp/generate_figures-contours-cons_512.dat' using 1:2:3,'/tmp/generate_figures-contours-cons_1024.dat' using 1:2:3,'/tmp/generate_figures-contours-cons_2048.dat' using 1:2:3" \
"set xtics 25;set ytics 25;set nokey;set xrange [0:100];set yrange [0:100];set view 0,0;unset surface;set cntrparam levels discrete 0.5;set contour" \
"size 1200,800" 


# Kill all of that hideous white-space
for i in `ls | grep png$`; do
    convert -trim $i $i;
done
mv *png ${PHYSBAM}/Papers/JCP/Fully_Conservative/figures/
rm -rf /tmp/generate_figures-*
