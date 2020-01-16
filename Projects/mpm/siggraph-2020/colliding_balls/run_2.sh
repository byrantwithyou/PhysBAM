#!/bin/bash

# Usage: ./run_2.sh out_ok
# Usage: ./run_2.sh out_ok out_pic out_no_c
# Usage: ./run_2.sh out_*

LF=200

(

for s in $@ ; do
    o=${s/out/render}
    echo RENDERING $s
    rm -rf render
    for i in {0..9} ; do
        rm -f sim
        ln -s $s sim
        gnome-terminal --disable-factory -e "bash -c 'pushd /opt/hfs18.0 >&/dev/null ; source houdini_setup -q ; popd >&/dev/null ; hrender -j 8 -d mantra_ipr colliding_balls.hipnc -e -f $i $LF -i 10'" &
    done
    wait
    echo FINISHED $s MOVING TO $o
    rm -rf $o
    mv render $o
    for i in $o/*exr ; do convert $i -gamma 2.2 -background '#cccccc' -flatten `echo $i | sed 's/\./_/g' | sed 's/_exr$/.png/'` & done
done

)

