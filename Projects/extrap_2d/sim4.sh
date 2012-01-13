#!/bin/bash

function work()
{
    d=$1
    shift
    name=`echo "$@" | sed 's/\.txt./-/g' | sed 's/txt$/png/'`
    eps-cap/compos4.pl $d/data "$@" > $d.tex
    latex $d.tex >/dev/null ; dvips $d.dvi ; convert $d.ps $d/frames/$name ; ps2pdf $d.ps
#    rm -f $d.tex $d.dvi $d.ps
}

for d in fig-mattress fig-mattress-rc2 ; do
    work $d 001.txt 080.txt 160.txt
done
