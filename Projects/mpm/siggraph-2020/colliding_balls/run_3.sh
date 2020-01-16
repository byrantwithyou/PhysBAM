#!/bin/bash

LF=200

Y0="9.1"
Y5="12"
SH="6"

(

    for s in $@ ; do
        o=${s/out/tex}
        i=${s/out/render}
        rm -rf $o
        mkdir $o
        ./parse_3.pl < $s/common/log.txt > $o/data.txt
        for f in `seq 0 $LF` ; do
            n=`printf %04d $f`
            head -n $(($f+2)) $o/data.txt > $o/data-$n.txt
            L=`tail -n 1 $o/data-$n.txt`
            DTV=$(perl -e "\$a=log(\$ARGV[1]+1e-30)/log(10)+$SH;if(\$a<-6){\$a=-6;}if(\$a>6){\$a=6;}\$a=\$a/5*($Y5-$Y0)+$Y0;print \$a;" $L )
            DTF=$(perl -e "\$a=log(\$ARGV[2]+1e-30)/log(10)+$SH;if(\$a<-6){\$a=-6;}if(\$a>6){\$a=6;}\$a=\$a/5*($Y5-$Y0)+$Y0;print \$a;" $L )
            DTC=$(perl -e "\$a=log(\$ARGV[3]+1e-30)/log(10)+$SH;if(\$a<-6){\$a=-6;}if(\$a>6){\$a=6;}\$a=\$a/5*($Y5-$Y0)+$Y0;print \$a;" $L )
            sed -e "s@XXXX@../$i/colliding_balls_mantra_ipr_$n.png@" -e "s@DATA@data-$n.txt@g" -e "s@DTV@$DTV@g" -e "s@DTF@$DTF@g" -e "s@DTC@$DTC@g" < template_3.tex > $o/comp-$n.tex
            ( cd $o ; pdflatex -halt-on-error comp-$n.tex ; convert comp-$n.pdf comp-$n.png ) &
        done
        wait
    done

)

#{9.1+log10(1)/5*(12.0-9.1)}
#render_ok/colliding_balls_mantra_ipr_XXXX.png
