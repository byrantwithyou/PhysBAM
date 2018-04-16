#!/bin/bash

REF=scans/M1-00819_FemaleRight.obj
PTS=${REF/.obj/-pts.txt}

LS_RES=20
LS_RAD=2
REG=1

#./align_bone_start -i $REF -right -`grep front $PTS` -`grep in $PTS` -`grep out $PTS` -levelset_radius $LS_RAD -levelset_res $LS_RES -ls reference_levelset.phi.gz -reg $REG -compute_levelset

echo A

LS_RAD=1.5

for S in scans/*.obj ; do
    PS=${S/.obj/-pts.txt}

    SIDE=-left
    echo $PS | grep -iq right && SIDE=-right
    echo start $S
    ./align_bone_start -i $S $SIDE -`grep front $PS` -`grep in $PS` -`grep out $PS` -levelset_radius $LS_RAD -levelset_res $LS_RES -ls reference_levelset.phi.gz -reg $REG
    echo end $S
    OUT=`echo $S | sed 's@.*/\(.*\).obj@\1@'`
    rm -rf out-$OUT
    mv output out-$OUT

done


echo B




#
#-compute_levelset  compute level set
#-front             index of vertex on front flat part (-1)
#-i                 file containing bone geometry ()
#-in                index of vertex on inside top round part (-1)
#-left              is left leg
#-levelset_radius   scale level set radius compared with fit (2)
#-levelset_res      level set resolution (10)
#-ls                file containing template level set ()
#-out               index of vertex on outside top round part (-1)
#-reg               regularization parameter (1)
#-right             is left leg
#

