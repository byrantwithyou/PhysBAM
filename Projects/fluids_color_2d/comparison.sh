#!/bin/bash

OLD=new_test_order
NEW=new_test_order-10
CMP=comparison

rm -rf $CMP
mkdir -p $CMP

for i in `ls $NEW` ; do
    ~/bin/src/compare $OLD/$i $NEW/$i $CMP/$i
done

wait
