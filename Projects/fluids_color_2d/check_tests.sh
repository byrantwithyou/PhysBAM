#!/bin/bash

NEW_DIRECTORY="new_test_order"
if [ ! -d "$NEW_DIRECTORY" ]; then
    echo "No images to compare...terminating. Please run test-order.sh with a version of the code you know is not broken."
    exit 1
fi

OLD_DIRECTORY="old_test_order"
if [ ! -d "$OLD_DIRECTORY" ]; then
    echo "control images have not yet been copied over, doing so now. Check will be trivially successful."
    cp -r $NEW_DIRECTORY $OLD_DIRECTORY
fi

COMPARE_DIRECTORY="compare_tests"
if [ ! -d "$COMPARE_DIRECTORY" ]; then
    mkdir $COMPARE_DIRECTORY
fi

cd $NEW_DIRECTORY

bad_tests=0
for f in *.png
    do
        g=`compare -metric AE $f ../$OLD_DIRECTORY/$f ../$COMPARE_DIRECTORY/$f 2>&1`
    if [ $g -gt 0 ]; then
        echo "Mismatch in plot $f"
        let "bad_tests=$bad_tests+1"
    fi
done
if [ $bad_tests -gt 0 ]; then
   echo "$bad_tests mismatched plots found."
else 
   echo "All tests passed."
fi
