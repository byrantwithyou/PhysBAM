#!/bin/bash

# Assumed type: f=float, d=double
# Usage: ./insert_float_bit.sh f file
# Usage: ./insert_float_bit.sh d file

FL=$1
shift
for IN in "$@" ; do
    if [ "$FL" = "f" ] ; then
        SUF="float"
        CH='\x00'
    else
        SUF="double"
        CH='\x01'
    fi
    mv "$IN" "$IN.$SUF"
    if `echo $IN | grep -q '\.gz$'` ; then
        (
            echo -n -e "$CH"
            gzip -d -c "$IN.$SUF"
        ) | gzip > "$IN"
    else
        (
            echo -n -e "$CH"
            cat "$IN.$SUF"
        ) > "$IN"
    fi
done
