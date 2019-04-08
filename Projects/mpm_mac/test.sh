#!/bin/bash

rm -rf  with_test no_test 

ARGS="./mpm_mac 16 -bc_periodic -resolution 32 -last_frame 10 -test_output_prefix test -regular_seeding"
$ARGS -o no_test  > n
$ARGS -test_periodic -o with_test  > w

octave -q --eval 'load "no_test/test-001.txt" ; XX=X;VV=V; load "with_test/test-001.txt" ; max(abs(X-XX)) , max(abs(V-VV))'
meld n w
