#!/bin/bash

rm -f conv-* cc.txt
rm -rf o-c-*
t.pl "$@" | bash
grep -h END conv-? conv-?? conv-??? 2>/dev/null | sed 's/W@Z-//' | sed 's/[()]//g' > cc.txt
order-of-accuracy 2 3 < cc.txt
