#!/usr/bin/perl -wnpi
my $v='[a-zA-Z_][a-zA-Z0-9_]*';

s/for\(($v)=1;\1<=([-*:+a-zA-Z0-9_.()]+);\1\+\+\)/for($1=0;$1<$2;$1++)/g;


