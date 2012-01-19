#!/usr/bin/perl -wnpi
my $v='[a-zA-Z_][.a-zA-Z0-9_]*';

s/for\(ID ($v)\(1\); *\1<=([-*<>:+a-z\/A-Z0-9_.()]+); *\1\+\+\)/for(ID $1(0);$1<$2;$1++)/g;


