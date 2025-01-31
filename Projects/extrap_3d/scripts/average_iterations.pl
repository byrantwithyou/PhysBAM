#!/usr/bin/perl -w
require strict;

my $total=0;
my $num=0;
my $max=0;

while(<>)
{
    if(/iterations.*value="([0-9]*)"/)
    {
        $total += $1;
        $num++;
        $max=$1 if $1>$max;
    }
}
printf "solves: $num    average iterations/solve: %.2f  max iterations/solve: $max\n", $total/$num;

