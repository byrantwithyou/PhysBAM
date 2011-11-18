#!/usr/bin/perl -w
require strict;

my $total=0;
my $num=0;

while(<>)
{
    if(/iterations.*value="([0-9]*)"/)
    {
        $total += $1;
        $num++;
    }
}
printf "solves: $num    average iteations: %.2f\n", $total/$num;

