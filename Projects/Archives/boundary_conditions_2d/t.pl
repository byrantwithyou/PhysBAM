#!/usr/bin/perl -w

require strict;
use Math::Trig 'acos';

my $pi=Math::Trig::acos(-1);

for my $n (1,3,5,9,13,17,23,27,35)
#for my $n (1,3,5,9,13,17)
{
    my $r=$n*4;
    my $N=$n*4;
    my $dt=$pi/4/$N;
#    print "test_bed -gibou -bc_types wwww -circle -projection -viscosity 10 -advection -steps $N -resolution $r -dt $dt -o o-c-$n > conv-$n &\n";
    print "test_bed -gibou -bc_types wwww @ARGV -steps $N -resolution $r -dt $dt -o o-c-$n > conv-$n &\n";
}
print "wait\n";


