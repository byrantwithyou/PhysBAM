#!/usr/bin/perl -w

require strict;

my $time = 0;
my $vort = 0;
my $vort_particle = 0;
my $ke_grid = 0;
my $ke_particle = 0;
my $te_grid = 0;
my $te_particle = 0;
print "time ke-grid ke-particle te-grid te-particle vort-particle\n";
while(<STDIN>)
{
    if(/substep dt: (.*)</){$time+=$1;}
    if(/l2 velocity (.*)  l2 vorticity (.*)</){$vort=$2;}
    if(/ke grid (.*)  ke particle (.*)</){$ke_grid=$1;$ke_particle=$2;}
    if(/te grid (.*)  te particle (.*)</){$te_grid=$1;$te_particle=$2;}
    if(/vort particle (.*)</){
        $vort_particle=$1;
        print "$time $ke_grid $ke_particle $te_grid $te_particle $vort_particle\n";}
}
