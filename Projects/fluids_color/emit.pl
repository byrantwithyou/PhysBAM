#!/usr/bin/perl -w
require strict;
my $PHYSBAM=$ENV{'PHYSBAM'};
for my $j (0..4)
{
    my $mu = 2**($j-2);
    my $mu1 = $mu * 3;
    for my $i (0..9)
    {
        my $dt = 2**($i-16);
        for my $k (0..9)
        {
            my $rho = 2**($k-4);
            my $rho1=$rho*2;
            print `$PHYSBAM/Tools/batch/slave -p 80 -- nice ./fluids_color -use_ls -resolution 64 -last_frame 1000 -dt $dt 108 -o relax-stab-$j-$k-$i-hr -rho0 $rho -rho1 $rho1 -mu0 $mu -mu1 $mu1`;
            sleep(1);
        }
    }
}

