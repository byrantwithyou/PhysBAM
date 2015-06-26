#!/usr/bin/perl -w
require strict;
my $PHYSBAM=$ENV{'PHYSBAM'};
#for my $j (0..4)
for my $j (2)
{
    my $mu = 2**($j-2);
    for my $i (0..9)
    {
        my $dt = 2**($i-13);
#        for my $k (0..9)
        for my $k (3..6)
        {
            my $rho = 2**($k-4);
#            for my $r (32,64)
            for my $r (128)
            {
                print `$PHYSBAM/Tools/batch/slave -p 80 -- nice ./fluids_color -resolution $r -last_frame 200 -dt $dt 112 -rho0 $rho -mu0 $mu -o shed-stab-$j-$i-$k-$r`;
                sleep(1);
            }
        }
    }
}


