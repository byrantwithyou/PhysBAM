#!/usr/bin/perl -w
require strict;
my $PHYSBAM=$ENV{'PHYSBAM'};
for my $j (0..4)
{
    my $mu = 2**($j-2);
    my $mu1 = $mu * 3;
    for my $i (0..9)
    {
        my $dt = 2**($i-13);
         for my $k (0..9)
        {
            my $rho = 2**($k-4);
            my $rho1 = $rho*2;
            for my $r (32,64)
            {
                print `$PHYSBAM/Tools/batch/slave -p 80 -- nice ./fluids_color_2d -resolution $r -last_frame 1000 -dt $dt 104 -rho0 $rho -rho1 $rho1 -mu0 $mu -mu1 $mu1 -o stab-at-ii-$j-$i-$k-$r`;
                sleep(1);
            }
        }
    }
}


