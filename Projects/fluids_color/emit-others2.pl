#!/usr/bin/perl -w
require strict;
my $PHYSBAM=$ENV{'PHYSBAM'};
for my $s ((0,1))
{
for my $j (0..0)
{
    my $mu = 2**($j-2);
    $mu=1;
    for my $i (0..9)
    {
        my $dt = 2**($i-17);
        for my $k ((3,6))
        {
            my $rho = 2**($k-4);
            for my $r (32,64)
            {
                print `$PHYSBAM/Tools/batch/slave -p 80 -- nice ./fluids_color_2d -resolution $r -last_frame 200 -dt $dt 33 -mode $s -rho0 $rho -mu0 $mu -o shed-stab-m$s-$j-$i-$k-$r`;
                sleep(1);
            }
        }
    }
}
}


