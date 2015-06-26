#!/usr/bin/perl -w
require strict;
my $PHYSBAM=$ENV{'PHYSBAM'};
#for my $s ((0,1,2,3,4,5,6,8,9,10,11,12,13,14,16,17,18,19,20,21,24,28,33))
for my $s ((3,4,5,11,17,18,28))
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
                my $pr=80-$s;
                print `$PHYSBAM/Tools/batch/slave -p $pr -- nice ./fluids_color_2d -resolution $r -last_frame 200 -dt $dt $s -rho0 $rho -mu0 $mu -o shed-stab-$s-$j-$i-$k-$r`;
                sleep(1);
            }
        }
    }
}
}


