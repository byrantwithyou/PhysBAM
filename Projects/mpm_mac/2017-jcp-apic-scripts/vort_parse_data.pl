#!/usr/bin/perl -w

require strict;

my $time = 0;
my $vel = 0;
my $vort = 0;
my $tay = 0;
my $init = 1;
my $init_vel = 0;
my $init_vort = 0;
my $init_tay = 0;
print "time vel vort tay\n";
while(<STDIN>)
{
    if(/substep dt: (.*)</){$time+=$1;}
    if(/l2 velocity (.*)  l2 vorticity (.*)</){$vel=$1;$vort=$2;}
    if(/taylor total (.*)</)
    {
        $tay=$1;
        if($init)
        {
            $init_vel = $vel;
            $init_vort = $vort;
            $init_tay = $tay;
            $init=0;
        }
        else
        {
            if($init_vel!=0 && $init_vort!=0 && $init_tay!=0)
            {
                $vel/=$init_vel;
                $vort/=$init_vort;
                $tay/=$init_tay;
            }
            print "$time $vel $vort $tay\n";
        }
    }
}
