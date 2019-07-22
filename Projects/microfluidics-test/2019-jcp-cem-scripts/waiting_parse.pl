#!/usr/bin/perl -w

require strict;

my $waiting = 0;
my $solving = 0;
while(<STDIN>)
{
    if(/^([^ ]*) ([^ ]*) ([^ ]*)$/)
    {
        $waiting+=$2;
        $solving+=$3;
    }
    if(/^([^ ]*)$/){
        $waiting+=$1;
    }
}
printf("%f %f %.1f\n",$waiting,$solving,$waiting/($waiting + $solving)*100);
