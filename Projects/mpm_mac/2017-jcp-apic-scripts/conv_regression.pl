#!/usr/bin/perl -w
require strict;

my $in=$ARGV[0];
my $out=$ARGV[1];
open(my $I,"<",$in);
open(my $O,">",$out);
my @r=();
my $head=<$I>;
print $O $head;

my $sx=0;
my $sy=0;
my $sz=0;
my $sx2=0;
my $sy2=0;
my $sz2=0;
my $sxy=0;
my $sxz=0;

while(<$I>)
{
    /(.*) (.*) (.*)/;
    my ($x,$y,$z)=(log($1),log($2),log($3));
    push @r, $1;
    $sx+=$x;
    $sy+=$y;
    $sz+=$z;
    $sx2+=$x*$x;
    $sy2+=$y*$y;
    $sz2+=$z*$z;
    $sxy+=$x*$y;
    $sxz+=$x*$z;
}
close($I);
my $n=@r;
my $ay = ($sxy*$n-$sx*$sy)/($sx2*$n-$sx*$sx);
my $by = ($sy*$sx2-$sxy*$sx)/($sx2*$n-$sx*$sx);
my $az = ($sxz*$n-$sx*$sz)/($sx2*$n-$sx*$sx);
my $bz = ($sz*$sx2-$sxz*$sx)/($sx2*$n-$sx*$sx);
for(@r)
{
    my $y=exp($ay*log($_)+$by);
    my $z=exp($az*log($_)+$bz);
    print $O "$_ $y $z\n";
}
close($O);
printf("%.2f %.2f\n",-$ay,-$az);
