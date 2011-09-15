#!/usr/bin/perl -w

my @heights=(1,2,3,4,5,4,3,3,4,5,6,2,1,-1);

for $j (0..6){
    my $y=.2*$j-.1;
    my $xmin=-10;
    my $xmax=10;
    my $ymin=$y+.05;
    my $ymax=$y+.25;
    my $zmin='';
    my $zmax='';
    my $cont=0;
    for $i (0..$#heights){
        if($heights[$i]>=$j){
            my $z=-8-2*$i;
            if($i=0 || $heights[$i-1]<$j){$zmax=$z+1;}
            $zmin=$z-1;}
        elsif($i>0 && $heights[$i-1]>=$j){
            print <<EOF
Object{
    Name="Steps$i$j"
    Type="Box"
    Xmin=$xmin
    Xmax=$xmax
    Ymin=$ymin
    Ymax=$ymax
    Zmin=$zmin
    Zmax=$zmax
    Shader="ShaderHeight$j"
}
EOF
        }
    }
}

