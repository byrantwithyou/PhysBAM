#!/usr/bin/perl -w

my @heights=(1,2,3,4,5,4,3,3,4,5,6,2,1);

for $i (0..$#heights){
    for $j (0..$heights[$i]){
        my $y=.2*$j-.1;
        my $z=-6-2*$i;
        my $xmin=-10;
        my $xmax=10;
        my $ymin=$y+.25;
        my $ymax=$y+.05;
        my $zmin=$z-1;
        my $zmax=$z+1;
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

