#!/usr/bin/perl -w

my $poles=25;
my $rsph=160;
my $dsph=160;
my $rtori=160;
my $dtori=160;
my $lathe=160;
my $list;
print <<EOF

List_Object{
    Name="Rigid_Body_List"
    Type="Rigid_Body_List"
    Prefix="input"
    Shader="BlackShader"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=true
    Flatten_Geometry=false
    Accel_By_Triangle=false
    Range=
}

EOF
;

# poles
for my $i (1..$poles){
    print <<EOF
List_Object{
    Name="Rigid$i"
    Type="Rigid_Body_Instance"
    List_Name="Rigid_Body_List"
    Shader="PoleShader"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=true
    Flatten_Geometry=false
    Accel_By_Triangle=false
    Range=$i
}

EOF
}

# 80 rigid spheres
for my $i ($poles+1..$poles+$rsph){
    print <<EOF
List_Object{
    Name="Rigid$i"
    Type="Rigid_Body_Instance"
    List_Name="Rigid_Body_List"
    Shader="SemireflectiveBlackShader3"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=false
    Flatten_Geometry=false
    Accel_By_Triangle=false
    Range=$i
}

EOF
}

# 80 rings
for my $i ($poles+$rsph+1..$poles+$rsph+$rtori){
    print <<EOF
List_Object{
    Name="Rigid$i"
    Type="Rigid_Body_Instance"
    List_Name="Rigid_Body_List"
    Shader="LatheShader"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=true
    Flatten_Geometry=false
    Accel_By_Triangle=false
    Range=$i
}

EOF
}

my @lathe_shaders=qw(LatheShaderR LatheShaderO LatheShaderY LatheShaderG LatheShaderB LatheShaderP);
#my @lathe_shaders=qw(LatheShader);
# lathe chains
for my $i ($poles+$rsph+$rtori+1..$poles+$rsph+$rtori+$lathe*6){
#    my $shader=$lathe_shaders[($i-1)%1];
    my $shader=$lathe_shaders[($i-1)%6];
    print <<EOF
List_Object{
    Name="Rigid$i"
    Type="Rigid_Body_Instance"
    List_Name="Rigid_Body_List"
    Shader="$shader"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=true
    Flatten_Geometry=false
    Accel_By_Triangle=false
    Range=$i
}

EOF
}

# deformable tori & spheres
my @ind=map {2*$_-1} (1..($dtori+$dsph));
$list=join(",",@ind);
print <<EOF
List_Object{
    Name="Deformable-tori"
    Type="Deformable_Object"
    Prefix="input"
    Shader="DefShader1"
    Smooth_Normals=true
    Subdivide_Geometry=true
    Preserve_Creases=false
    Flatten_Geometry=false
    Accel_By_Triangle=false
EOF
;
    for(@ind){
        my $color=$lathe_shaders[$_%3];
        print "    Shader$_=$color\n";}
    print <<EOF
    Range=$list
}

EOF
;
