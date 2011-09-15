#!/usr/bin/perl -w

my $poles=25;
my $rsph=80;
my $dsph=80;
my $rtori=80;
my $dtori=80;
my $lathe=80;

# poles
my $list=join(',',(1..$poles));
print <<EOF
List_Object{
    Name="Rigid-poles"
    Type="Rigid_Body_List"
    Prefix="input"
    Shader="RedShader"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=true
    Range=$list
}

EOF
;

# 80 rigid spheres
$list=join(',',($poles+1..$poles+$rsph));
print <<EOF
List_Object{
    Name="Rigid-spheres"
    Type="Rigid_Body_List"
    Prefix="input"
    Shader="BlueShader"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=false
    Range=$list
}

EOF
;

# 80 rings
$list=join(',',($poles+$rsph+1..$poles+$rsph+$rtori));
print <<EOF
List_Object{
    Name="Rigid-rings"
    Type="Rigid_Body_List"
    Prefix="input"
    Shader="YellowShader"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=true
    Range=$list
}

EOF
;

# lathe chains
$list=join(',',($poles+$rsph+$rtori+1..$poles+$rsph+$rtori+$lathe*6));
print <<EOF
List_Object{
    Name="Rigid-lathe"
    Type="Rigid_Body_List"
    Prefix="input"
    Shader="LatheShaderG"
    Smooth_Normals=true
    Subdivide_Geometry=false
    Preserve_Creases=true
    Range=$list
}

EOF
;

# deformable tori
for my $i (1..$dtori){
    print <<EOF
List_Object{
    Name="Deformable$i"
    Type="Deformable_Object"
    Prefix="input"
    Shader="LatheShaderO"
    Smooth_Normals=true
    Preserve_Creases=false
    Range=$i
}

EOF
}

# deformable spheres
for my $i ($dtori+1..$dtori+$dsph){
    print <<EOF
List_Object{
    Name="Deformable$i"
    Type="Deformable_Object"
    Prefix="input"
    Shader="LatheShaderP"
    Smooth_Normals=true
    Preserve_Creases=false
    Range=$i
}

EOF
}


