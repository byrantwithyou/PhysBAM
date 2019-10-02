#!/usr/bin/perl -w

my $file = $ARGV[0];

open(my $F,"<",$file);
my @lines=<$F>;
close($F);

for(@lines)
{
    my $start_non_empty=/\S/;
    s/LEVESET/LEVELSET/g;

    s/GRID_2D<T>::/T_GRID::/g;
    s/typedef GRID_[1-3]D<T> T_GRID;//;
    s/typedef typename T_GRID::SCALAR T;\b/typedef typename TV::SCALAR T;/g;
    s/typedef typename T_GRID::VECTOR_T TV;//g;
    s/typedef typename T_GRID::VECTOR_INT TV_INT;/typedef VECTOR<int,TV::m> TV_INT;/g;
    s/typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;//;
    s/typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;//;
    s/typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;//;
    s/typedef typename ARRAY<T,TV_INT>::template REBIND<int>::TYPE T_ARRAYS_INT;//;
    s/typedef typename ARRAY<T,TV_INT>::template REBIND<bool>::TYPE T_ARRAYS_BOOL;//;
    s/typedef typename T_GRID::ARRAYS_SCALAR T_ARRAYS_SCALAR;//;
    s/typedef typename T_GRID::FLUID_COLLISION_BODY_LIST T_FLUID_COLLISION_BODY_LIST;//;
    s/typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;//;
    s/typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;//;
    s/TV::dimension/TV::m/g;
    s/typedef typename MESH_POLICY<TV::m-1>::MESH T_BOUNDARY_MESH;/typedef typename BASIC_SIMPLEX_POLICY<TV::m-1>::MESH T_BOUNDARY_MESH;/;
    s/typedef typename SIMPLEX_POLICY<TV,TV::m>::SIMPLEX T_SIMPLEX;/typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SIMPLEX;/;
    s/typedef typename T_GRID::FAST_LEVELSET T_FAST_LEVELSET;//;
    s/typedef typename RIGID_BODY_POLICY<TV>::RIGID_BODY T_RIGID_BODY;//;
    s/typedef typename SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;/typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;/;
    s/typedef typename SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_BOUNDARY_OBJECT;/typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_BOUNDARY_OBJECT;/;
    s/typedef typename T_GRID::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;//;
    s/typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;//;
    s/typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<ARRAY<int> >::TYPE T_FACE_ARRAYS_INT_ARRAY;//;

    s/typedef typename GRID_POLICY<TV>::UNIFORM_GRID T_GRID;//;
    s/typedef typename TV::template REBIND<int>::TYPE TV_INT;/typedef VECTOR<int,TV::m> TV_INT;/;
    s/typedef typename ARRAY<T,TV_INT>::template REBIND_LENGTH<GRID<TV>::dimension+2>::TYPE T_ARRAYS_DIMENSION_SCALAR;//;
    s/T_ARRAYS_DIMENSION_SCALAR/ARRAY<VECTOR<T,TV::m+2>,TV_INT>/g;
    
    s/typedef ARRAY<TV,TV_INT> T_ARRAYS_VECTOR;//;
    s/T_ARRAYS_VECTOR/ARRAY<TV,TV_INT>/g;
    s/typedef typename T_GRID::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;//;
    
    s/template *<class T_GRID>/template<class TV>/g;
    s/\bBOX\b/RANGE/g;
    s/\bCELL_ITERATOR \b/CELL_ITERATOR<TV> /g;
    s/\bFACE_ITERATOR \b/FACE_ITERATOR<TV> /g;
    s/\bNODE_ITERATOR \b/NODE_ITERATOR<TV> /g;
    s/\bCOLLISION_BODY<TV>/COLLISION_GEOMETRY<TV>/g;
    s/\bT_RIGID_BODY\b/RIGID_BODY<TV>/g;
    s/<T_GRID>/<TV>/g;
    s/<T_GRID,/<TV,/g;
    s/<GRID_[1-3]D<T> >/<TV>/g;
    s/(public.*)<GRID_([1-3])D<T_input> >/$1<VECTOR<T_input,$2> >/g;
    s/\bT_GRID&/GRID<TV>&/g;
    s/\bT_GRID::/GRID<TV>::/g;
    s/\bT_ARRAYS_SCALAR\b/ARRAY<T,TV_INT>/g;
    s/\bT_ARRAYS_INT\b/ARRAY<int,TV_INT>/g;
    s/\bT_ARRAYS_BOOL\b/ARRAY<bool,TV_INT>/g;
    s/\bT_FACE_ARRAYS_SCALAR\b/ARRAY<T,FACE_INDEX<TV::m> >/g;
    s/\bT_FACE_ARRAYS_INT\b/ARRAY<int,FACE_INDEX<TV::m> >/g;
    s/\bT_FACE_ARRAYS_BOOL\b/ARRAY<bool,FACE_INDEX<TV::m> >/g;
    s/\bT_FACE_ARRAYS_INT_ARRAY\b/ARRAY<ARRAY<int>,FACE_INDEX<TV::m> >/g;
    s/\bFLUID_COLLISION_BODY_LIST_UNIFORM\b/GRID_BASED_COLLISION_GEOMETRY_UNIFORM/g;
    s/\bT_FAST_LEVELSET\b/LEVELSET<TV>/g;
    s/rigid_body_particle_index/particle_index/g;
    s/\bLIST_ARRAY\b/ARRAY/g;
    s/IDENTITY_MAP/IDENTITY_ARRAY/g;
    s/RAW_ARRAY/ARRAY_VIEW/g;
    s/\bVECTOR_ND\b/ARRAY/g;
    s/\bRIGID_BODY_PARTICLE\b/RIGID_BODY_PARTICLES/g;
    s/\bRIGID_BODY_LIST\b/RIGID_BODY_COLLECTION/g;
    
    s/FILE_UTILITIES:://g;
    s/typename GRID<TV>::VECTOR_T\b/TV/g;

    s/\bSOLIDS_PARTICLE\b/GEOMETRY_PARTICLES/;

    while(s/str\(boost::format\(([^()]+?)\)%/@<Z>$1##/)
    {
        while(s/##([^@;()]+)%/##$1,/){}
        s/##/,/g;
        s/@<Z>/LOG::sprintf(/g;
    }

    while(s/str\(boost::format\(("[^"]+?")\)%/@<Z>$1##/)
    {
        while(s/##([^@;()]+)%/##$1,/){}
        s/##/,/g;
        s/@<Z>/LOG::sprintf(/g;
    }

    s/template class ([A-Za-z0-1_]+)<GRID_([123])D<(float|double)> >;/template class $1<VECTOR<$3,$2> >;/;
    s/PHYSBAM_OVERRIDE/override/g;

    s/boost::is_fundamental/is_fundamental/g;

    s/typename ARRAY<T,TV_INT>::template REBIND<([^;]+)>::TYPE/ARRAY<$1,TV_INT>/;

    s/GRID_[123]D<T> /GRID<TV> /;
    s/GRID_[123]D<T>&/GRID<TV>&/;
    s/\bIS_SAME\b/is_same/g;
    s/> TV::SCALAR /> typename TV::SCALAR /g;

    s/<GRID_2D<T>,/<TV,/g;

    s/typedef ([^;]+) \1;//g;
    s/typedef [^;]+>;//g;
    
    if(/#include/)
    {
        s/UNIFORM_GRID_ITERATOR_CELL_[123]D/CELL_ITERATOR/g;
        s/UNIFORM_GRID_ITERATOR_FACE_[123]D/FACE_ITERATOR/g;
        s/UNIFORM_GRID_ITERATOR_NODE_[123]D/NODE_ITERATOR/g;
        s/MINRES2/MINRES/g;
        s/RIGID_BODY_[123]D/RIGID_BODY/g;
        s/GRID_[123]D/GRID/g;
        s/FACE_ARRAYS_[123]D/FACE_ARRAYS/g;
        s/ARRAYS_[123]D/ARRAYS_ND/g;
        s/\bFAST_LEVELSET\b/LEVELSET/g;
        s/_PARTICLE\.h/_PARTICLES.h/g;

        s/#include.*\bNONCOPYABLE.h.*//;
        s/#include.*\bVOF_ADVECTION.h.*//;
        s/#include.*\bARRAY_DIVISION.h.*//;
        s/#include.*\bARRAY_ELEMENTWISE_MAGNITUDE_SQUARED.h.*//;
        s/#include.*\bARRAY_EXPRESSION_ENABLER.h.*//;
        s/#include.*\bARRAY_PRODUCT.h.*//;
    }

    my $end_non_empty=/\S/;
    if($end_non_empty!=$start_non_empty){$_='';}
}

open(my $G,">",$file);

for(@lines){print $G $_;}

close($G);

`$ENV{'PHYSBAM'}/Scripts/misc/fix_headers_file.sh $file`;
