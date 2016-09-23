//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

using namespace PhysBAM;

typedef float RW;
typedef float T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;
typedef VECTOR<int,2> IV2;

// r = 1/4
// R = 1/2
// z axis

T r = .5;
T R = 1;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();

    TETRAHEDRALIZED_VOLUME<T> tet;
    FILE_UTILITIES::Read_From_File(STREAM_TYPE(1.0f),"/home/craig/PhysBAM/Public_Data/Tetrahedralized_Volumes/adaptive_torus_float.tet.gz",tet);
    tet.Update_Number_Nodes();
    tet.Initialize_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>* ts=tet.triangulated_surface->Create_Compact_Copy();
    ts->Loop_Subdivide();

    int n=ts->particles.X.m;
    ARRAY<VECTOR<T,2> > pts(n*4);
    ARRAY<TV_INT> elem=ts->mesh.elements;
    for(int i=0;i<n;i++){
        TV X=ts->particles.X(i);
        T theta=atan2(X.y,X.x);
        T phi=atan2(X.z,sqrt(sqr(X.x)+sqr(X.y))-R);
        VECTOR<T,2> a(theta,phi),b(theta+2*pi,phi),c(theta,phi+2*pi),d(theta+2*pi,phi+2*pi);
        pts(i)=a;
        pts(i+n)=b;
        pts(i+2*n)=c;
        pts(i+3*n)=d;}

    LOG::printf("range %P\n",RANGE<VECTOR<T,2> >::Bounding_Box(pts));

    for(int i=0;i<elem.m;i++){
        VECTOR<VECTOR<T,2>,3> uv(pts.Subset(elem(i)));
        RANGE<VECTOR<T,2> > bb=RANGE<VECTOR<T,2> >::Bounding_Box(uv);
        VECTOR<T,2> range=bb.Edge_Lengths();
        for(int a=0;a<2;a++)
            if(range(a)>pi)
                for(int j=0;j<3;j++)
                    if(uv(j)(a)<0){
                        elem(i)(j)+=(a+1)*n;}
        RANGE<VECTOR<T,2> > bb2=RANGE<VECTOR<T,2> >::Bounding_Box(pts.Subset(elem(i)));}

    ARRAY<int> map(pts.m);
    map.Subset(elem.Flattened()).Fill(1);
    int k=0;
    for(int i=0;i<map.m;i++){
        pts(k)=pts(i);
        map(i)=map(i)?k++:-1;}
    pts.Resize(k);
    LOG::printf("range %P\n",RANGE<VECTOR<T,2> >::Bounding_Box(pts));
    elem.Flattened()=map.Subset(elem.Flattened());
    RANGE<VECTOR<T,2> > texture_range=RANGE<VECTOR<T,2> >::Bounding_Box(pts);
    for(int i=0;i<pts.m;i++) pts(i)=(pts(i)-texture_range.min_corner)/texture_range.Edge_Lengths();
    FILE_UTILITIES::Write_To_File(STREAM_TYPE(1.0f),"adaptive_torus_float.uv",pts,0,elem);

    IV2 size(500,500);
    GRID<VECTOR<T,2> > grid(size,RANGE<VECTOR<T,2> >::Unit_Box());
    ARRAY<TV,IV2> image(size);
    for(NODE_ITERATOR<VECTOR<T,2> > it(grid);it.Valid();it.Next()){
        VECTOR<T,2> X=it.Location();
        VECTOR<T,2> Y=X*texture_range.Edge_Lengths()+texture_range.min_corner+2*pi;
        VECTOR<T,2> Z=Y/(2*pi)*VECTOR<T,2>(8,4);
        int parity=IV2(Z).Sum()%2;
        TV color(parity,0,1-parity);
        image(it.index)=color;}

    PNG_FILE<T>::Write("adaptive_torus_float_texture.png",image);
    EPS_FILE<T> eps("test.eps",RANGE<VECTOR<T,2> >::Unit_Box()*VECTOR<T,2>(size));
    eps.Draw_Object(RANGE<VECTOR<T,2> >::Unit_Box());
    for(int i=0;i<elem.m;i++)
        eps.Draw_Object(VECTOR<VECTOR<T,2>,3>(pts.Subset(elem(i))));

    return 0;
}
