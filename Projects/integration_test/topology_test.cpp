//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
using namespace PhysBAM;

typedef double T;
typedef float RW;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,3> TV_INT;

#define rdtscll(val) do {  \
        unsigned int __a,__d;                        \
        asm volatile("rdtsc" : "=a" (__a), "=d" (__d));                 \
        (val) = ((unsigned long long)__a) | (((unsigned long long)__d)<<32); \
    } while(0)
inline unsigned long long rdtsc(){unsigned long long x;rdtscll(x);return x;}

std::string output_directory="output";

template<class TV> DEBUG_PARTICLES<TV>& Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}

void Dump_Frame(const char* title)
{
    static int frame=0;
    char buff[100];
    sprintf(buff, "%s/%i", output_directory.c_str(), frame);
    FILE_UTILITIES::Create_Directory(buff);
//    FILE_UTILITIES::Write_To_File<RW>((std::string)buff+"/mac_velocities.gz",u);
    if(title) FILE_UTILITIES::Write_To_Text_File((std::string)buff+"/frame_title",title);
    Get_Debug_Particles<TV>().Write_Debug_Particles(STREAM_TYPE((RW())),output_directory,frame);
    frame++;
}

int main(int argc, char* argv[])
{
    typedef typename MARCHING_CUBES_COLOR<TV>::INTERFACE_ELEMENT INTERFACE_ELEMENT;
    typedef typename MARCHING_CUBES_COLOR<TV>::BOUNDARY_ELEMENT BOUNDARY_ELEMENT;
    srand(time(0));
    Get_Debug_Particles<TV>();
    std::string colors="abbbbbbb";
    int n=5,s=5;
    bool use_seed=false;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-o",&output_directory,"output","output directory");
    parse_args.Add("-c",&colors,"color","color case");
    parse_args.Add("-n",&n,"num","extra random tests");
    parse_args.Add("-s",&s,&use_seed,"seed","random seed");
    parse_args.Parse();

    ARRAY<TV> color_map;
    color_map.Append(TV(1,0,0));
    color_map.Append(TV(1,(T).5,0));
    color_map.Append(TV(1,1,0));
    color_map.Append(TV(0,1,0));
    color_map.Append(TV((T).5,1,(T).5));
    color_map.Append(TV(0,1,1));
    color_map.Append(TV(0,0,1));
    color_map.Append(TV(1,0,1));
    color_map.Append(TV(1,1,1));
    color_map.Append(TV((T).5,(T).5,(T).5));

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");

    RANDOM_NUMBERS<T> random;
    if(use_seed) random.Set_Seed(s);
    GRID<TV> grid(TV_INT()+n,RANGE<TV>::Unit_Box());
    ARRAY<T,TV_INT> phi(grid.Node_Indices());
    phi.Fill(1);
    ARRAY<int,TV_INT> color(grid.Node_Indices());
    color.Fill(0);
    for(NODE_ITERATOR<TV> it(grid,-2);it.Valid();it.Next()){
        phi(it.index)=random.Get_Uniform_Number(0.1,1);
        color(it.index)=random.Get_Uniform_Integer(0,0);}

    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();

    ARRAY<INTERFACE_ELEMENT> surface;
    ARRAY<BOUNDARY_ELEMENT> boundary;

    HASHTABLE<VECTOR<TV_INT,2>,ARRAY<VECTOR<int,2> > > hash;

    unsigned long long t0=rdtsc();
    const VECTOR<TV_INT,8>& bits=GRID<TV>::Binary_Counts(TV_INT());
    for(NODE_ITERATOR<TV> it(grid,-1);it.Valid();it.Next()){
        VECTOR<T,8> p;
        VECTOR<int,8> c;
        for(int i=0;i<8;i++){
            p(i)=phi(it.index+bits(i));
            c(i)=color(it.index+bits(i));}

        surface.Remove_All();
        boundary.Remove_All();
        MARCHING_CUBES_COLOR<TV>::Get_Elements_For_Cell(surface,boundary,c,p);

        for(int i=0;i<surface.m;i++){
            surface(i).face.X*=grid.dX;
            surface(i).face.X+=it.Location();
//            Add_Debug_Object(surface(i).face.X,color_map(surface(i).y),color_map(surface(i).z));

            VECTOR<int,2> col(surface(i).color_pair.x,surface(i).color_pair.y);
            TV_INT a(surface(i).face.X(0)*100000000);
            TV_INT b(surface(i).face.X(1)*100000000);
            TV_INT c(surface(i).face.X(2)*100000000);
            hash.Get_Or_Insert(VECTOR<TV_INT,2>(a,b)).Append(col);
            hash.Get_Or_Insert(VECTOR<TV_INT,2>(b,c)).Append(col);
            hash.Get_Or_Insert(VECTOR<TV_INT,2>(c,a)).Append(col);}}
    unsigned long long t1=rdtsc();
    printf("setup: %.2f\n", (t1-t0)/3059.107e6);

    for(HASHTABLE<VECTOR<TV_INT,2>,ARRAY<VECTOR<int,2> > >::ITERATOR it(hash);it.Valid();it.Next()){
        ARRAY<int> a,b;
        ARRAY<VECTOR<int,2> >& c=it.Data();
        for(int i=0;i<c.m;i++){
            a.Append(c(i).x);
            b.Append(c(i).y);}
        ARRAY<VECTOR<int,2> >* d=hash.Get_Pointer(it.Key().Reversed());
        if(d){
            for(int i=0;i<d->m;i++){
                b.Append((*d)(i).x);
                a.Append((*d)(i).y);}}
        a.Sort();
        b.Sort();
        if(a!=b) LOG::cout<<it.Key()<<"   "<<a<<"    "<<b<<std::endl;}

    Dump_Frame("A");


    return 0;
}
