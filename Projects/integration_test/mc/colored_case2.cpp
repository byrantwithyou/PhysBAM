//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_EDGE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
using namespace PhysBAM;

typedef double T;
typedef float RW;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,2> TV_INT;
typedef VECTOR<T,3> TV3;

std::string output_directory;

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
    srand(time(0));
    Get_Debug_Particles<TV>();
    PARSE_ARGS parse_args;
    parse_args.Add_String_Argument("-o","output","output directory");
    parse_args.Add_String_Argument("-c","aaaa","color case");
    parse_args.Add_Integer_Argument("-n",5,"extra random tests");
    parse_args.Parse(argc,argv);
    output_directory=parse_args.Get_String_Value("-o");
    std::string colors=parse_args.Get_String_Value("-c");
//    int n=parse_args.Get_Integer_Value("-n");

    ARRAY<TV3> color_map(128);
    color_map(0)=TV3(1,0,0);
    color_map(1)=TV3(0,1,0);
    color_map(2)=TV3(0,0,1);
    color_map(3)=TV3(1,1,0);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");

    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();

    RANDOM_NUMBERS<T> random;

    GRID<TV> grid(TV_INT()+40,RANGE<TV>::Centered_Box(),true);
    ARRAY<T,TV_INT> phi(grid.Node_Indices());
    ARRAY<int,TV_INT> color(grid.Node_Indices());

    struct ANALYTIC_TEST
    {
        T r;
        TV n; // rotate angle with respect to e_y
        VECTOR<TV,3> centers;
        VECTOR<TV,3> normals;
        VECTOR<VECTOR<int,3>,3> sectors;
        ARRAY<T> mu;
        ANALYTIC_TEST()
        {
            mu.Append(1);mu.Append(2);mu.Append(3);mu.Append(4);
            r=(T).4;
            n=TV::Axis_Vector(1);
            centers(0)=TV::Axis_Vector(1);
            centers(1).x=(T)sqrt(3)/2;centers(1).y=-(T)1/2;
            centers(2).x=-(T)sqrt(3)/2;centers(2).y=-(T)1/2;
            for(int i=0;i<3;i++){
                normals(i).x=-centers(i).y;
                normals(i).y=centers(i).x;
                centers(i)*=r;
                for(int j=0;j<3;j++) sectors(i)(j)=(i+j)%3;}
        }
        T Phi(const TV& X) const
        {
//            return abs(X.Max_Abs()-.5111+X.Magnitude()/10);
            



            TV x=X-.02;
            int i;
            for(i=0;i<2;i++) if(x.Dot(normals(sectors(i)(1)))>=0 && x.Dot(normals(sectors(i)(2)))<0) break;
            T d=(x-centers(i)).Magnitude();
            if(d>r && x.Magnitude()>r/100) return d-r;
            else return min(abs(d-r),abs(x.Dot(normals(sectors(i)(1)))),abs(x.Dot(normals(sectors(i)(2)))));
        }
        int Color(const TV& X) const
        {
//            return (X.Max_Abs()-.5111+X.Magnitude()/10)<0;
            TV x=X-.02;
            int i;
            for(i=0;i<2;i++) if(x.Dot(normals(sectors(i)(1)))>=0 && x.Dot(normals(sectors(i)(2)))<0) break;
            T d=(x-centers(i)).Magnitude();
            if(d>r && x.Magnitude()>r/100) return 0;
            else return i+1;
        }
    } test;

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        phi(it.index)=test.Phi(it.Location());
        color(it.index)=test.Color(it.Location());}

    HASHTABLE<VECTOR<int,2>,SEGMENTED_CURVE_2D<T>*> surface;
    HASHTABLE<int,SEGMENTED_CURVE_2D<T>*> boundary;

    MARCHING_CUBES_COLOR<TV>::Get_Elements(grid,surface,boundary,color,phi);

    for(HASHTABLE<VECTOR<int,2>,SEGMENTED_CURVE_2D<T>*>::ITERATOR it(surface);it.Valid();it.Next()){
        const SEGMENTED_CURVE_2D<T>& curve=*it.Data();
        VECTOR<int,2> color=it.Key();
        for(int e=0;e<curve.mesh.elements.m;e++){
            SEGMENT_2D<T> seg=curve.Get_Element(e);
            Add_Debug_Object(seg.X-seg.Normal()*(T).004,color_map(color.y));
            Add_Debug_Object(seg.X+seg.Normal()*(T).004,color_map(color.x));}}

    for(UNIFORM_GRID_ITERATOR_EDGE<TV> it(grid);it.Valid();it.Next()){
        TV A=grid.Node(it.Full_Index().First_Node_Index());
        TV B=grid.Node(it.Full_Index().Second_Node_Index());
        VECTOR<TV,2> C(A,B);
        Add_Debug_Object(C,TV3()+(T).5);
    }

    Dump_Frame("Surface");

    return 0;
}
