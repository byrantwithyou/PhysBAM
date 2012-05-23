//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
using namespace PhysBAM;

typedef double T;
typedef float RW;
typedef VECTOR<T,2> TV;
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
    int n=parse_args.Get_Integer_Value("-n");

    ARRAY<TV3> color_map(128);
    color_map('a')=TV3(1,0,0);
    color_map('b')=TV3(0,1,0);
    color_map('c')=TV3(0,0,1);
    color_map('d')=TV3(1,1,0);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");

    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();

    VECTOR<T,4> phi;
    RANDOM_NUMBERS<T> random;

    for(int r=0;r<n;r++){
        ARRAY<TRIPLE<SEGMENT_2D<T>,int,int> > surface;
        ARRAY<PAIR<SEGMENT_2D<T>,int> > boundary;
        VECTOR<int,4> color_vector;
        for(int i=0;i<4;i++) color_vector(i)=colors[i];
        random.Fill_Uniform(phi,0.01,1);
        MARCHING_CUBES_COLOR<TV>::Get_Elements_For_Cell(surface,boundary,color_vector,phi);

        for(int i=0;i<surface.m;i++){
            Add_Debug_Object(surface(i).x.X-surface(i).x.Normal()*(T).01,color_map(surface(i).z));
            Add_Debug_Object(surface(i).x.X+surface(i).x.Normal()*(T).01,color_map(surface(i).y));}
        for(int i=0;i<boundary.m;i++){
            Add_Debug_Object(boundary(i).x.X-boundary(i).x.Normal()*(T).01,color_map(boundary(i).y));
            Add_Debug_Object(boundary(i).x.X+boundary(i).x.Normal()*(T).01,TV3((T).5,(T).5,(T).5));}
        for(int i=0;i<4;i++){
            TV X(i&1,i/2&1);
            Add_Debug_Particle(X,color_map(colors[i]));}
        if(0)
        for(int a=0,k=0;a<2;a++){
            int mask=1<<a;
            for(int v=0;v<4;v++)
                if(!(v&mask)){
                    int w=v|mask;
                    Add_Debug_Object(VECTOR<TV,2>(TV(v&1,v/2&1),TV(w&1,w/2&1)),(T).4*TV3(1,1,1));
                    k++;}}
        Dump_Frame(colors.c_str());
        for(int i=0;i<4;i++)
            colors[i]=rand()%4+'a';}
    
    Dump_Frame("Surface");

    return 0;
}
