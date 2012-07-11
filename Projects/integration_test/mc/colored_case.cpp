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
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_COLOR.h>
using namespace PhysBAM;

typedef double T;
typedef float RW;
typedef VECTOR<T,3> TV;

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
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_String_Argument("-o","output","output directory");
    parse_args.Add_String_Argument("-c","abbbbbbb","color case");
    parse_args.Add_Integer_Argument("-n",5,"extra random tests");
    parse_args.Parse();
    output_directory=parse_args.Get_String_Value("-o");
    std::string colors=parse_args.Get_String_Value("-c");
    int n=parse_args.Get_Integer_Value("-n");

    ARRAY<TV> color_map(128);
    color_map('a')=TV(1,0,0);
    color_map('b')=TV(1,1,0);
    color_map('c')=TV(0,1,0);
    color_map('d')=TV(0,1,1);
    color_map('e')=TV(0,0,1);
    color_map('f')=TV(1,0,1);
    color_map('g')=TV(1,1,1);
    color_map('h')=TV((T).5,(T).5,(T).5);

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");

    MARCHING_CUBES_COLOR<TV>::Initialize_Case_Table();

    VECTOR<T,8> phi;
    RANDOM_NUMBERS<T> random;

    for(int r=0,rr=0;r<n;r++,rr++){
        if(rr>0 && rr%1000==0) printf("====================================================== %i ======================================================\n",rr);
        ARRAY<TRIPLE<TRIANGLE_3D<T>,int,int> > surface;
        ARRAY<PAIR<TRIANGLE_3D<T>,int> > boundary;
        VECTOR<int,8> color_vector;
        for(int i=0;i<8;i++) color_vector(i)=colors[i];
        random.Fill_Uniform(phi,0.01,1);
        MARCHING_CUBES_COLOR<TV>::Get_Elements_For_Cell(surface,boundary,color_vector,phi);

        bool show=false;
        if(0)
        for(int i=0;i<surface.m;i++)
            for(int j=0;j<surface.m;j++)
                for(int k=0;k<3;k++){
                    SEGMENT_3D<T> seg(surface(j).x.X(k),surface(j).x.X((k+1)%3));
                    bool touch=false;
                    for(int r=0;r<3;r++)
                        for(int s=0;s<2;s++)
                            if((surface(i).x.X(r)-seg.X(s)).Magnitude()<1e-12)
                                touch=true;
                    if(touch) continue;
                    T a;
                    TV weights;
                    if(INTERSECTION::Intersects(seg,surface(i).x,a,weights,(T)0)){
//                        printf("weights: %g %g %g   %g\n", weights.x, weights.y, weights.z, a);
                        if(weights.Min()>1e-12 && a>1e-12 && 1-a>1e-12){
                            printf("weights: %g %g %g   %g\n", weights.x, weights.y, weights.z, a);
                            show=true;}}
                }
        if(0){
            if(!show){
                r--;
                for(int i=0;i<8;i++)
                    colors[i]=rand()%8+'a';
                continue;}
            else printf("FOUND: %i (%i)\n", r, rr);}

        for(int i=0;i<surface.m;i++)
            Add_Debug_Object(surface(i).x.X,color_map(surface(i).y),color_map(surface(i).z));
        for(int i=0;i<boundary.m;i++)
            Add_Debug_Object(boundary(i).x.X,color_map(boundary(i).y)*(T).2+(T).8,color_map(boundary(i).y));
        for(int i=0;i<8;i++){
            TV X(i&1,i/2&1,i/4&1);
            Add_Debug_Particle(X,color_map(colors[i]));}
        for(int a=0,k=0;a<3;a++){
            int mask=1<<a;
            for(int v=0;v<8;v++)
                if(!(v&mask)){
                    int w=v|mask;
                    Add_Debug_Object(VECTOR<TV,2>(TV(v&1,v/2&1,v/4&1),TV(w&1,w/2&1,w/4&1)),(T).4*TV(1,1,1));
                    k++;}}
        Dump_Frame(colors.c_str());
        for(int i=0;i<8;i++)
            colors[i]=rand()%8+'a';}
    
    Dump_Frame("Surface");

    return 0;
}
