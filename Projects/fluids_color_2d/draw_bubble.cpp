//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// TBD
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Tools/Images/PNG_FILE.h>
#include <PhysBAM_Geometry/Images/EPS_FILE.h>


using namespace PhysBAM;

template<class TV>
void Draw_Bubble(PARSE_ARGS& parse_args)
{
    
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    STREAM_TYPE stream_type((T()));
    GRID<TV> grid;
    ARRAY<T,VECTOR<int,2> > phi;
    int frame=1;
    int expand_factor=1;
    int number_of_colors=2;
    std::string sim_dir,filename,filename_image,filename_eps_image;

    parse_args.Add("-frame",&frame,"frame","Frame to draw");
    parse_args.Add("-colors",&frame,"colors","Number of colors in sim");
    parse_args.Extra_Optional(&sim_dir,"sim dir","simulation directory");
    parse_args.Parse();

    filename=STRING_UTILITIES::string_sprintf("%s/%d/levelset_0.gz",sim_dir.c_str(),frame);
    filename_image=STRING_UTILITIES::string_sprintf("russ_image",sim_dir.c_str(),frame);
    filename_eps_image=STRING_UTILITIES::string_sprintf("russ.eps",sim_dir.c_str(),frame);
    

    // eps_writer.Draw_Object(TV(1,2),TV(40,30),TV(22,5));

    ARRAY<LEVELSET<TV>*> levelsets(number_of_colors);
    ARRAY<T,TV_INT> pressure;
    ARRAY<VECTOR<T,3>, VECTOR<int,2> > pressure_image;
    
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/pressure",sim_dir.c_str(),frame),pressure);

    EPS_FILE<T> eps_writer(filename_eps_image,RANGE<TV>(pressure.domain-pressure.domain.min_corner));
    std::cout<<eps_writer.output_box<<" "<<eps_writer.bounding_box<<std::endl;

//         std::cout<<pressure.Min()<<" "<<pressure.Max()<<std::endl;
         pressure-=pressure.Min();
         //       std::cout<<pressure.Min()<<" "<<pressure.Max()<<std::endl;
            pressure/=pressure.Max();
            //std::cout<<pressure.Min()<<" "<<pressure.Max()<<std::endl;
            pressure_image.Resize(pressure.domain*expand_factor);
            pressure_image.Fill(VECTOR<T,3>(1,0,0));
            //std::cout<<pressure_image.domain<<std::endl;
             pressure_image*=pressure;
             // for(int i=-5;i<37;i++)for(int j=-5;j<69;j++)std::cout<<i<<" "<<j<<" "<<pressure_image(i,j)<<std::endl;

            PNG_FILE<T>::Write(filename_image,pressure_image);    
            std::cout<<RANGE<TV>(pressure.domain)<<std::endl;
    std::cout<<eps_writer.output_box<<" "<<eps_writer.bounding_box<<std::endl;

    for(int i=0;i<number_of_colors;i++)
    {
        levelsets(i)=new LEVELSET<TV>(grid,phi);
        FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/%d/levelset_%d.gz",sim_dir.c_str(),frame,i),*(levelsets(i)));
        SEGMENTED_CURVE_2D<T>& sc=*SEGMENTED_CURVE_2D<T>::Create();
        MARCHING_CUBES<TV>::Create_Surface(sc,levelsets(i)->grid,levelsets(i)->phi);
        
        std::cout<<"This levelset has "<<sc.mesh.elements.m<<" elements"<<std::endl;
        for (int t=0;t<sc.mesh.elements.m;t++){
        int node1=sc.mesh.elements(t)(0),node2=sc.mesh.elements(t)(1);
//        std::cout<<sc.particles.X(node1)<<" "<<sc.particles.X(node2)<<std::endl;
        eps_writer.Draw_Object(sc.particles.X(node1),sc.particles.X(node2));
        }
            

        

        eps_writer.bounding_box=RANGE<TV>(TV(-1,0),TV(1,6));


    }


}

int main(int argc,char *argv[])
{
    bool use_double=false;
    PARSE_ARGS parse_args(argc,argv);
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;
    parse_args.Add("-double",&use_double,"use doubles");
    parse_args.Parse(true);
    
    if(use_double) Draw_Bubble<VECTOR<double,2> >(parse_args);
    else  Draw_Bubble<VECTOR<float,2> >(parse_args);
}

