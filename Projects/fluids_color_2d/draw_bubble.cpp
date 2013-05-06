//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Images/PNG_FILE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Images/EPS_FILE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_MULTIPLE_UNIFORM.h>

using namespace PhysBAM;

template<class TV>
void Draw_Bubble(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    STREAM_TYPE stream_type((T()));
    ARRAY<T,TV_INT> phi;
    int frame=1;
    std::string sim_dir,base_filename;
    TV_INT size(500,500);

    parse_args.Add("-frame",&frame,"frame","Frame to draw");
    parse_args.Add("-size_x",&size.x,"size","Image size");
    parse_args.Add("-size_y",&size.y,"size","Image size");
    parse_args.Extra(&sim_dir,"sim dir","simulation directory");
    parse_args.Extra(&base_filename,"filename","Base filename for output images");
    parse_args.Parse();

    GRID<TV> grid;
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/common/grid",sim_dir.c_str()),grid);

    ARRAY<T,TV_INT> pressure;
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/pressure",sim_dir.c_str(),frame),pressure);

    INTERPOLATED_COLOR_MAP<T> cm;
    cm.mn=pressure.Min();
    cm.mx=pressure.Max();
    cm.colors.Add_Control_Point(cm.mn,VECTOR<T,3>(0,0,0));
    cm.colors.Add_Control_Point(cm.mx,VECTOR<T,3>(1,1,1));

    ARRAY<VECTOR<T,3>,TV_INT> pressure_image(size);
    GRID<TV> image_grid(size,grid.domain,true);
    CUBIC_MN_INTERPOLATION_UNIFORM<GRID<TV>,T> interp;
    for(CELL_ITERATOR<TV> it(image_grid);it.Valid();it.Next())
        pressure_image(it.index)=cm(interp.Periodic(grid,pressure,it.Location()));
    PNG_FILE<T>::Write(base_filename+".png",pressure_image);

    EPS_FILE<T> eps_writer(base_filename+".eps",grid.domain);
    ARRAY<LEVELSET<TV>*> levelsets;
    for(int i=0;;i++){
        LEVELSET<TV>* ls=new LEVELSET<TV>(grid,phi);
        try{FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/%d/levelset_%d.gz",sim_dir.c_str(),frame,i),*ls);}
        catch(...){delete ls;break;}
        levelsets.Append(ls);

        SEGMENTED_CURVE_2D<T>& sc=*SEGMENTED_CURVE_2D<T>::Create();
        MARCHING_CUBES<TV>::Create_Surface(sc,levelsets(i)->grid,levelsets(i)->phi);

        for(int t=0;t<sc.mesh.elements.m;t++)
            eps_writer.Draw_Object(VECTOR<TV,2>(sc.particles.X.Subset(sc.mesh.elements(t))));}
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
    else Draw_Bubble<VECTOR<float,2> >(parse_args);
}

