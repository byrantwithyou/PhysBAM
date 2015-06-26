//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Images/EPS_FILE.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_MULTIPLE_UNIFORM.h>

using namespace PhysBAM;

template<class TV>
void Draw_Bubble(PARSE_ARGS& parse_args)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    STREAM_TYPE stream_type((T()));
    int frame=1;
    T line_width=.01;
    bool depressurize=false;
    std::string sim_dir,base_filename;
    TV_INT size(500,500);
    INTERVAL<T> pressure_interval=INTERVAL<T>::Empty_Box();
    TV normalize_point;
    T rho=(T)2;
    T dx=(T)2/(T)64;
    bool normalize_pressure=false,force_dims=true;

    parse_args.Add("-frame",&frame,"frame","Frame to draw");
    parse_args.Add("-size",&size,"size","Image size");
    parse_args.Add("-p_min",&pressure_interval.min_corner,"value","Pressure minimum for image");
    parse_args.Add("-p_max",&pressure_interval.max_corner,"value","Pressure maximum for image");
    parse_args.Add("-norm_point",&normalize_point,&normalize_pressure,"location","Point to center pressure");
    parse_args.Add("-depressurize",&depressurize,"Subtract out hydrostatic pressure");
    parse_args.Add("-rho",&rho,"rho","Density");
    parse_args.Add("-line_width",&line_width,"line_width","Line Width");
    parse_args.Add("-dx",&dx,"dx","Spatial Resolution");
    parse_args.Add("-force_dims",&force_dims,"Force y dimension so image is proportional to grid");
    parse_args.Extra(&sim_dir,"sim dir","simulation directory");
    parse_args.Extra(&base_filename,"filename","Base filename for output images");
    parse_args.Parse();

    GRID<TV> grid;
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/common/grid",sim_dir.c_str()),grid);
    std::cout<<grid.domain<<size<<std::endl;
    if(force_dims) size.y=(int)(grid.domain.Edge_Lengths().y/grid.domain.Edge_Lengths().x*size.x);

    ARRAY<T,TV_INT> pressure;
    FILE_UTILITIES::Read_From_File(stream_type,LOG::sprintf("%s/%d/pressure",sim_dir.c_str(),frame),pressure);

    if(depressurize)
    {
        for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next())
        {
            pressure(it.index)+=rho*dx*it.index.y*9.8;
        }
    }


    T max_phi=grid.domain.Edge_Lengths().Magnitude();
    ARRAY<T,TV_INT> best_phi(grid.Domain_Indices(),true,max_phi);
    ARRAY<int,TV_INT> best_color(grid.Domain_Indices(),true,-1);
    ARRAY<LEVELSET<TV>*> levelsets;
    ARRAY<ARRAY<T,TV_INT>*> phis;

    {
        EPS_FILE<T> eps_writer(base_filename+".eps",RANGE<TV>(TV(),TV(size)));
        eps_writer.Use_Fixed_Bounding_Box(grid.domain);
        eps_writer.cur_format.line_width=line_width;

        for(int i=0;;i++){
            ARRAY<T,TV_INT>* phi=new ARRAY<T,TV_INT>;
            LEVELSET<TV>* ls=new LEVELSET<TV>(grid,*phi);
            try{FILE_UTILITIES::Read_From_File<T>(LOG::sprintf("%s/%d/levelset_%d.gz",sim_dir.c_str(),frame,i),*ls);}
            catch(...){delete ls;break;}
            phis.Append(phi);
            levelsets.Append(ls);

            for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next())
                if(ls->phi(it.index)<best_phi(it.index)){
                    best_phi(it.index)=ls->phi(it.index);
                    best_color(it.index)=i;}

            SEGMENTED_CURVE_2D<T>& sc=*SEGMENTED_CURVE_2D<T>::Create();
            MARCHING_CUBES<TV>::Create_Surface(sc,levelsets(i)->grid,levelsets(i)->phi);

            for(int t=0;t<sc.mesh.elements.m;t++)
                eps_writer.Draw_Object(sc.particles.X(sc.mesh.elements(t)(0)),sc.particles.X(sc.mesh.elements(t)(1)));}
    }

    CUBIC_MN_INTERPOLATION_UNIFORM<TV,T> interp;
    ARRAY<ARRAY<T,TV_INT> > color_pressure(levelsets.m);
    for(int i=0;i<levelsets.m;i++){
        color_pressure(i)=pressure;
        EXTRAPOLATION_HIGHER_ORDER<TV,T> eho(grid,*levelsets(i),20,3,4);
        eho.periodic=true;
        eho.Extrapolate_Cell([&](const TV_INT& index){return best_color(index)==i && best_phi(index)<grid.dX.Max()*(T).01;},color_pressure(i));}

    INTERPOLATED_COLOR_MAP<T> cm;
    if(pressure_interval.Empty()) pressure_interval=INTERVAL<T>(pressure.Min(),pressure.Max());
    cm.mn=pressure_interval.min_corner;
    cm.mx=pressure_interval.max_corner;
    cm.colors.Add_Control_Point(cm.mn,VECTOR<T,3>(0,0,0));
    cm.colors.Add_Control_Point(cm.mx,VECTOR<T,3>(1,1,1));

    INTERVAL<T> p_range=INTERVAL<T>::Empty_Box();
    ARRAY<VECTOR<T,3>,TV_INT> pressure_image(size);
    GRID<TV> image_grid(size,grid.domain,true);
    T p_shift=0;
    if(normalize_pressure){
        T best=max_phi;
        int best_index=-1;
        for(int i=0;i<levelsets.m;i++){
            T p=levelsets(i)->Extended_Phi(normalize_point);
            if(p<best){best=p;best_index=i;}}
        T p=interp.Periodic(grid,color_pressure(best_index),normalize_point);
        p_shift=-p;}

    for(CELL_ITERATOR<TV> it(image_grid);it.Valid();it.Next()){
        T best=max_phi;
        int best_index=-1;
        TV X=it.Location();
        for(int i=0;i<levelsets.m;i++){
            T p=levelsets(i)->Extended_Phi(X);
            if(p<best){best=p;best_index=i;}}
        T p=interp.Periodic(grid,color_pressure(best_index),X)+p_shift;
        if(best<=0) p_range.Enlarge_To_Include_Point(p);else p=cm.mx;
        pressure_image(it.index)=cm(p);}
    PNG_FILE<T>::Write(base_filename+".png",pressure_image);
    LOG::cout<<"pressure range: "<<p_range<<std::endl;

    int ret=system(LOG::sprintf("convert %s.png %s.eps -composite %s-full.png",base_filename.c_str(),base_filename.c_str(),base_filename.c_str()).c_str());
    PHYSBAM_ASSERT(!ret);
    phis.Delete_Pointers_And_Clean_Memory();
    levelsets.Delete_Pointers_And_Clean_Memory();
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


