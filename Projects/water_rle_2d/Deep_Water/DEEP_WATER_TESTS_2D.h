#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEEP_WATER_TESTS_2D
//#####################################################################
#ifndef __DEEP_WATER_TESTS_2D__
#define __DEEP_WATER_TESTS_2D__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
namespace PhysBAM{

template<class T_input>
class DEEP_WATER_TESTS_2D:public WATER_STANDARD_TESTS_2D<RLE_GRID_2D<T_input> >
{
    typedef T_input T;
    typedef RLE_GRID_2D<T> T_GRID;
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::HORIZONTAL_GRID T_HORIZONTAL_GRID;typedef typename GRID_ARRAYS_POLICY<T_HORIZONTAL_GRID>::ARRAYS_SCALAR T_ARRAYS_HORIZONTAL_T;
    typedef typename T_GRID::FACE_Y_ITERATOR FACE_Y_ITERATOR;typedef typename T_HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;
    typedef typename TV_HORIZONTAL::template REBIND<int>::TYPE TV_HORIZONTAL_INT;
public:
    typedef WATER_STANDARD_TESTS_2D<T_GRID> BASE;
    using BASE::grid;using BASE::test_number;using BASE::rigid_body_collection;
    
    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example;
    bool use_deep_water;
    T deep_water_level;
    int sphere;

    DEEP_WATER_TESTS_2D(SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example_input,FLUIDS_PARAMETERS<T_GRID>& fluids_parameters_input,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const int test_number_input,const int resolution)
        :WATER_STANDARD_TESTS_2D<T_GRID>(example_input,fluids_parameters_input,rigid_body_collection_input,Shallow_Test_Number(test_number_input),resolution),
        example(example_input),use_deep_water(true),deep_water_level(0),sphere(0)
    {
        test_number=test_number_input;
        LOG::cout<<"Running deep water test number "<<test_number<<" at resolution "<<resolution<<std::endl;

        if(test_number==2){}
        else if(test_number==7){
            example.first_frame=0;example.last_frame=100;
            grid.Initialize(50*resolution+1,10*resolution+1,0,5,0,1);}
        else if(test_number==8){
            example.first_frame=0;example.last_frame=500;
            grid.Initialize(50*resolution+1,10*resolution+1,0,5,0,1);}
        else{
            LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}

        //use_deep_water=false;

        example.output_directory=STRING_UTILITIES::string_sprintf("Deep_Water/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));
    }

    static int Shallow_Test_Number(const int test_number)
    {return test_number>6?1:test_number;}

    void Initialize_Advection(const bool always_use_objects=false) PHYSBAM_OVERRIDE
    {bool use_objects=always_use_objects || test_number==7;
    BASE::Initialize_Advection(use_objects);}

    T Initial_Phi(const TV& X) const
    {if(test_number==7) return X.y-(T).512;
    else if(test_number==8) return X.y-(T).512;
    return BASE::Initial_Phi(X);}

    T Initial_Phi_Object(const TV& X) const
    {if(test_number==7) return rigid_body_collection.Rigid_Body(sphere).Implicit_Geometry_Extended_Value(X);
    return BASE::Initial_Phi_Object(X);}

    void Initialize_Bodies() PHYSBAM_OVERRIDE
    {if(test_number==7){
        sphere=rigid_body_collection.Add_Rigid_Body(stream_type,example.data_directory+"/Rigid_Bodies_2D/circle",(T).3,true,true,false);
        rigid_body_collection.rigid_body_particle.X(sphere)=VECTOR<T,2>((T)2,(T).5);
        rigid_body_collection.Rigid_Body(sphere).is_kinematic=true;}}

    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
    {if(test_number==7 && id==sphere){
        TV start(2,.5),V(1,0);
        T stop_time=.5;// 2.5;
        INTERPOLATION_CURVE<T,TV> motion_curve;
        motion_curve.Add_Control_Point(0,start);
        motion_curve.Add_Control_Point(stop_time,start+stop_time*V);
        frame.t=motion_curve.Value(time);}}

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
    {return false;}

    void Limit_Dt(T& dt,const T time)
    {for(int r=1;r<=rigid_body_collection.rigid_body_particle.array_collection->Size();r++){
        TV velocity=rigid_body_collection.rigid_body_particle.twist(r).linear;
        T rigid_dt_denominator=Dot_Product(abs(velocity),grid.One_Over_DX());
        dt=min(dt,Robust_Inverse(rigid_dt_denominator));}
    BASE::Limit_Dt(dt,time);}

    void Get_Surface_Pressure(const T dt,const T time)
    {const T_GRID& grid=*example.fluids_parameters.grid;
    const ARRAY<T>& phi=example.particle_levelset.phi;
    ARRAY<T>& u_interface=example.incompressible.projection.laplace.u_interface;
    if(test_number==8){
        T strength=2,size=1,center=2+time;
        T grid_length=grid.horizontal_grid.Domain().Edge_Lengths().x;
        ARRAY<T,VECTOR<int,1> > surface_pressure(grid.horizontal_grid.Get_MAC_Grid().Domain_Indices(grid.number_of_ghost_cells));
        for(HORIZONTAL_CELL_ITERATOR column(grid.horizontal_grid,1);column.Valid();column.Next()){
            T x=column.Location().x;
            T offset=fmod(x-center,grid_length);
            if(offset<0) offset+=grid_length;
            if(offset>grid_length/2) offset-=grid_length;
            surface_pressure(column.Cell_Index())=strength*sqrt(max((T)0,1-sqr(offset/size)));}
        surface_pressure*=dt;
        u_interface.Fill(0);
        for(typename T_GRID::FACE_X_ITERATOR face(grid,0);face;face++)if(LEVELSET_UTILITIES<T>::Interface(phi(face.cell1.Cell()),phi(face.cell2.Cell()))){
            T theta=LEVELSET_UTILITIES<T>::Theta(phi(face.cell1.Cell()),phi(face.cell2.Cell()));
            u_interface(face.Face())=(1-theta)*surface_pressure(face.cell1.i)+theta*surface_pressure(face.cell2.i);}
        for(FACE_Y_ITERATOR face(grid,0);face;face++)
            u_interface(face.Face())=surface_pressure(face.i());}}

//#####################################################################
};
}
#endif
#endif
