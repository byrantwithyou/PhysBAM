//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VELOCITY_TESTS__
#define __VELOCITY_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>

namespace PhysBAM{

template<class TV>
class VELOCITY_TESTS:public INCOMPRESSIBLE_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef INCOMPRESSIBLE_EXAMPLE<TV> BASE;

    int test_number;
    bool use_conservative_advection;
    T kinetic_energy_gained,momentum_gained;
    STREAM_TYPE stream_type;

public:
    using BASE::mac_grid; using BASE::incompressible;using BASE::projection;using BASE::output_directory;using BASE::mpi_grid;using BASE::domain_boundary;
    using BASE::face_velocities;using BASE::last_frame;using BASE::write_substeps_level;using BASE::rigid_geometry_collection;using BASE::boundary_scalar;using BASE::density;
    using BASE::boundary;using BASE::restart;using BASE::temperature;

    VELOCITY_TESTS(const STREAM_TYPE stream_type_input,const PARSE_ARGS& parse_args)
        :INCOMPRESSIBLE_EXAMPLE<TV>(stream_type_input),use_conservative_advection(false),kinetic_energy_gained(0),momentum_gained(0),stream_type(stream_type_input)
    {
        last_frame=10;
        write_substeps_level=parse_args.Get_Integer_Value("-substep");;
        test_number=parse_args.Get_Integer_Value("-test_number");
        int scale=parse_args.Get_Integer_Value("-scale");
        restart=parse_args.Get_Integer_Value("-restart");
        use_conservative_advection=parse_args.Is_Value_Set("-conservative");
        output_directory=STRING_UTILITIES::string_sprintf("Velocity_Tests/Test_%d_%d_%d%s",test_number,scale,TV::dimension,use_conservative_advection?"_conservative":"");
        TV_INT counts=TV_INT::All_Ones_Vector()*scale;
        RANGE<TV> range=RANGE<TV>(TV(),TV::All_Ones_Vector());
        mac_grid.Initialize(counts,range,true);
        if(test_number==1) incompressible.Set_Viscosity(100);
        if(test_number==2) incompressible.Set_Viscosity(100);
        if(test_number==1) incompressible.Set_Gravity();
        if(use_conservative_advection){
            incompressible.Set_Custom_Advection(*new ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>());}
        boundary_scalar.Set_Fixed_Boundary(true,0);
        if(test_number==1){
            VECTOR<VECTOR<bool,2>,TV::dimension> extrapolation;
            for(int i=0;i<TV::dimension;i++){
                if(i==2){extrapolation(i)(1)=true;extrapolation(i)(2)=true;}
                else{extrapolation(i)(1)=false;extrapolation(i)(2)=false;}}
            boundary_scalar.Set_Constant_Extrapolation(extrapolation);}
    }

    ~VELOCITY_TESTS()
    {
        if(use_conservative_advection) delete incompressible.advection;
    }

    void Initialize_Bodies()
    {
    }
    
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
    {
    }

    void Write_Output_Files(const int frame)
    {
        std::cout<<std::endl<<"Energy Gained at f "<<frame<<" is "<<kinetic_energy_gained<<std::endl;
        // std::cout<<std::endl<<"Energy Gained from visc at f "<<frame<<" is "<<incompressible.kinetic_energy_gained<<std::endl;
        // kinetic_energy_gained+=incompressible.kinetic_energy_gained;incompressible.kinetic_energy_gained=0;
        if(frame==0) kinetic_energy_gained=0;
        T total_energy=0;
        for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            total_energy+=.5*face_velocities(iterator.Full_Index())*face_velocities(iterator.Full_Index())*mac_grid.min_dX;}
        std::cout<<std::endl<<"Total Energy is "<<total_energy<<std::endl;
        std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/energy_gained",kinetic_energy_gained);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/momentum_gained",momentum_gained);
        BASE::Write_Output_Files(frame);
    }

    void Set_Boundary_Conditions(const T time)
    {
        projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
        for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*axis+axis_side;
            if(domain_boundary(axis)(axis_side)){
                TV_INT interior_cell_offset=axis_side==0?TV_INT():-TV_INT::Axis_Vector(axis);    
                for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                    TV_INT boundary_face=axis_side==0?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                    if((axis!=2 && test_number==1) || ((axis==2 || (axis!=2 && iterator.Location()(2)<.9)) && test_number==2)){if(face_velocities.Component(axis).Valid_Index(boundary_face)) projection.elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(axis,boundary_face))=true;}
                    else {projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
                for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT boundary_face=axis_side==0?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                    if(((axis!=2 && test_number==1) || ((axis==2 || (axis!=2 && iterator.Location()(2)<.9)) && test_number==2)) && face_velocities.Component(axis).Valid_Index(boundary_face)) face_velocities(FACE_INDEX<TV::dimension>(axis,boundary_face))=0;}}}
        for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            if(test_number==2 && iterator.Location()(2)>.9){
                projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
                if(iterator.Axis()==1){
                    momentum_gained+=(1-face_velocities(iterator.Full_Index()))*mac_grid.dX(iterator.Axis());
                    kinetic_energy_gained+=.5*(1-face_velocities(iterator.Full_Index())*face_velocities(iterator.Full_Index()))*mac_grid.dX(iterator.Axis());
                    face_velocities(iterator.Full_Index())=1;}
                else{
                    momentum_gained-=face_velocities(iterator.Full_Index())*mac_grid.dX(iterator.Axis());
                    kinetic_energy_gained+=.5*(0-face_velocities(iterator.Full_Index())*face_velocities(iterator.Full_Index()))*mac_grid.dX(iterator.Axis());
                    face_velocities(iterator.Full_Index())=0;}}}
    }

    void Initialize_Fields()
    {
        for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=temperature(iterator.Cell_Index())=0;
    }

    void Get_Scalar_Field_Sources(const T time)
    {
    }

    void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const ARRAY<T,TV_INT>& density_ghost,const T dt,const T time)
    {
    }

//#####################################################################
};
}
#endif
