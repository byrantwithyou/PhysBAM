//#####################################################################
// Copyright 2009-2010, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_TESTS__
#define __SMOKE_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>

namespace PhysBAM{

template<class TV>
class SMOKE_TESTS:public INCOMPRESSIBLE_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef INCOMPRESSIBLE_EXAMPLE<TV> BASE;

    RANGE<TV> source_box;
    int test_number;
    bool use_conservative_advection;

    GRID<TV> upsampled_mac_grid;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;

public:
    using BASE::mac_grid; using BASE::incompressible;using BASE::projection;using BASE::output_directory;using BASE::incompressible;using BASE::mpi_grid;using BASE::domain_boundary;using BASE::face_velocities;
    using BASE::last_frame;using BASE::write_substeps_level;using BASE::rigid_geometry_collection;using BASE::boundary_scalar;using BASE::density;using BASE::boundary;using BASE::restart;

    SMOKE_TESTS(const STREAM_TYPE stream_type_input,const PARSE_ARGS& parse_args)
        :INCOMPRESSIBLE_EXAMPLE<TV>(stream_type_input),use_conservative_advection(false)
    {
        last_frame=200;
        write_substeps_level=parse_args.Get_Integer_Value("-substep");;
        test_number=parse_args.Get_Integer_Value("-test_number");
        int scale=parse_args.Get_Integer_Value("-scale");
        restart=parse_args.Get_Integer_Value("-restart");
        use_conservative_advection=parse_args.Is_Value_Set("-conservative");
        incompressible.conserve_kinetic_energy=parse_args.Is_Value_Set("-energy");
        output_directory=STRING_UTILITIES::string_sprintf("Smoke_Tests/Test_%d_%d_%d%s",test_number,scale,TV::dimension,use_conservative_advection?"_conservative":"");
        source_box.min_corner=TV::Constant_Vector(0.45);
        source_box.max_corner=TV::Constant_Vector(0.55);
        source_box.min_corner(2)=T(0);
        source_box.max_corner(2)=T(0.05);
        TV_INT counts=TV_INT::All_Ones_Vector()*scale;
        RANGE<TV> range=RANGE<TV>(TV(),TV::All_Ones_Vector());
        mac_grid.Initialize(counts,range,true);
        if(use_conservative_advection){
            incompressible.Set_Custom_Advection(*new ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>());}
    }

    ~SMOKE_TESTS()
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
        BASE::Write_Output_Files(frame);
        T total_energy=0;
        for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            if(!density.Valid_Index(iterator.First_Cell_Index()))
                total_energy+=.5*density(iterator.Second_Cell_Index())*face_velocities(iterator.Full_Index())*face_velocities(iterator.Full_Index())*mac_grid.dX(iterator.Axis());
            else if(!density.Valid_Index(iterator.Second_Cell_Index()))
                total_energy+=.5*density(iterator.First_Cell_Index())*face_velocities(iterator.Full_Index())*face_velocities(iterator.Full_Index())*mac_grid.dX(iterator.Axis());
            else
                total_energy+=.25*(density(iterator.First_Cell_Index())+density(iterator.Second_Cell_Index()))*face_velocities(iterator.Full_Index())*face_velocities(iterator.Full_Index())*mac_grid.dX(iterator.Axis());}
        std::cout<<std::endl<<"Total Energy is "<<total_energy<<std::endl;
    }

    template<class TV2>
    bool Source_Box_Lazy_Inside(const TV2& location){PHYSBAM_NOT_IMPLEMENTED();return false;}

    template<class T2>
    bool Source_Box_Lazy_Inside(const VECTOR<T2,2>& location)
    {return source_box.Lazy_Inside(location);}

    template<class T2>
    bool Source_Box_Lazy_Inside(const VECTOR<T2,3>& location)
    {return source_box.Lazy_Inside(location);}

    void Set_Boundary_Conditions(const T time)
    {projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
        if(domain_boundary(axis)(axis_side)){
            TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);    
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                TV_INT boundary_face=axis_side==1?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(axis!=2){ if(face_velocities.Component(axis).Valid_Index(boundary_face)) projection.elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(axis,boundary_face))=true;}
                else {projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT boundary_face=axis_side==1?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if((axis!=2) && face_velocities.Component(axis).Valid_Index(boundary_face)) face_velocities(FACE_INDEX<TV::dimension>(axis,boundary_face))=0;}}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
        if(Source_Box_Lazy_Inside(iterator.Location())){
            projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==2) face_velocities(iterator.Full_Index())=1;
            else face_velocities(iterator.Full_Index())=0;}}}

    void Initialize_Fields()
    {for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=temperature(iterator.Cell_Index())=0;}

    void Get_Scalar_Field_Sources(const T time)
    {for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
        if(Source_Box_Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;}
    
    void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const ARRAY<T,TV_INT>& density_ghost,const T dt,const T time)
    {
    }

//#####################################################################
};
}
#endif
