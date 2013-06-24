//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_TESTS__
#define __SMOKE_TESTS__

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <Dynamics/Drivers/PLS_EXAMPLE.h>

namespace PhysBAM{

template<class TV>
class WATER_TESTS:public PLS_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef PLS_EXAMPLE<TV> BASE;

    CYLINDER<T> source_cyl;
    RANGE<TV> source;

public:
    using BASE::mac_grid; using BASE::incompressible;using BASE::projection;using BASE::output_directory;using BASE::mpi_grid;using BASE::domain_boundary;using BASE::face_velocities;
    using BASE::particle_levelset_evolution;using BASE::write_substeps_level;using BASE::restart;using BASE::last_frame;

    WATER_TESTS(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args)
        :PLS_EXAMPLE<TV>(stream_type)
    {
        int test_number=1;
        last_frame=200;
        int scale=128;

        parse_args.Add("-restart",&restart,"frame","restart frame");
        parse_args.Add("-scale",&scale,"scale","fine scale grid resolution");
        parse_args.Add("-substeps",&write_substeps_level,"frame","output-substep level");
        parse_args.Parse();

        mac_grid.Initialize(TV_INT::All_Ones_Vector()*scale,RANGE<TV>(TV(),TV::All_Ones_Vector()),true);
        output_directory=STRING_UTILITIES::string_sprintf("Water_Tests/Test_%d_%d",test_number,scale);
        if(TV::dimension==3){
            VECTOR<T,3> point1,point2;
            point1=VECTOR<T,3>::All_Ones_Vector()*(T).6;point1(0)=.4;point1(1)=.95;point2=VECTOR<T,3>::All_Ones_Vector()*(T).6;point2(0)=.4;point2(1)=1;
            source_cyl.Set_Endpoints(point1,point2);source_cyl.radius=.1;}
        else{
            TV point1,point2;
            point1=TV::All_Ones_Vector()*(T).5;point1(0)=.4;point1(1)=.95;point2=TV::All_Ones_Vector()*(T).65;point2(0)=.55;point2(1)=1;
            source.min_corner=point1;source.max_corner=point2;}
    }

    void Write_Output_Files(const int frame)
    {BASE::Write_Output_Files(frame);}

    void Set_Boundary_Conditions(const T time)
    {projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*axis+axis_side;
        TV_INT interior_cell_offset=axis_side==0?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==0?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==0?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(domain_boundary(axis)(axis_side)){
            for(FACE_ITERATOR<TV> iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                if(particle_levelset_evolution.phi(face+interior_cell_offset)<=0){
                    if(face_velocities.Component(axis).Valid_Index(face)){projection.elliptic_solver->psi_N.Component(axis)(face)=true;face_velocities.Component(axis)(face)=0;}}
                else{TV_INT cell=face+exterior_cell_offset;projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}}
        else for(FACE_ITERATOR<TV> iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
            projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        if(time<=3 && Lazy_Inside_Source(iterator.Location())){
            projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==1) face_velocities(iterator.Full_Index())=-1;
            else face_velocities(iterator.Full_Index())=0;}}}

    bool Lazy_Inside_Source(const VECTOR<T,2> X)
    {
        return source.Lazy_Inside(X);
    }

    bool Lazy_Inside_Source(const VECTOR<T,3> X)
    {
        return source_cyl.Lazy_Inside(X);
    }

    void Set_Phi_Inside_Source(const VECTOR<int,2>& index,const VECTOR<T,2> X)
    {
        if(source.Lazy_Inside(X)) particle_levelset_evolution.phi(index)=min(particle_levelset_evolution.phi(index),source.Signed_Distance(X));
    }

    void Set_Phi_Inside_Source(const VECTOR<int,3>& index,const VECTOR<T,3> X)
    {
        if(source_cyl.Lazy_Inside(X)) particle_levelset_evolution.phi(index)=min(particle_levelset_evolution.phi(index),source_cyl.Signed_Distance(X));
    }

    void Adjust_Phi_With_Sources(const T time)
    {
        if(time>3) return;
        for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){Set_Phi_Inside_Source(iterator.Cell_Index(),iterator.Location());}
    }

    void Initialize_Phi()
    {
        ARRAY<T,TV_INT>& phi=particle_levelset_evolution.phi;
        for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
            //TV center=TV::All_Ones_Vector()*.5;center(1)=.75;
            //static SPHERE<TV> circle(center,(T).2);
            const TV &X=iterator.Location();
            //phi(iterator.Cell_Index())=min(circle.Signed_Distance(X),X.y-(T).412134);}
            phi(iterator.Cell_Index())=X.y-(T)mac_grid.dX.Min()*5;}
            //phi(iterator.Cell_Index())=X.y-(T).412134;}
    }

//#####################################################################
};
}

#endif
