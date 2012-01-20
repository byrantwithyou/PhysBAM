//#####################################################################
// Copyright 2010, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADVECTION_TESTS__
#define __ADVECTION_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_CONSERVATIVE_ENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>

namespace PhysBAM{

template<class TV>
class ADVECTION_TESTS:public INCOMPRESSIBLE_EXAMPLE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef INCOMPRESSIBLE_EXAMPLE<TV> BASE;

    bool use_conservative_advection,use_eno_advection;
    int test_number;
    int local_order;

public:
    using BASE::mac_grid; using BASE::incompressible;using BASE::projection;using BASE::output_directory;using BASE::incompressible;using BASE::mpi_grid;using BASE::domain_boundary;using BASE::face_velocities;
    using BASE::last_frame;using BASE::write_substeps_level;using BASE::rigid_geometry_collection;using BASE::boundary_scalar;using BASE::density;using BASE::boundary;using BASE::restart;using BASE::analytic_test;
    using BASE::cfl;using BASE::frame_rate;using BASE::stream_type;using BASE::temperature;

    ADVECTION_TESTS(const STREAM_TYPE stream_type_input,const PARSE_ARGS& parse_args)
        :INCOMPRESSIBLE_EXAMPLE<TV>(stream_type_input),use_conservative_advection(false)
    {
        last_frame=250;analytic_test=true;
        write_substeps_level=parse_args.Get_Integer_Value("-substep");;
        test_number=parse_args.Get_Integer_Value("-test_number");
        if(test_number==3) last_frame=24;
        int scale=parse_args.Get_Integer_Value("-scale");
        cfl=parse_args.Get_Double_Value("-cfl");
        restart=parse_args.Get_Integer_Value("-restart");
        use_conservative_advection=parse_args.Is_Value_Set("-conservative");
        use_eno_advection=parse_args.Is_Value_Set("-eno");
        local_order=parse_args.Get_Integer_Value("-order");
        output_directory=STRING_UTILITIES::string_sprintf("Advection_Tests/Test_%d_%d_%d_%1.2f%s%s",test_number,scale,TV::dimension,cfl,use_conservative_advection?"_conservative":(use_eno_advection?"_eno":""),local_order>1?"_high_order":"");
        TV_INT counts=TV_INT::All_Ones_Vector()*scale*(test_number==3?1:3);
        RANGE<TV> range=RANGE<TV>(TV(),TV::All_Ones_Vector()*((test_number==3)?100:15));
        mac_grid.Initialize(counts,range,true);
        if(use_conservative_advection){
            if(local_order==3){
                ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV>,FACE_LOOKUP_UNIFORM<GRID<TV> > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > >* advection_conservative=
                    new ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV>,FACE_LOOKUP_UNIFORM<GRID<TV> > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > >();
                advection_conservative->use_second_order=true;
                //advection_conservative->clamp_weights=false;
                //ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV>,FACE_LOOKUP_UNIFORM<GRID<TV> > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > >* advection_conservative=
                //    new ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV>,FACE_LOOKUP_UNIFORM<GRID<TV> > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > >();
                incompressible.Set_Custom_Advection(*advection_conservative);}
            else if(local_order==2){
                ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV>,FACE_LOOKUP_UNIFORM<GRID<TV> > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > >* advection_conservative=
                    new ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV>,FACE_LOOKUP_UNIFORM<GRID<TV> > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > >();
                advection_conservative->use_second_order=true;
                advection_conservative->clamp_weights=false;
                //ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV>,FACE_LOOKUP_UNIFORM<GRID<TV> > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > >* advection_conservative=
                //    new ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T,AVERAGING_UNIFORM<GRID<TV>,FACE_LOOKUP_UNIFORM<GRID<TV> > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > >();
                incompressible.Set_Custom_Advection(*advection_conservative);}
            else if(local_order==1){
                ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=new ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>();
                advection_conservative->clamp_weights=false;
                incompressible.Set_Custom_Advection(*advection_conservative);}
            else PHYSBAM_FATAL_ERROR();}
        else if(use_eno_advection){
            ADVECTION_CONSERVATIVE_ENO<GRID<TV>,T>* eno=new ADVECTION_CONSERVATIVE_ENO<GRID<TV>,T>();
            eno->Set_Order(local_order);
            incompressible.Set_Custom_Advection(*eno);}
        else if(local_order>3) PHYSBAM_FATAL_ERROR();
    }

    ~ADVECTION_TESTS()
    {
        if(use_conservative_advection || local_order>1) delete incompressible.advection;
    }

    void Initialize_Bodies()
    {
    }
    
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
    {
    }

    void Write_Output_Files(const int frame)
    {
        T total_density=0;
        BASE::Write_Output_Files(frame);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/time",output_directory.c_str(),frame),BASE::Time_At_Frame(frame));
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            total_density+=density(iterator.Cell_Index())*mac_grid.min_dX;}
        //for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
        //    if(frame==50 && density(iterator.Cell_Index())>.1&&density(iterator.Cell_Index())<.9) std::cout<<"Pos is "<<iterator.Location()<<std::endl;}
        ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(incompressible.advection);
        if(advection_conservative) total_density+=advection_conservative->total_mass_lost;
        if(advection_conservative) total_density-=advection_conservative->total_mass_gained;
        T max_error=0,total_error=0;int count=0;
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            TV X_o=iterator.Location()-((T)frame)*1/(T)frame_rate;
            T analytic=0;if(X_o.x<0.75 && X_o.x>0.25) analytic=1;
            T error=abs(analytic-(density(iterator.Cell_Index())));
            count++;max_error=max(error,max_error);total_error+=error;}
        LOG::cout<<std::endl;
        LOG::cout<<"Total density is "<<total_density<<std::endl;
        TV location=TV::All_Ones_Vector()*2.78;
        LOG::cout<<"Frame "<<frame<<" has error "<<total_error/count<<" max "<<max_error<<std::endl;
        if(frame==50) std::cout<<"Density is "<<density(mac_grid.Index(location))<<std::endl;
        if(frame==50) std::cout<<"Density Error is "<<(1-density(mac_grid.Index(location)))<<std::endl;
    }

    void Set_Boundary_Conditions(const T time)
    {
    }

    void Initialize_Fields()
    {
        if(test_number==1){    
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=1;
            for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
                //density(iterator.Cell_Index())=.5;
                if(iterator.Location().x<=0.75 && iterator.Location().x>=0.25) density(iterator.Cell_Index())=.5*(sin(2*pi/(.5)*(iterator.Location().x-0.25)-pi/2.)+1); else density(iterator.Cell_Index())=0;
                //if(iterator.Location().x<0.75 && iterator.Location().x>0.25) density(iterator.Cell_Index())=1.5;else density(iterator.Cell_Index())=.5;
                //if(iterator.Location().x<0.75 && iterator.Location().x>0.25) density(iterator.Cell_Index())=1;else density(iterator.Cell_Index())=0;
                temperature(iterator.Cell_Index())=0;}}
        else if(test_number==2){
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=sin(pi/(5.)*iterator.Location().x);
            //for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=1+.5*sin(2*pi/(5)*iterator.Location().x);
            for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
                if(iterator.Location().x<=0.75 && iterator.Location().x>=0.25) density(iterator.Cell_Index())=1;else density(iterator.Cell_Index())=0;
                temperature(iterator.Cell_Index())=0;}}
        else if(test_number==3){
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
                if(iterator.Axis()==1) face_velocities(iterator.Full_Index())=628*(pi/314.*(50.-iterator.Location()(2)));
                else if(iterator.Axis()==2) face_velocities(iterator.Full_Index())=628*(pi/314.*(iterator.Location()(1)-50.));}
            TV center=TV::All_Ones_Vector()*75;center(1)=50;
            SPHERE<TV> sphere(center,15);BOX<TV> box;box.min_corner=TV::All_Ones_Vector()*50;box.min_corner(1)=47.5;box.max_corner=TV::All_Ones_Vector()*75;box.max_corner(1)=52.5;
            for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
                if(sphere.Lazy_Inside(iterator.Location()) && !box.Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;else density(iterator.Cell_Index())=0;
                temperature(iterator.Cell_Index())=0;}}
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
