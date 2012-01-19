//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
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

    CYLINDER<T> source_box_3d;
    RANGE<TV> source_box_2d;
    T vc;
    int test_number;
    STREAM_TYPE stream_type;
    int rigid_particle_id;
    bool use_collisions;
    int upsample;
    T buoyancy_clamp;
    bool use_conservative_advection;

    GRID<TV> upsampled_mac_grid;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;

public:
    using BASE::mac_grid; using BASE::incompressible;using BASE::projection;using BASE::output_directory;using BASE::incompressible;using BASE::mpi_grid;using BASE::domain_boundary;
    using BASE::face_velocities;using BASE::last_frame;using BASE::write_substeps_level;using BASE::rigid_geometry_collection;using BASE::boundary_scalar;using BASE::density;
    using BASE::boundary;using BASE::restart;using BASE::temperature;

    SMOKE_TESTS(const STREAM_TYPE stream_type_input,const PARSE_ARGS& parse_args)
        :INCOMPRESSIBLE_EXAMPLE<TV>(stream_type_input),stream_type(stream_type_input),use_collisions(false),use_conservative_advection(false)
    {
        last_frame=200;
        write_substeps_level=parse_args.Get_Integer_Value("-substep");;
        test_number=parse_args.Get_Integer_Value("-test_number");
        int scale=parse_args.Get_Integer_Value("-scale");
        upsample=parse_args.Get_Integer_Value("-upsample");
        restart=parse_args.Get_Integer_Value("-restart");
        vc=parse_args.Get_Double_Value("-vc");
        buoyancy_clamp=parse_args.Get_Double_Value("-buoyancy");
        output_directory=STRING_UTILITIES::string_sprintf("Smoke_Tests/Test_%d_%d_%d_vc_%1.2f_%dd_%1.2f",test_number,scale,upsample,vc,TV::dimension,buoyancy_clamp);
        T source_radius=parse_args.Get_Double_Value("-source_radius");
        use_conservative_advection=parse_args.Is_Value_Set("-conservative");
        if(TV::dimension==2){
            source_box_2d.min_corner=TV::Constant_Vector(0.25-source_radius);
            source_box_2d.max_corner=TV::Constant_Vector(0.25+source_radius);
            source_box_2d.min_corner.y=T(0);
            source_box_2d.max_corner.y=T(0.05);}
        else{
            source_box_3d.radius=T(source_radius);
            VECTOR<T,3> endpoint1=VECTOR<T,3>::Constant_Vector(0.25),endpoint2=VECTOR<T,3>::Constant_Vector(0.25);
            endpoint1(2)=T(0);endpoint2(2)=T(0.05);
            source_box_3d.Set_Endpoints(endpoint1,endpoint2);}
        TV_INT counts=TV_INT::All_Ones_Vector()*scale;
        RANGE<TV> range=RANGE<TV>(TV(),TV::Constant_Vector(0.5));
        counts(2)*=2;
        range.max_corner(2)*=2;
        mac_grid.Initialize(counts,range,true);
        upsampled_mac_grid.Initialize(counts*upsample,range,true);
        if(test_number>1){
            incompressible.Initialize_Grids(mac_grid);
            incompressible.Set_Body_Force(true);}
        if(use_conservative_advection) incompressible.Set_Custom_Advection(*new ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>());
    }

    ~SMOKE_TESTS()
    {
        if(use_conservative_advection) delete incompressible.advection;
    }

    void Initialize_Confinement()
    {
        incompressible.Set_Vorticity_Confinement(vc);
    }

    void Initialize_Bodies()
    {
        std::string model_file_name;
        switch(test_number){
            case 1:
                break;
            case 2:
                use_collisions=true;
                if(TV::dimension==2) model_file_name="../../Public_Data/Rigid_Bodies_2D/circle";
                else model_file_name="../../Public_Data/Rigid_Bodies/sphere";
                rigid_particle_id=rigid_geometry_collection.Add_Rigid_Geometry(stream_type,model_file_name,T(.07),true,true,true,true);
                rigid_geometry_collection.particles.X(rigid_particle_id)=TV::Constant_Vector(T(.25));
                rigid_geometry_collection.particles.X(rigid_particle_id)(2)=0.5;
                rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->is_static=true;
                break;
            case 3:
                use_collisions=true;
                if(TV::dimension==2) model_file_name="../../Public_Data/Rigid_Bodies_2D/circle";
                else model_file_name="../../Public_Data/Rigid_Bodies/sphere";
                rigid_particle_id=rigid_geometry_collection.Add_Rigid_Geometry(stream_type,model_file_name,T(.1),true,true,true,true);
                rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->is_static=false;
                break;}
        rigid_geometry_collection.Update_Kinematic_Particles();
    }
    
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
    {
        if(test_number==3){
            frame.t=TV::Constant_Vector(T(0.25));frame.t(2)=0.5;
            if(time<4.0) frame.t(2)=sin(2*pi*time/4.)*0.2+0.5;
            else frame.t(1)=sin(2*pi*time/4.)*0.1+0.25;}
    }

    void Write_Output_Files(const int frame)
    {
        BASE::Write_Output_Files(frame);
        if(upsample==1) return;
        ARRAY<T,TV_INT> upsampled_density(upsampled_mac_grid.Domain_Indices());
        ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(mac_grid,3,false);
        boundary->Fill_Ghost_Cells_Face(mac_grid,face_velocities,face_velocities_ghost,0,3);
        for(typename GRID<TV>::CELL_ITERATOR iterator(upsampled_mac_grid);iterator.Valid();iterator.Next()){
            upsampled_density(iterator.Cell_Index())=interpolation.Clamped_To_Array(mac_grid,density,iterator.Location());}
        std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_mac_velocities",face_velocities_ghost);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",upsampled_mac_grid);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/coarse_grid",mac_grid);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",upsampled_density);
    }

    template<class TV2>
    bool Source_Box_Lazy_Inside(const TV2& location){PHYSBAM_NOT_IMPLEMENTED();return false;}

    template<class T2>
    bool Source_Box_Lazy_Inside(const VECTOR<T2,2>& location)
    {return source_box_2d.Lazy_Inside(location);}

    template<class T2>
    bool Source_Box_Lazy_Inside(const VECTOR<T2,3>& location)
    {return source_box_3d.Lazy_Inside(location);}

    void Set_Boundary_Conditions(const T time)
    {projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*(axis-1)+axis_side;
        if(domain_boundary(axis)(axis_side)){
            TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);    
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                TV_INT boundary_face=axis_side==1?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(axis!=2){ if(face_velocities.Component(axis).Valid_Index(boundary_face)) projection.elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(axis,boundary_face))=true;}
                else {projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT boundary_face=axis_side==1?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(axis!=2 && face_velocities.Component(axis).Valid_Index(boundary_face)) face_velocities(FACE_INDEX<TV::dimension>(axis,boundary_face))=0;}}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
        if(test_number==1 && Source_Box_Lazy_Inside(iterator.Location())){
            projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==2)face_velocities(iterator.Full_Index())=1;
            else face_velocities(iterator.Full_Index())=0;}
        if(use_collisions){
            int first_cell_in_solid=rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(iterator.First_Cell_Center()),second_cell_in_solid=rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(iterator.Second_Cell_Center());
            if(first_cell_in_solid+second_cell_in_solid==1 || (test_number==3 && first_cell_in_solid+second_cell_in_solid==2)){
                projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
                TV rigid_velocity=RIGID_GEOMETRY<TV>::Pointwise_Object_Velocity(rigid_geometry_collection.particles.twist(rigid_particle_id),rigid_geometry_collection.particles.X(rigid_particle_id),iterator.Location());
                TV fluid_velocity;
                fluid_velocity(iterator.Axis())=face_velocities.Component(iterator.Axis())(iterator.Face_Index());
                for(int i=1;i<TV::dimension;i++){int axis=(iterator.Axis()-1+i)%TV::dimension+1;
                    fluid_velocity(axis)=(
                        face_velocities.Component(axis)(iterator.Face_Index())
                        +face_velocities.Component(axis)(iterator.Face_Index()-TV_INT::Axis_Vector(iterator.Axis()))
                        +face_velocities.Component(axis)(iterator.Face_Index()+TV_INT::Axis_Vector(axis))
                        +face_velocities.Component(axis)(iterator.Face_Index()+TV_INT::Axis_Vector(axis)-TV_INT::Axis_Vector(iterator.Axis())))*0.25;}
                TV rigid_normal=rigid_geometry_collection.Rigid_Geometry(rigid_particle_id).Implicit_Geometry_Normal(iterator.Location());
                T projected_velocity=Dot_Product(rigid_normal,fluid_velocity-rigid_velocity);
                if(projected_velocity<0) face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(fluid_velocity-rigid_normal*projected_velocity)(iterator.Axis());}
            else if(first_cell_in_solid+second_cell_in_solid==2 && test_number==2){
                projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
                face_velocities.Component(iterator.Axis())(iterator.Face_Index())=RIGID_GEOMETRY<TV>::Pointwise_Object_Velocity(rigid_geometry_collection.particles.twist(rigid_particle_id),rigid_geometry_collection.particles.X(rigid_particle_id),iterator.Location())(iterator.Axis());}}}}

    void Initialize_Fields()
    {for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=temperature(iterator.Cell_Index())=0;}

    void Get_Scalar_Field_Sources(const T time)
    {for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
        if(Source_Box_Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
        if(use_collisions && rigid_geometry_collection.particles.rigid_geometry(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=0;}
    
    void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const ARRAY<T,TV_INT>& density_ghost,const T dt,const T time)
    {
        T density_buoyancy_constant=T(2)/buoyancy_clamp;
        for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,0,GRID<TV>::WHOLE_REGION,0,2);iterator.Valid();iterator.Next()){ // y-direction forces only
            T face_density=min((density_ghost(iterator.First_Cell_Index())+density_ghost(iterator.Second_Cell_Index()))*T(.5),buoyancy_clamp);
            T density_difference=face_density;
            if(density_difference>0) force.Component(2)(iterator.Face_Index())=density_buoyancy_constant*density_difference;}
    }

//#####################################################################
};
}
#endif
