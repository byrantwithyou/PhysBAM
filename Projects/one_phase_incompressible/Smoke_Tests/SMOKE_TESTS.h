//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_TESTS__
#define __SMOKE_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VORTICITY_CONFINEMENT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>

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

    GRID<TV> upsampled_mac_grid;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;

public:
    using BASE::mac_grid; using BASE::incompressible;using BASE::projection;using BASE::output_directory;using BASE::mpi_grid;using BASE::domain_boundary;
    using BASE::face_velocities;using BASE::last_frame;using BASE::write_substeps_level;using BASE::rigid_body_collection;using BASE::boundary_scalar;using BASE::density;
    using BASE::boundary;using BASE::restart;using BASE::temperature;

    SMOKE_TESTS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :INCOMPRESSIBLE_EXAMPLE<TV>(stream_type_input),vc((T).06),test_number(1),stream_type(stream_type_input),use_collisions(false),upsample(1),buoyancy_clamp(1)
    {
        int scale=64;
        T source_radius=(T).05;
        parse_args.Add("-scale",&scale,"scale","grid resolution");
        parse_args.Add("-restart",&restart,"frame","restart");
        parse_args.Add("-vc",&vc,"scalar","vorticity confinement");
        parse_args.Add("-source_radius",&source_radius,"radius","radius of source");
        parse_args.Add("-buoyancy",&buoyancy_clamp,"const","buoyancy constant");
        parse_args.Add("-test_number",&test_number,"test","test number");
        parse_args.Add("-substeps",&write_substeps_level,"level","level of writing sub-steps");
        parse_args.Add("-upsample",&upsample,"level","level of refinement");
        parse_args.Parse();

        last_frame=200;
        output_directory=STRING_UTILITIES::string_sprintf("Smoke_Tests/Test_%d_%d_%d_vc_%1.2f_%dd_%1.2f",test_number,scale,upsample,vc,TV::dimension,buoyancy_clamp);
        if(TV::dimension==2){
            source_box_2d.min_corner=TV::Constant_Vector(0.25-source_radius);
            source_box_2d.max_corner=TV::Constant_Vector(0.25+source_radius);
            source_box_2d.min_corner.y=T(0);
            source_box_2d.max_corner.y=T(0.05);}
        else{
            source_box_3d.radius=T(source_radius);
            VECTOR<T,3> endpoint1=VECTOR<T,3>::Constant_Vector(0.25),endpoint2=VECTOR<T,3>::Constant_Vector(0.25);
            endpoint1(1)=T(0);endpoint2(1)=T(0.05);
            source_box_3d.Set_Endpoints(endpoint1,endpoint2);}
        TV_INT counts=TV_INT::All_Ones_Vector()*scale;
        RANGE<TV> range=RANGE<TV>(TV(),TV::Constant_Vector(0.5));
        counts(1)*=2;
        range.max_corner(1)*=2;
        mac_grid.Initialize(counts,range,true);
        upsampled_mac_grid.Initialize(counts*upsample,range,true);
        if(test_number>1){
            incompressible.Initialize_Grids(mac_grid);
            incompressible.Set_Body_Force(true);}
    }

    ~SMOKE_TESTS()
    {
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
                rigid_particle_id=rigid_body_collection.Add_Rigid_Body(stream_type,model_file_name,T(.07),true,true,true,true);
                rigid_body_collection.rigid_body_particles.frame(rigid_particle_id).t=TV::Constant_Vector(T(.25));
                rigid_body_collection.rigid_body_particles.frame(rigid_particle_id).t(1)=0.5;
                rigid_body_collection.rigid_body_particles.rigid_body(rigid_particle_id)->is_static=true;
                break;
            case 3:
                use_collisions=true;
                if(TV::dimension==2) model_file_name="../../Public_Data/Rigid_Bodies_2D/circle";
                else model_file_name="../../Public_Data/Rigid_Bodies/sphere";
                rigid_particle_id=rigid_body_collection.Add_Rigid_Body(stream_type,model_file_name,T(.1),true,true,true,true);
                rigid_body_collection.rigid_body_particles.rigid_body(rigid_particle_id)->is_static=false;
                break;}
        rigid_body_collection.Update_Kinematic_Particles();
    }
    
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
    {
        if(test_number==3){
            frame.t=TV::Constant_Vector(T(0.25));frame.t(1)=0.5;
            if(time<4.0) frame.t(1)=sin(2*pi*time/4.)*0.2+0.5;
            else frame.t(0)=sin(2*pi*time/4.)*0.1+0.25;}
    }

    void Write_Output_Files(const int frame)
    {
        BASE::Write_Output_Files(frame);
        if(upsample==1) return;
        ARRAY<T,TV_INT> upsampled_density(upsampled_mac_grid.Domain_Indices());
        ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(mac_grid,3,false);
        boundary->Fill_Ghost_Faces(mac_grid,face_velocities,face_velocities_ghost,0,3);
        for(CELL_ITERATOR<TV> iterator(upsampled_mac_grid);iterator.Valid();iterator.Next()){
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
    for(int axis=0;axis<TV::dimension;axis++) for(int axis_side=0;axis_side<2;axis_side++){int side=2*axis+axis_side;
        if(domain_boundary(axis)(axis_side)){
            TV_INT interior_cell_offset=axis_side==0?TV_INT():-TV_INT::Axis_Vector(axis);    
            for(FACE_ITERATOR<TV> iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                TV_INT boundary_face=axis_side==0?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(axis!=2){ if(face_velocities.Component(axis).Valid_Index(boundary_face)) projection.elliptic_solver->psi_N(FACE_INDEX<TV::dimension>(axis,boundary_face))=true;}
                else {projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
            for(FACE_ITERATOR<TV> iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT boundary_face=axis_side==0?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
                if(axis!=2 && face_velocities.Component(axis).Valid_Index(boundary_face)) face_velocities(FACE_INDEX<TV::dimension>(axis,boundary_face))=0;}}}
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        if(test_number==1 && Source_Box_Lazy_Inside(iterator.Location())){
            projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
            if(iterator.Axis()==1)face_velocities(iterator.Full_Index())=1;
            else face_velocities(iterator.Full_Index())=0;}
        if(use_collisions){
            int first_cell_in_solid=rigid_body_collection.rigid_body_particles.rigid_body(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(iterator.First_Cell_Center()),second_cell_in_solid=rigid_body_collection.rigid_body_particles.rigid_body(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(iterator.Second_Cell_Center());
            if(first_cell_in_solid+second_cell_in_solid==1 || (test_number==3 && first_cell_in_solid+second_cell_in_solid==2)){
                projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
                TV rigid_velocity=RIGID_BODY<TV>::Pointwise_Object_Velocity(rigid_body_collection.rigid_body_particles.twist(rigid_particle_id),rigid_body_collection.rigid_body_particles.frame(rigid_particle_id).t,iterator.Location());
                TV fluid_velocity;
                fluid_velocity(iterator.Axis())=face_velocities.Component(iterator.Axis())(iterator.Face_Index());
                for(int i=1;i<TV::dimension;i++){int axis=(iterator.Axis()-1+i)%TV::dimension+1;
                    fluid_velocity(axis)=(
                        face_velocities.Component(axis)(iterator.Face_Index())
                        +face_velocities.Component(axis)(iterator.Face_Index()-TV_INT::Axis_Vector(iterator.Axis()))
                        +face_velocities.Component(axis)(iterator.Face_Index()+TV_INT::Axis_Vector(axis))
                        +face_velocities.Component(axis)(iterator.Face_Index()+TV_INT::Axis_Vector(axis)-TV_INT::Axis_Vector(iterator.Axis())))*0.25;}
                TV rigid_normal=rigid_body_collection.Rigid_Body(rigid_particle_id).Implicit_Geometry_Normal(iterator.Location());
                T projected_velocity=Dot_Product(rigid_normal,fluid_velocity-rigid_velocity);
                if(projected_velocity<0) face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(fluid_velocity-rigid_normal*projected_velocity)(iterator.Axis());}
            else if(first_cell_in_solid+second_cell_in_solid==2 && test_number==2){
                projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
                face_velocities.Component(iterator.Axis())(iterator.Face_Index())=RIGID_BODY<TV>::Pointwise_Object_Velocity(rigid_body_collection.rigid_body_particles.twist(rigid_particle_id),rigid_body_collection.rigid_body_particles.frame(rigid_particle_id).t,iterator.Location())(iterator.Axis());}}}}

    void Initialize_Fields()
    {for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()) face_velocities(iterator.Full_Index())=0;
    for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()) density(iterator.Cell_Index())=temperature(iterator.Cell_Index())=0;}

    void Get_Scalar_Field_Sources(const T time)
    {for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next())
        if(Source_Box_Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=1;
    for(CELL_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next())
        if(use_collisions && rigid_body_collection.rigid_body_particles.rigid_body(rigid_particle_id)->Implicit_Geometry_Lazy_Inside(iterator.Location())) density(iterator.Cell_Index())=0;}
    
    void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::dimension> >& force,const ARRAY<T,TV_INT>& density_ghost,const T dt,const T time)
    {
        T density_buoyancy_constant=T(2)/buoyancy_clamp;
        for(FACE_ITERATOR<TV> iterator(mac_grid,0,GRID<TV>::WHOLE_REGION,-1,1);iterator.Valid();iterator.Next()){ // y-direction forces only
            T face_density=min((density_ghost(iterator.First_Cell_Index())+density_ghost(iterator.Second_Cell_Index()))*T(.5),buoyancy_clamp);
            T density_difference=face_density;
            if(density_difference>0) force.Component(1)(iterator.Face_Index())=density_buoyancy_constant*density_difference;}
    }

//#####################################################################
};
}
#endif
