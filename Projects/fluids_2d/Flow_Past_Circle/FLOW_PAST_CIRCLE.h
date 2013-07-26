//#####################################################################
// Copyright 2001-2006, Ron Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOW_PAST_CIRCLE
//#####################################################################
#ifndef __FLOW_PAST_CIRCLE__
#define __FLOW_PAST_CIRCLE__

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template <class T>
class FLOW_PAST_CIRCLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T,2> >
{
public:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;
    using BASE::output_directory;using BASE::data_directory;using BASE::stream_type;using BASE::fluid_collection;using BASE::solid_body_collection;using BASE::resolution;
    using BASE::Mark_Outside;using BASE::Get_Boundary_Along_Ray;using BASE::parse_args;

    RIGIDS_STANDARD_TESTS<TV> solids_tests;
    T rho,rho_bottom,rho_top,buoyancy_constant;
    RANGE<TV> source_domain;
    ARRAY<T,TV_INT> phi_object;
    LEVELSET<TV> levelset_object;
    SPHERE<TV> circle;
    GEOMETRY_PARTICLES<TV> debug_particles;
    ARRAY<TV> sample_points;
    bool shed,opt_enlarge;
    
    FLOW_PAST_CIRCLE(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type,0,fluids_parameters.SMOKE),solids_tests(stream_type,data_directory,solid_body_collection.rigid_body_collection),
        levelset_object(*fluids_parameters.grid,phi_object),circle(TV((T)2,(T)2),(T).5),shed(false),opt_enlarge(false)
    {
        //fluids_parameters.cfl=0.75;
        fluids_parameters.cfl=.9;
        fluids_parameters.gravity=0;
        fluids_parameters.density=1;
        //fluids_parameters.cfl=1.75;
        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_levelset_viscosity=true;
        write_output_files=true;fluids_parameters.write_debug_data=true;
        source_domain=RANGE<TV>(TV((T).45,(T)0),TV((T).55,(T).1));
//        fluids_parameters.use_maccormack_semi_lagrangian_advection=true;

        sample_points.Append(TV(1.25,2.5));
        sample_points.Append(TV(3,2.25));
        sample_points.Append(TV(2,3));

        output_directory="Flow_Past_Circle/output";
        debug_particles.template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    }
    
    ~FLOW_PAST_CIRCLE()
    {}

    T Constraints_CFL() PHYSBAM_OVERRIDE {return FLT_MAX;}
    bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE {return false;}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_SPH_Particles_For_Sources(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Adjust_Phi_With_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_Begin_Callback(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_End_Callback(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {} 
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Initialize_SPH_Particles() PHYSBAM_OVERRIDE {}
    void Initialize_Velocities() PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {} // time at start of substep
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Read_Output_Files_Fluids(const int frame) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    void Set_PD_Targets(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Setup_Initial_Refinement() PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add("-viscosity",&fluids_parameters.viscosity,"value","viscosity");
    parse_args->Add("-enlarge",&opt_enlarge,"value","Enlarge");
    parse_args->Add("-shed",&shed,"shed");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    solids_tests.data_directory=data_directory;
    if(shed) fluids_parameters.grid->Initialize(TV_INT((int)(2.5*resolution)+1,resolution+1),RANGE<TV>(TV(-(T)2.5,-(T)3.5),TV(15,(T)3.5)));
    else fluids_parameters.grid->Initialize(TV_INT(resolution+1,resolution+1),RANGE<TV>(TV(0,0),TV(4,4)));
    if(fluids_parameters.viscosity) fluids_parameters.implicit_viscosity=true;
    if(opt_enlarge) circle.radius+=fluids_parameters.grid->dX.Min();
    if(shed) fluids_parameters.viscosity=(T).01;
    else fluids_parameters.viscosity=(T).1;
    if(shed) circle.center=TV();
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Analytic_Sphere(circle.radius,fluids_parameters.density,7);
    rigid_body.Frame().t=circle.center;
    rigid_body.is_static=true;
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(solid_body_collection.rigid_body_collection);
    fluids_parameters.incompressible_iterations=1000;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<2> >& face_velocities,ARRAY<bool,FACE_INDEX<2> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV> u_grid=fluids_parameters.grid->Get_Face_Grid(0),v_grid=fluids_parameters.grid->Get_Face_Grid(1);
    T flow_speed=(T)1;
    for(int j=0;j<u_grid.counts.y;j++){
        FACE_INDEX<2> a(1,TV_INT(1,j)),b(1,TV_INT(u_grid.counts.x,j));
        face_velocities(a)=flow_speed;
        psi_N(a)=true;}
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Construct_Levelsets_For_Objects(time);
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time)
{
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        if(circle.Inside(iterator.Location(),(T).1*fluids_parameters.grid->dX.Max()))
            fluids_parameters.incompressible->projection.elliptic_solver->psi_D(iterator.index)=true;
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(debug_particles.Size()){
        FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),frame));
        FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),frame),debug_particles);
        debug_particles.Delete_All_Elements();}
    if(frame==1){
        FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),0));
        FILE_UTILITIES::Write_To_File(this->stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),0),debug_particles);}

    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    LINEAR_INTERPOLATION_UNIFORM<TV,T> interp;
    for(int i=0;i<sample_points.m;i++){
        TV X=sample_points(i),V;
        for(int d=0;d<V.m;d++)
            V(d)=interp.Clamped_To_Array(fluids_parameters.grid->Get_Face_Grid(d),face_velocities.Component(d),X);
        LOG::cout<<"velocity at "<<X<<" : "<<V<<std::endl;}
}
//#####################################################################
// Function Mark_Outside
//#####################################################################
void Mark_Outside(ARRAY<bool,FACE_INDEX<TV::m> >& outside) PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,0);iterator.Valid();iterator.Next()) outside(iterator.Full_Index())=true;
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,2);iterator.Valid();iterator.Next()) outside(iterator.Full_Index())=true;
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,0,GRID<TV>::BOUNDARY_REGION,3);iterator.Valid();iterator.Next()) outside(iterator.Full_Index())=true;
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(circle.Lazy_Inside(iterator.Location())){
        VECTOR<FACE_INDEX<2>,4> faces;
        GRID<TV>::Neighboring_Faces(faces,iterator.index);
        for(int i=0;i<faces.m;i++) outside(faces(i))=true;}
}
//#####################################################################
// Function Get_Viscosity_Boundary_Along_Ray
//#####################################################################
typename BOUNDARY_CONDITIONS_CALLBACKS<TV>::RAY_TYPE Get_Boundary_Along_Ray(const FACE_INDEX<TV::m>& f1,const FACE_INDEX<TV::m>& f2,T& theta,T& value) PHYSBAM_OVERRIDE
{
    typename BOUNDARY_CONDITIONS_CALLBACKS<TV>::RAY_TYPE type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::unused;
    TV X1=fluids_parameters.grid->Face(f1);
    TV X2=fluids_parameters.grid->Face(f2);
    if(circle.Inside(X2,-fluids_parameters.grid->dX.Max()*(T)2)){ // circle
        theta=1;
        value=0;
        type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::noslip;}
    else if(f1.index.y!=f2.index.y){ // top or bottom
        theta=1;
        value=0;
        type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::slip;}
    else if(f1.index.x<f2.index.x){ // right
        theta=(T).5;
        value=0;
        type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::free;}
    else if(f1.index.x>f2.index.x){ // left
        theta=1;
        value=TV(1,0)(f2.axis);
        type=BOUNDARY_CONDITIONS_CALLBACKS<TV>::slip;}
    else PHYSBAM_FATAL_ERROR("Did not expect a boundary condition here.");

    static VECTOR<T,3> color_map[]={VECTOR<T,3>(1,0,0),VECTOR<T,3>(1,.5,0),VECTOR<T,3>(1,0,1),VECTOR<T,3>(0,.5,0),VECTOR<T,3>(0,1,1),VECTOR<T,3>(1,1,0)};

    if(ARRAY_VIEW<VECTOR<T,3> >* color_attribute=debug_particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR)){
        int p=debug_particles.Add_Element();
        debug_particles.X(p)=X1+theta*(X2-X1);
        (*color_attribute)(p)=color_map[type];}

    return type;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Write_Output_Files(frame);
    ARRAY<T,FACE_INDEX<2> > face_velocities_ghost(*fluids_parameters.grid,3,false);
    fluids_parameters.incompressible->boundary->Fill_Ghost_Faces(*fluids_parameters.grid,fluid_collection.incompressible_fluid_collection.face_velocities,face_velocities_ghost,0,3);
    ARRAY<VECTOR<T,1>,TV_INT> grid_vorticity(fluids_parameters.grid->Domain_Indices(3),false);
    ARRAY<T,TV_INT> grid_vorticity_magnitude(fluids_parameters.grid->Domain_Indices(3),false);
    VORTICITY_UNIFORM<TV>::Vorticity(*fluids_parameters.grid,FACE_LOOKUP_UNIFORM<TV>(face_velocities_ghost),grid_vorticity,grid_vorticity_magnitude);
    //CELL_ITERATOR<TV> fuckyou(*fluids_parameters.grid,3,GRID<TV>::GHOST_REGION);
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid,3,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()) grid_vorticity(iterator.Cell_Index())=VECTOR<T,1>();
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/grid_vorticity.%d",output_directory.c_str(),frame),grid_vorticity);
}
//#####################################################################
};
}
#endif


