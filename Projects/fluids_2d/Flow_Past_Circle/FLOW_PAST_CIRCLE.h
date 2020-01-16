//#####################################################################
// Copyright 2001-2006, Ron Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOW_PAST_CIRCLE
//#####################################################################
#ifndef __FLOW_PAST_CIRCLE__
#define __FLOW_PAST_CIRCLE__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Computations/VORTICITY_UNIFORM.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Poisson/PROJECTION_UNIFORM.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Projection/BOUNDARY_CONDITION_DOUBLE_FINE.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{

template <class T>
class FLOW_PAST_CIRCLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T,2> >
{
public:
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;
    using BASE::viewer_dir;using BASE::data_directory;using BASE::stream_type;using BASE::fluid_collection;using BASE::solid_body_collection;using BASE::resolution;
    using BASE::Mark_Outside;using BASE::Get_Boundary_Along_Ray;

    RIGIDS_STANDARD_TESTS<TV> solids_tests;
    T rho,rho_bottom,rho_top,buoyancy_constant;
    RANGE<TV> source_domain;
    SPHERE<TV> circle;
    IMPLICIT_OBJECT<TV> * io;
    ARRAY<TV> sample_points;
    SEGMENTED_CURVE_2D<T>* sc=0;
    
    bool shed,opt_enlarge,use_sc=false;
    
    FLOW_PAST_CIRCLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,0,fluids_parameters.SMOKE),solids_tests(stream_type_input,data_directory,solid_body_collection.rigid_body_collection),
        circle(TV((T)2,(T)2),(T).5),shed(false),opt_enlarge(false)
    {
        //fluids_parameters.cfl=0.75;
        fluids_parameters.cfl=.9;
        fluids_parameters.gravity=TV();
        fluids_parameters.density=1;
        //fluids_parameters.cfl=1.75;
        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=false;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_levelset_viscosity=false;
        write_output_files=true;fluids_parameters.write_debug_data=true;
        source_domain=RANGE<TV>(TV((T).45,(T)0),TV((T).55,(T).1));
//        fluids_parameters.use_maccormack_semi_lagrangian_advection=true;

        sample_points.Append(TV(1.25,2.5));
        sample_points.Append(TV(3,2.25));
        sample_points.Append(TV(2,3));

        if(!this->user_output_directory)
            viewer_dir.output_directory="Flow_Past_Circle/output";

        int sc_n=20;
        parse_args.Add("-viscosity",&fluids_parameters.viscosity,&fluids_parameters.implicit_viscosity,"value","viscosity");
        parse_args.Add("-enlarge",&opt_enlarge,"value","Enlarge");
        parse_args.Add("-shed",&shed,"shed");
        parse_args.Add("-sc",&use_sc,"use segmented curve");
        parse_args.Add("-sc_n",&sc_n,"n","segments on curve");
        parse_args.Parse();

        if(opt_enlarge) circle.radius+=fluids_parameters.grid->dX.Min();
        if(!fluids_parameters.implicit_viscosity)
        {
            if(shed) fluids_parameters.viscosity=(T).01;
            else fluids_parameters.viscosity=(T).1;
            fluids_parameters.implicit_viscosity=true;
        }
        if(shed) circle.center=TV();
        solids_tests.data_directory=data_directory;
        if(shed) fluids_parameters.grid->Initialize(TV_INT((int)(2.5*resolution),resolution),RANGE<TV>(TV(-(T)2.5,-(T)3.5),TV(15,(T)3.5)),true);
        else fluids_parameters.grid->Initialize(TV_INT(resolution,resolution),RANGE<TV>(TV(0,0),TV(4,4)),true);
        fluids_parameters.Use_Unified_Boundary_Conditions();
        io=Make_IO(circle);
        sc=TESSELLATION::Tessellate_Boundary(circle,sc_n);

        this->get_unified_boundary_conditions=[this](BOUNDARY_CONDITION_DOUBLE_FINE<TV>* bc_fine,const T time)
            {
                if(use_sc) bc_fine->Set(*sc,bc_fine->bc_noslip,[](const auto& data){Add_Debug_Particle(data.X,VECTOR<T,3>(1,0,0));return 0;});
                else bc_fine->Set(io,bc_fine->bc_noslip,0);
                bc_fine->Set_Domain_Walls(1,bc_fine->bc_noslip,[](const auto& data){return TV(1,0)(data.face.axis);});
                bc_fine->Set_Domain_Walls(2,bc_fine->bc_free,0);
                bc_fine->Set_Domain_Walls(12,bc_fine->bc_slip,0);
            };
    }
    
    ~FLOW_PAST_CIRCLE()
    {
        delete io;
        delete sc;
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    RIGID_BODY<TV>& rigid_body=solids_tests.Add_Analytic_Sphere(circle.radius,fluids_parameters.density,7);
    rigid_body.Frame().t=circle.center;
    rigid_body.is_static=true;
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(solid_body_collection.rigid_body_collection);
    fluids_parameters.incompressible_iterations=1000;
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
void Initialize_Advection()    override
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) override
{
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    LINEAR_INTERPOLATION_UNIFORM<TV,T> interp;
    for(int i=0;i<sample_points.m;i++){
        TV X=sample_points(i),V;
        for(int d=0;d<V.m;d++)
            V(d)=interp.Clamped_To_Array(fluids_parameters.grid->Get_Face_Grid(d),face_velocities.Component(d),X);
        LOG::cout<<"velocity at "<<X<<" : "<<V<<std::endl;}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files() const override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Write_Output_Files();
    ARRAY<T,FACE_INDEX<2> > face_velocities_ghost(*fluids_parameters.grid,3,no_init);
    fluids_parameters.incompressible->boundary->Fill_Ghost_Faces(*fluids_parameters.grid,fluid_collection.incompressible_fluid_collection.face_velocities,face_velocities_ghost,0,3);
    ARRAY<VECTOR<T,1>,TV_INT> grid_vorticity(fluids_parameters.grid->Domain_Indices(3),no_init);
    ARRAY<T,TV_INT> grid_vorticity_magnitude(fluids_parameters.grid->Domain_Indices(3),no_init);
    VORTICITY_UNIFORM<TV>::Vorticity(*fluids_parameters.grid,FACE_LOOKUP_UNIFORM<TV>(face_velocities_ghost),grid_vorticity,grid_vorticity_magnitude);
    //CELL_ITERATOR<TV> fuckyou(*fluids_parameters.grid,3,GRID<TV>::GHOST_REGION);
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid,3,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()) grid_vorticity(iterator.Cell_Index())=VECTOR<T,1>();
    Write_To_File(stream_type,viewer_dir.current_directory+"/grid_vorticity",grid_vorticity);
    if(use_sc)
    {
        Dump_Surface(*sc,VECTOR<T,3>(.5,1,0));
    }
}
//#####################################################################
};
}
#endif


