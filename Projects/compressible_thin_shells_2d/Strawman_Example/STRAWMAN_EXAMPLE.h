#ifndef __STRAWMAN_EXAMPLE__
#define __STRAWMAN_EXAMPLE__

#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_CONSERVATIVE_ENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/TYPED_STREAM.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_LINEAR_EXTRAPOLATION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>

namespace PhysBAM{

template<class TV>
class STRAWMAN_EXAMPLE : public EXAMPLE<TV>,LEVELSET_CALLBACKS<GRID<TV> >
{
    typedef EXAMPLE<TV> BASE;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,2> TV_INT;
    typedef GRID<TV> T_GRID;
    typedef ARRAY<T,TV_INT> T_ARRAY_SCALAR;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef BOUNDARY_PHI_WATER<TV> T_BOUNDARY_PHI_WATER;

public:
    using BASE::Time_At_Frame;using BASE::parse_args;
private:
    using BASE::frame_title;using BASE::stream_type;
    using BASE::last_frame;using BASE::frame_rate;
    using BASE::output_directory;

    int resolution,order,frame;
    T_GRID grid;
    T_ARRAY_SCALAR rho,rho_tmp;
    T_ARRAY_SCALAR phi,phi_tmp;
public:
    ARRAY<TV,TV_INT> velocity;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<bool,FACE_INDEX<TV::dimension> > face_valid_mask;

    enum UNCOVERED_CELL_METHOD {ANALYTIC,EXTRAPOLATION,GFM,NEW_GFM};
    UNCOVERED_CELL_METHOD method;
    bool fill_ghost_region;

    ADVECTION_CONSERVATIVE_ENO<T_GRID,T> advection_scheme;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T> passive_advection_scheme;
    BOUNDARY_LINEAR_EXTRAPOLATION<TV,T> boundary;

    // Initial Conditions
    T initial_distance;          // d0
    T solid_velocity;            // vS
    T fluid_tangential_velocity; // u

    RIGID_GEOMETRY_COLLECTION<TV> rigid_geometry_collection;
    T_GRID_BASED_COLLISION_GEOMETRY collision_bodies_affecting_fluid;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID> *pls_evolution;
    T_BOUNDARY_PHI_WATER pls_boundary;
    bool opt_analytic,opt_extrapolation,opt_gfm,opt_new_gfm;

    inline TV fluid_velocity_field(const TV& X,const T& time) const
    {return TV(fluid_tangential_velocity,-X(2)/(initial_distance+solid_velocity*time));}

    inline T f(const T& xbar) const
    {return (T)-1*xbar;}

    inline T analytic_solution(const TV& X,const T& time) const
    {return f(X.x - fluid_tangential_velocity*time) * X.y / ((solid_velocity*time+initial_distance)*(solid_velocity*time+initial_distance));}

    T initial_phi(const TV& X) const
    {TV interface_normal((T)-.342,(T).94),interface_point((T)1,(T)-.4);return max(solid_position(0)-X.y,TV::Dot_Product(interface_normal,interface_point-X));}

    inline T solid_position(const T& time) const
    {return initial_distance+solid_velocity*time;}

    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE {
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
            V_levelset(iterator.Full_Index())=fluid_velocity_field(iterator.Location(),time)(iterator.Axis());
    }
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<GRID<TV> >& levelset_multiple,ARRAY<T,FACE_INDEX<TV::dimension> >& V_levelset,const T time) const PHYSBAM_OVERRIDE {
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
            V_levelset(iterator.Full_Index())=fluid_velocity_field(iterator.Location(),time)(iterator.Axis());
    }

  public:
    STRAWMAN_EXAMPLE()
        :BASE(STREAM_TYPE(T())),resolution(100),order(1),frame(0),fill_ghost_region(true),initial_distance((T)-.5),solid_velocity((T)-1),
        fluid_tangential_velocity(0),rigid_geometry_collection(0),collision_bodies_affecting_fluid(grid),pls_evolution(0),
        opt_analytic(false),opt_extrapolation(false),opt_gfm(false),opt_new_gfm(false)
    {}

    inline T CFL(const T time) const
    {return .6/((T)sqrt(4*solid_velocity*solid_velocity + fluid_tangential_velocity*fluid_tangential_velocity)*(T)2*(T)resolution);}

//#####################################################################
// Function Register_Options
//#####################################################################
virtual void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add("-analytic",&opt_analytic,"analytic");
    parse_args->Add("-extrapolation",&opt_extrapolation,"extrapolation");
    parse_args->Add("-gfm",&opt_gfm,"gfm");
    parse_args->Add("-new_gfm",&opt_new_gfm,"new gfm");
    parse_args->Add("-fluid_velocity",&fluid_tangential_velocity,"value","fluid tangential velocity");
    parse_args->Add("-resolution",&resolution,"value","resolution");
    parse_args->Add("-order",&order,"value","order");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
virtual void Parse_Options() PHYSBAM_OVERRIDE
{
    last_frame=200;
    frame_rate=(T)400;
    BASE::Parse_Options();
    if(opt_analytic) method=ANALYTIC;
    else if(opt_extrapolation) method=EXTRAPOLATION;
    else if(opt_gfm) method=GFM;
    else method=NEW_GFM;

    advection_scheme.Set_Order(order);
    output_directory=STRING_UTILITIES::string_sprintf("Strawman_Example/Solution_%f_Resolution_%d_Order_%d",fluid_tangential_velocity,resolution,order);
    switch(method){
      case ANALYTIC:output_directory+="_Analytic";break;
      case EXTRAPOLATION:output_directory+="_Extrapolated";break;
      case GFM:output_directory+="_GFM";break;
      case NEW_GFM:output_directory+="_New";break;}

    grid.Initialize(TV_INT()+resolution,RANGE<TV>(TV(0,(T)-1),TV((T)1,0)),true);
    rho.Resize(grid.Domain_Indices(3)); rho_tmp.Resize(grid.Domain_Indices(3));
    phi.Resize(grid.Domain_Indices(3)); phi_tmp.Resize(grid.Domain_Indices(3));

    velocity.Resize(grid.Domain_Indices(3));
    face_velocities.Resize(grid.Domain_Indices(3));
    face_valid_mask.Resize(grid.Domain_Indices(3));

    pls_evolution=new PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>(grid,3,false);
}
//#####################################################################
// Function Initialize
//#####################################################################
void Initialize()
{
    { // All this just to get a line...
        RIGID_GEOMETRY<TV>& rigid_geometry=*new RIGID_GEOMETRY<TV>(rigid_geometry_collection,true);
        ANALYTIC_IMPLICIT_OBJECT<LINE_2D<T> >& implicit_structure=*new ANALYTIC_IMPLICIT_OBJECT<LINE_2D<T> >(LINE_2D<T>(TV(0,1),TV(0,0)));
        rigid_geometry.Add_Structure(implicit_structure);

        GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
        SEGMENTED_CURVE_2D<T>& segmented_curve=*SEGMENTED_CURVE_2D<T>::Create(particles);
        particles.Add_Elements(2);particles.X(1)=TV(-.1,0);particles.X(2)=TV((T)1.1,0);
        segmented_curve.mesh.number_nodes=2;segmented_curve.mesh.elements.Preallocate(1);
        segmented_curve.mesh.elements.Append(TV_INT(1,2));
        segmented_curve.Update_Segment_List();
        rigid_geometry.Add_Structure(segmented_curve);
        rigid_geometry_collection.Add_Rigid_Geometry(&rigid_geometry,stream_type,"",(T)1,false,false,false,false);

        rigid_geometry_collection.particles.structure_ids(1)=VECTOR<int,3>();
        rigid_geometry_collection.particles.structure_ids(1)(1)=rigid_geometry_collection.structure_list.Add_Element(rigid_geometry.structures(1));
        rigid_geometry_collection.particles.structure_ids(1)(2)=rigid_geometry_collection.structure_list.Add_Element(rigid_geometry.structures(2));
        rigid_geometry_collection.Update_Kinematic_Particles();
        rigid_geometry_collection.particles.frame(1).t.y=solid_position((T)0);

        collision_bodies_affecting_fluid.Initialize_Grids();
        collision_bodies_affecting_fluid.Add_Bodies(rigid_geometry_collection);

        collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        collision_bodies_affecting_fluid.Rasterize_Objects();
        collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*grid.dX.Min(),5);
        collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    }

    { // PLS
        T time=Time_At_Frame(0);
        collision_bodies_affecting_fluid.Compute_Psi_N_Zero_Velocity(face_valid_mask);
        for(FACE_ITERATOR<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
            face_valid_mask(iterator.Full_Index())=!face_valid_mask(iterator.Full_Index());
            face_velocities(iterator.Full_Index())=fluid_velocity_field(iterator.Location(),time)(iterator.Axis());}
        pls_boundary.Set_Velocity_Pointer(face_velocities);
        {
            pls_evolution->Initialize_Domain(grid);
            pls_evolution->Particle_Levelset(0).Set_Band_Width(6);
        }
        pls_evolution->time=time;
        pls_evolution->Set_CFL_Number((T).9);
        pls_evolution->Particle_Levelset(1).mpi_grid=0;
        pls_evolution->Levelset_Advection(1).Set_Custom_Advection(passive_advection_scheme);
        {
            pls_evolution->Initialize_Domain(grid);
            pls_evolution->Particle_Levelset(0).Set_Band_Width(6);
        }
        pls_evolution->Set_Number_Particles_Per_Cell(16);
        pls_evolution->Set_Levelset_Callbacks(*this);
        pls_evolution->Initialize_FMM_Initialization_Iterative_Solver(true);
        pls_evolution->Particle_Levelset(0).levelset.Set_Custom_Boundary(pls_boundary); // TODO: better boundary
        pls_evolution->Bias_Towards_Negative_Particles(false);
        pls_evolution->Particle_Levelset(1).Store_Unique_Particle_Id();
        pls_evolution->use_particle_levelset=true;
        pls_evolution->Particle_Levelset(0).levelset.Set_Face_Velocities_Valid_Mask(&face_valid_mask);
        //pls_evolution->Particle_Levelset(0).levelset.Set_Collision_Body_List(collision_bodies_affecting_fluid); // TODO: ?
        pls_evolution->Particle_Levelset(0).Set_Collision_Distance_Factors(.1,1);
        {
            pls_evolution->Initialize_Domain(grid);
            pls_evolution->Particle_Levelset(0).Set_Band_Width(6);
        }

        for(CELL_ITERATOR<TV> iterator(grid,3);iterator.Valid();iterator.Next())
            pls_evolution->Particle_Levelset(0).levelset.phi(iterator.Cell_Index())=initial_phi(iterator.Location());
        pls_evolution->Make_Signed_Distance();
        pls_evolution->Fill_Levelset_Ghost_Cells(time);
        pls_evolution->Set_Seed(2606);
        pls_evolution->Seed_Particles(time);
        pls_evolution->Delete_Particles_Outside_Grid();

        pls_boundary.Use_Extrapolation_Mode(false);
    }

    for(CELL_ITERATOR<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
        rho(iterator.Cell_Index())=analytic_solution(iterator.Location(),(T)0);
        phi(iterator.Cell_Index())=initial_phi(iterator.Location());
        velocity(iterator.Cell_Index())=fluid_velocity_field(iterator.Location(),(T)0);}
    Fill_Solid_Cells((T)0,solid_position(0),rho,phi,velocity);
}
//#####################################################################
// Function Fill_Solid_Cells
//#####################################################################
void Fill_Solid_Cells(const T time,const T solid_boundary,T_ARRAY_SCALAR& rho,T_ARRAY_SCALAR& phi,ARRAY<TV,TV_INT>& velocity){
    TV_INT solid_clamped_tn = grid.Index(TV((T).5,solid_boundary));
    int num_ghost_cells=3;
    if(fill_ghost_region) num_ghost_cells=solid_clamped_tn.y+3;

    // Fill the ghost region
    for(CELL_ITERATOR<TV> iterator(grid,3,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();TV location=iterator.Location();
        rho(cell_index)=analytic_solution(location,time);
        velocity(cell_index)=fluid_velocity_field(location,time);
        phi(cell_index)=initial_phi(location-time*velocity(cell_index));}

    if(solid_clamped_tn.y < 0) return;  // No solid cells outside the ghost region.

    // Extrapolate into the solid
    for(int offset=0;offset>=-(1-num_ghost_cells);--offset) for(int i=1;i <= resolution; ++i){
        rho(TV_INT(i,solid_clamped_tn.y+offset))=rho(TV_INT(i,solid_clamped_tn.y+1));
        velocity(TV_INT(i,solid_clamped_tn.y+offset))=velocity(TV_INT(i,solid_clamped_tn.y+1));}

    // Perform higher-order interpolation
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> rho_interpolation;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,TV> velocity_interpolation;
    for(int offset=0;offset>=(1-num_ghost_cells);--offset) for(int i=1;i <= resolution; ++i){
        TV location=grid.X(TV_INT(i,solid_clamped_tn.y+offset)),
           reflected_location=TV(location.x,(T)2*solid_boundary-location.y);
        rho(TV_INT(i,solid_clamped_tn.y+offset))=rho_interpolation.Clamped_To_Array(grid,rho,reflected_location);
        phi(TV_INT(i,solid_clamped_tn.y+offset))=rho_interpolation.Clamped_To_Array(grid,phi,reflected_location);}
        // We're given an analytic velocity field, so don't reflect it across the surface.
        //  velocity(TV_INT(i,solid_clamped_tn.y+offset)).y=(T)2*solid_velocity-velocity_interpolation.Clamped_To_Array(grid,velocity,reflected_location).y;}
}
//#####################################################################
// Function Fill_Uncovered_Cells
//#####################################################################
void Fill_Uncovered_Cells(const T& dt,const T& time,T_ARRAY_SCALAR& scalar_field,T_ARRAY_SCALAR& scalar_field_tmp){
    int num_ghost_cells=3;

    const T solid_boundary=solid_position(time);
    TV_INT solid_clamped_tn = grid.Index(TV((T).5,solid_boundary));
    if(solid_clamped_tn != grid.Index(TV((T).5,solid_position(time+dt))) && solid_clamped_tn.y > 0){
        switch(method){
        case ANALYTIC:
        case GFM:
            break;
        case EXTRAPOLATION:
            for(int i=1-num_ghost_cells;i<=resolution+num_ghost_cells;++i){TV_INT cell_index(i,solid_clamped_tn.y);
                scalar_field(cell_index)=scalar_field(TV_INT(i,solid_clamped_tn.y+1));}
            break;
        case NEW_GFM:{
            LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
            //for(int verticle_offset=0;verticle_offset>-3;--verticle_offset)
                for(int i=1-num_ghost_cells;i<=resolution+num_ghost_cells;++i){TV_INT cell_index(i,solid_clamped_tn.y);
                    TV stencil_center=grid.Center(cell_index)-(T)2*(solid_boundary-grid.Center(cell_index).y)*TV(fluid_tangential_velocity,solid_velocity);

                    T_GRID stencil_grid(TV_INT::All_Ones_Vector(),RANGE<TV>(-(T).5*grid.DX(),(T).5*grid.DX()),true);
                    T_ARRAY_SCALAR stencil_scalar_field(stencil_grid.Domain_Indices()),stencil_scalar_field_ghost(stencil_grid.Domain_Indices(3));
                    ARRAY<TV,TV_INT> stencil_velocity(stencil_grid.Domain_Indices(3));

                    TV surface_normal=TV(0,(T)1);
                    SYMMETRIC_MATRIX<T,2> transform=SYMMETRIC_MATRIX<T,2>::Identity_Matrix() -
                        (T)2*SYMMETRIC_MATRIX<T,2>(surface_normal.x*surface_normal.x,surface_normal.x*surface_normal.y,surface_normal.y*surface_normal.y);
                    for(CELL_ITERATOR<TV> iterator(stencil_grid,3);iterator.Valid();iterator.Next()){
                        TV location=stencil_center+transform*iterator.Location();
                        stencil_velocity(iterator.Cell_Index())=fluid_velocity_field(location,time);
                        stencil_scalar_field_ghost(iterator.Cell_Index())=interpolation.Clamped_To_Array(grid,scalar_field_tmp,location);}
                        stencil_scalar_field(TV_INT::All_Ones_Vector())=stencil_scalar_field_ghost(TV_INT::All_Ones_Vector());

                    advection_scheme.Update_Advection_Equation_Cell(stencil_grid,stencil_scalar_field,stencil_scalar_field_ghost,stencil_velocity,boundary,dt,time);
                    scalar_field(cell_index)=stencil_scalar_field(TV_INT::All_Ones_Vector());}}
            break;}}
}
//#####################################################################
// Function Euler_Step
//#####################################################################
void Euler_Step(const T dt, const T time){
    const T solid_boundary=solid_position(time);
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
        velocity(iterator.Cell_Index())=fluid_velocity_field(iterator.Location(),time);

    rigid_geometry_collection.particles.frame(1).t.y=solid_position(time);collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);
    rigid_geometry_collection.particles.frame(1).t.y=solid_position(time+dt);collision_bodies_affecting_fluid.Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,time+dt);
    collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE);
    collision_bodies_affecting_fluid.Restore_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);
    collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
    collision_bodies_affecting_fluid.Rasterize_Objects();

    //T maximum_fluid_speed=ARRAYS_COMPUTATIONS::Componentwise_Maxabs(velocity).Max();
    T maximum_fluid_speed=max(fluid_tangential_velocity,(T)2);
    T max_particle_collision_distance=pls_evolution->Particle_Levelset(1).max_collision_distance_factor*grid.dX.Max();
    collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*grid.dX.Max(),10);
    collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    face_valid_mask.Fill(false); collision_bodies_affecting_fluid.Compute_Psi_N_Zero_Velocity(face_valid_mask);
    for(FACE_ITERATOR<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
        face_valid_mask(iterator.Full_Index())=!face_valid_mask(iterator.Full_Index());
        face_velocities(iterator.Full_Index())=fluid_velocity_field(iterator.Location(),time)(iterator.Axis());}
    pls_boundary.Use_Extrapolation_Mode(false);
    {
        boundary.Fill_Ghost_Cells(grid,pls_evolution->Particle_Levelset(0).levelset.phi,phi_tmp,dt,time,3);
        LOG::Time("advancing levelset");
        pls_evolution->Advance_Levelset(dt);
        Fill_Uncovered_Cells(dt,time,pls_evolution->Particle_Levelset(0).levelset.phi,phi_tmp);
        // TODO(jontg): Fill uncovered cells new way here
        //  pls_evolution->Particle_Levelset(0).levelset.phi(iterator.Cell_Index())=initial_phi(iterator.Location());
        LOG::Time("advancing particles");
        pls_evolution->Advance_Particles(face_velocities,dt);
        pls_evolution->Modify_Levelset_And_Particles(&face_velocities);
    }

    if(method==ANALYTIC){
        for(CELL_ITERATOR<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
            rho(iterator.Cell_Index())=analytic_solution(iterator.Location(),time);
            phi(iterator.Cell_Index())=initial_phi(iterator.Location()-time*velocity(iterator.Cell_Index()));}
        rigid_geometry_collection.particles.frame(1).t.y=solid_position(time+dt);
        return;}

    boundary.Fill_Ghost_Cells(grid,rho,rho_tmp,dt,time,3);
    boundary.Fill_Ghost_Cells(grid,phi,phi_tmp,dt,time,3);
    Fill_Solid_Cells(time,solid_boundary,rho_tmp,phi_tmp,velocity);
    advection_scheme.Update_Advection_Equation_Cell(grid,rho,rho_tmp,velocity,boundary,dt,time);
    passive_advection_scheme.Update_Advection_Equation_Cell(grid,phi,phi_tmp,velocity,boundary,dt,time);

    Fill_Uncovered_Cells(dt,time,rho,rho_tmp);
    Fill_Uncovered_Cells(dt,time,phi,phi_tmp);
    for(CELL_ITERATOR<TV> iterator(grid,3,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();TV location=iterator.Location(),t0_location=location-(time+dt)*velocity(cell_index);
        rho(cell_index)=analytic_solution(location,time+dt);
        phi(cell_index)=initial_phi(t0_location);}
}
//#####################################################################
// Function Compute_Dt
//#####################################################################
T Compute_Dt(const T time,const T target_time,bool& done) const
{
    T dt=.6/((T)sqrt(4*solid_velocity*solid_velocity + fluid_tangential_velocity*fluid_tangential_velocity)*(T)2*(T)resolution);
    EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);
    return dt;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    std::string frame_folder=output_directory+STRING_UTILITIES::string_sprintf("/%d/",frame);
    if(frame==0) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/transverse_velocity",fluid_tangential_velocity);

    FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/center_velocities",velocity);
    FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/mac_velocities",face_velocities);
    rigid_geometry_collection.Write(stream_type,output_directory,frame);

    FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/density",rho);
    FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/phi",phi);

    { // PLS
        const PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& particle_levelset=pls_evolution->Particle_Levelset(0);
        FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"levelset",particle_levelset.levelset);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s",frame_folder.c_str(),"positive_particles"),particle_levelset.positive_particles);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s",frame_folder.c_str(),"negative_particles"),particle_levelset.negative_particles);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s",frame_folder.c_str(),"removed_positive_particles"),particle_levelset.removed_positive_particles);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s",frame_folder.c_str(),"removed_negative_particles"),particle_levelset.removed_negative_particles);
        FILE_UTILITIES::Write_To_Text_File(frame_folder+"last_unique_particle_id",particle_levelset.last_unique_particle_id);

        FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/psi_N",face_valid_mask);
    }

    FILE_UTILITIES::Write_To_Text_File(output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),frame_title);
}
//#####################################################################
};
}

#endif
