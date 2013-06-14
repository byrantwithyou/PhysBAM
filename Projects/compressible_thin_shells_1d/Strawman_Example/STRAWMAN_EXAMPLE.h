#ifndef __STRAWMAN_EXAMPLE__
#define __STRAWMAN_EXAMPLE__

#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_CONSERVATIVE_ENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/TYPED_STREAM.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_LINEAR_EXTRAPOLATION.h>

namespace PhysBAM{

template<class T,class RW>
class STRAWMAN_EXAMPLE
{
    typedef VECTOR<T,1> TV;
    typedef VECTOR<int,1> TV_INT;

    STREAM_TYPE stream_type;
    int resolution,frame;
    GRID<TV> grid;
    ARRAY<TV,TV_INT> velocity;
    ARRAY<T,TV_INT> rho_with_extrapolation,rho_fixed,rho_tmp;
    std::string output_directory;

    // Initial Conditions
    T initial_distance;          // d0
    T solid_velocity;            // vS

    RIGID_BODY_COLLECTION<TV> rigid_body_collection;

    inline TV fluid_velocity_field(const TV& X,const T& time) const
    {return TV(-X(2)/(initial_distance+solid_velocity*time));}

    inline T analytic_solution(const TV& X,const T& time) const
    {return -X.x / ((solid_velocity*time+initial_distance)*(solid_velocity*time+initial_distance));}

    inline T solid_position(const T& time) const
    {return initial_distance+solid_velocity*time;}

  public:
    STRAWMAN_EXAMPLE(const int resolution=100)
        : stream_type(T()),resolution(resolution),frame(0),grid(TV_INT(resolution),RANGE<TV>(TV((T)-1),TV(0)),true),
        velocity(grid.Domain_Indices(3)),rho_with_extrapolation(grid.Domain_Indices(3)),rho_fixed(grid.Domain_Indices(3)),
         rho_tmp(grid.Domain_Indices(3)),initial_distance((T)-.5),solid_velocity((T)-1),rigid_body_collection(0)
    {output_directory=STRING_UTILITIES::string_sprintf("Strawman_Example/Solution_Resolution_%d",resolution);}

//#####################################################################
// Function Initialize
//#####################################################################
void Initialize()
{
    // All this just to get a stupid line...
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    ANALYTIC_IMPLICIT_OBJECT<POINT_SIMPLEX_1D<T> >& implicit_structure=
        *new ANALYTIC_IMPLICIT_OBJECT<POINT_SIMPLEX_1D<T> >(POINT_SIMPLEX_1D<T>(TV(0),true));
    rigid_body.Add_Structure(implicit_structure);

    GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
    POINT_SIMPLICES_1D<T>& segmented_curve=*POINT_SIMPLICES_1D<T>::Create(particles);
    particles.Add_Elements(2);particles.X(1)=TV(0);particles.X(2)=TV(0);

    segmented_curve.mesh.number_nodes=2;segmented_curve.mesh.directions.Preallocate(2);
    segmented_curve.mesh.elements.Append(TV_INT(1));segmented_curve.mesh.directions(1)=false;
    segmented_curve.mesh.elements.Append(TV_INT(2));segmented_curve.mesh.directions(2)=true;
    segmented_curve.Update_Point_Simplex_List();
    rigid_body.Add_Structure(segmented_curve);
    rigid_body_collection.Add_Rigid_Body(&rigid_body,stream_type,"",(T)1,false,false,false,false);

    rigid_body_collection.rigid_body_particles.structure_ids(1)=VECTOR<int,3>();
    rigid_body_collection.rigid_body_particles.structure_ids(1)(1)=rigid_body_collection.structure_list.Add_Element(rigid_body.structures(1));
    rigid_body_collection.rigid_body_particles.structure_ids(1)(2)=rigid_body_collection.structure_list.Add_Element(rigid_body.structures(2));
    rigid_body_collection.Update_Kinematic_Particles();

    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/first_frame",0,"\n");
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);

    for(CELL_ITERATOR<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
        rho_fixed(iterator.Cell_Index())=analytic_solution(iterator.Location(),(T)0);
        rho_with_extrapolation(iterator.Cell_Index())=analytic_solution(iterator.Location(),(T)0);
        velocity(iterator.Cell_Index())=fluid_velocity_field(iterator.Location(),(T)0);}

    Fill_Solid_Cells((T)0,solid_position(0),rho_fixed,velocity);
    Fill_Solid_Cells((T)0,solid_position(0),rho_with_extrapolation,velocity);
}
//#####################################################################
// Function Fill_Solid_Cells
//#####################################################################
void Fill_Solid_Cells(const T time,const T solid_boundary,ARRAY<T,TV_INT>& rho,ARRAY<TV,TV_INT>& velocity){
    TV_INT solid_clamped_tn = grid.Index(TV(solid_boundary));

    // Fill the ghost region
    for(CELL_ITERATOR<TV> iterator(grid,3,GRID<TV>::GHOST_REGION);iterator.Valid();iterator.Next()){
        if(iterator.Cell_Index().x < 0){TV_INT reference_point=TV_INT(1);  // Outflow boundary conditions on the left side
            rho_fixed(iterator.Cell_Index())=rho_fixed(reference_point);
            rho_with_extrapolation(iterator.Cell_Index())=rho_fixed(reference_point);}
        else{
            rho_fixed(iterator.Cell_Index())=analytic_solution(iterator.Location(),time);
            rho_with_extrapolation(iterator.Cell_Index())=analytic_solution(iterator.Location(),time);}
        velocity(iterator.Cell_Index())=fluid_velocity_field(iterator.Location(),time);}

    if(solid_clamped_tn.x < 0) return;

    // Extrapolate into the solid
    for(int offset=0;offset>=-2;--offset){
        rho(TV_INT(solid_clamped_tn.x+offset))=rho(TV_INT(solid_clamped_tn.x+1));
        velocity(TV_INT(solid_clamped_tn.x+offset))=velocity(TV_INT(solid_clamped_tn.x+1));}

    // Perform higher-order interpolation
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> rho_interpolation;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> velocity_interpolation;
    for(int offset=0;offset>=-2;--offset){
        TV location=grid.X(TV_INT(solid_clamped_tn.x+offset)),
           reflected_location=TV((T)2*solid_boundary-location.x);
        rho(TV_INT(solid_clamped_tn.x+offset))=rho_interpolation.Clamped_To_Array(grid,rho,reflected_location);
        velocity(TV_INT(solid_clamped_tn.x+offset)).x=(T)2*solid_velocity-velocity_interpolation.Clamped_To_Array(grid,velocity,reflected_location).x;}
}
//#####################################################################
// Function Euler_Step
//#####################################################################
void Euler_Step(const T dt, const T time){
    ADVECTION_CONSERVATIVE_ENO<GRID<TV>,T> advection_scheme;advection_scheme.Set_Order(1);
    BOUNDARY_LINEAR_EXTRAPOLATION<TV,T> boundary;

    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
        velocity(iterator.Cell_Index())=fluid_velocity_field(iterator.Location(),time);
    const T solid_boundary=solid_position(time);

    boundary.Fill_Ghost_Cells(grid,rho_with_extrapolation,rho_tmp,dt,time,3);
    Fill_Solid_Cells(time,solid_boundary,rho_tmp,velocity);
    advection_scheme.Update_Advection_Equation_Cell(grid,rho_with_extrapolation,rho_tmp,velocity,boundary,dt,time);

    boundary.Fill_Ghost_Cells(grid,rho_fixed,rho_tmp,dt,time,3);
    Fill_Solid_Cells(time,solid_boundary,rho_tmp,velocity);
    advection_scheme.Update_Advection_Equation_Cell(grid,rho_fixed,rho_tmp,velocity,boundary,dt,time);

    TV_INT solid_clamped_tn = grid.Index(TV(solid_boundary));
    if(solid_clamped_tn != grid.Index(TV(solid_position(time+dt))) && solid_clamped_tn.x > 0){
        //Write_Output_Files(time+dt,"Before filling uncovered cells");
        LOG::cout<<"Crossover detected... "<<solid_clamped_tn<<", "<<grid.Index(TV(solid_position(time+dt)))<<std::endl;
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> rho_interpolation;
        rho_with_extrapolation(solid_clamped_tn)=rho_with_extrapolation(TV_INT(solid_clamped_tn.x+1));
        {
            TV stencil_center=grid.Center(solid_clamped_tn)-(T)2*(solid_boundary-grid.Center(solid_clamped_tn).x)*TV(-1);
            GRID<TV> stencil_grid(TV_INT::All_Ones_Vector(),RANGE<TV>(-(T).5*grid.DX(),(T).5*grid.DX()),true);
            ARRAY<T,TV_INT> stencil_data(stencil_grid.Domain_Indices()),stencil_data_ghost(stencil_grid.Domain_Indices(3));
            ARRAY<TV,TV_INT> stencil_velocity(stencil_grid.Domain_Indices(3));

            TV surface_normal=TV((T)1);
            MATRIX<T,1> transform=MATRIX<T,1>::Identity_Matrix() - (T)2*MATRIX<T,1>(surface_normal.x*surface_normal.x);
            for(CELL_ITERATOR<TV> iterator(stencil_grid,3);iterator.Valid();iterator.Next()){
                TV location=stencil_center+transform*iterator.Location();
                stencil_velocity(iterator.Cell_Index())=fluid_velocity_field(location,time);
                stencil_data_ghost(iterator.Cell_Index())=rho_interpolation.Clamped_To_Array(grid,rho_tmp,location);}
            stencil_data(TV_INT(1))=stencil_data_ghost(TV_INT(1));

            advection_scheme.Update_Advection_Equation_Cell(stencil_grid,stencil_data,stencil_data_ghost,stencil_velocity,boundary,dt,time);
            rho_fixed(solid_clamped_tn)=stencil_data(TV_INT(1));
        }}
}
//#####################################################################
// Function Run
//#####################################################################
void Run(){
    LOG::cout<<"Run() STARTING"<<std::endl;
    const T final_time=(T).5;
    T time=0,dt=.6/((T)sqrt(4*solid_velocity*solid_velocity)*(T)2*(T)resolution);

    Write_Output_Files(time,"Initial Setup");
    while (time < final_time){
        if(time+(T)2*dt>final_time && time+dt<final_time-(T)1e-12) dt = (final_time-time)/(T)2;
        Euler_Step(dt,time);
        Write_Output_Files(time+dt,STRING_UTILITIES::string_sprintf("Frame %d",frame));
        time+=dt;
    }
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const T& time,const std::string& frame_title)
{
    std::string frame_folder=output_directory+STRING_UTILITIES::string_sprintf("/%d/",frame);
    FILE_UTILITIES::Create_Directory(frame_folder);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/time",output_directory.c_str(),frame),time);

    FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/density",rho_with_extrapolation);
    FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/rho_fixed",rho_fixed);
    FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/center_velocities",velocity);

    rigid_body_collection.rigid_body_particles.frame(1).t.x=solid_position(time);
    rigid_body_collection.Write(stream_type,output_directory,frame);

    ARRAY<T,TV_INT> rho_analytic(grid.Domain_Indices(3));
    for(CELL_ITERATOR<TV> iterator(grid,3);iterator.Valid();iterator.Next()){
        rho_analytic(iterator.Cell_Index())=analytic_solution(iterator.Location(),time);}
    FILE_UTILITIES::Write_To_File(stream_type,frame_folder+"/analytic_rho",rho_analytic);

    TV_INT solid_clamped_tn = grid.Index(TV(solid_position(time)));
    rho_with_extrapolation(solid_clamped_tn)=rho_with_extrapolation(solid_clamped_tn+1);
    rho_fixed(solid_clamped_tn)=rho_fixed(solid_clamped_tn+1);
    rho_analytic(solid_clamped_tn)=rho_analytic(solid_clamped_tn+1);

    FILE_UTILITIES::Write_To_Text_File(output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),frame_title);

    frame++;
}
//#####################################################################
};
}

#endif
