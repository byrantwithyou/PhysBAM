//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_VELOCITY_TESTS 
//##################################################################### 
#ifndef __BOUNDARY_VELOCITY_TESTS__
#define __BOUNDARY_VELOCITY_TESTS__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template<class T,class RW=T>
class BOUNDARY_VELOCITY_TESTS:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Adjust_Phi_With_Sources;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Get_Source_Reseed_Mask;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Get_Source_Velocities;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;

    int test_number;
    
    bool use_source[2];
    BOX_2D<T> source[2];
    MATRIX<T,3> world_to_source[2];
    VECTOR_2D<T> source_velocity[2];

    bool use_object;
    ARRAY<T,VECTOR<int,2> > phi_object;
    ARRAY<VECTOR_2D<T> ,VECTOR<int,2> > V_object;
    LEVELSET_2D<T> object_levelset;

    PARAMETER_LIST parameter_list;
    bool custom_face_to_node;

    BOUNDARY_VELOCITY_TESTS(int test_number_input,int resolution)
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.WATER),object_levelset(fluids_parameters.grid,phi_object)
    {
        std::cout<<"Running Boundary Velocity Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;

        // set up the standard fluid environment
        first_frame=0;last_frame=1000;frame_rate=24;
        restart=false;restart_frame=0;
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[3][1]=true;fluids_parameters.domain_walls[3][2]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_debug_data=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.incompressible_iterations=1000;

        output_directory=STRING_UTILITIES::string_sprintf("Boundary_Velocity_Tests/Test_%d__Resolution_%d",test_number,resolution);

        for(int i=0;i<2;i++)use_source[i]=false;
        use_object=false;
        custom_face_to_node=false;
   
        fluids_parameters.grid.Initialize(48*(int)pow(2,resolution)+1,48*(int)pow(2,resolution)+1,0,1,0,1);

        if(test_number==1){
            use_object=true;
            BOX_2D<T> glass(.1,.9,0.1,2);MATRIX<T,2> xform=MATRIX<T,2>::Identity_Matrix();
            phi_object.Resize(fluids_parameters.grid,3);V_object.Resize(fluids_parameters.grid,3);
            for(int i=-2;i<=fluids_parameters.grid.m+3;i++) for(int j=-2;j<=fluids_parameters.grid.n+3;j++)
                phi_object(i,j)=-glass.Signed_Distance(xform*fluids_parameters.grid.X(i,j));
            object_levelset.Compute_Cell_Minimum_And_Maximum();}
        else if(test_number==2){
            use_object=true;
            BOX_2D<T> glass(.1,.9,0.1,2);MATRIX<T,2> xform=MATRIX<T,2>::Identity_Matrix();
            phi_object.Resize(fluids_parameters.grid,3);V_object.Resize(fluids_parameters.grid,3);
            for(int i=-2;i<=fluids_parameters.grid.m+3;i++) for(int j=-2;j<=fluids_parameters.grid.n+3;j++)
                phi_object(i,j)=-glass.Signed_Distance(xform*fluids_parameters.grid.X(i,j));
            object_levelset.Compute_Cell_Minimum_And_Maximum();}
        else if(test_number==3){
            use_object=true;
            BOX_2D<T> glass(-.4,.4,-.4,.4);MATRIX<T,3> xform=MATRIX<T,3>::Translation_Matrix(VECTOR_2D<T>(.5,.5))*MATRIX<T,3>::Rotation_Matrix_Z_Axis(pi/4);xform=xform.Inverse();
            phi_object.Resize(fluids_parameters.grid,3);V_object.Resize(fluids_parameters.grid,3);
            for(int i=-2;i<=fluids_parameters.grid.m+3;i++) for(int j=-2;j<=fluids_parameters.grid.n+3;j++)
                phi_object(i,j)=-glass.Signed_Distance(xform*fluids_parameters.grid.X(i,j));
            object_levelset.Compute_Cell_Minimum_And_Maximum();}
        else if(test_number==4){
            custom_face_to_node=true;
            fluids_parameters.grid.Initialize(221,221,-.05,1.05,-.05,1.05);
            use_object=true;T tol=1e-5;
            BOX_2D<T> glass(0+tol,1-tol,0+tol,1-tol);MATRIX<T,2> xform=MATRIX<T,2>::Identity_Matrix();
            phi_object.Resize(fluids_parameters.grid,3);V_object.Resize(fluids_parameters.grid,3);
            for(int i=-2;i<=fluids_parameters.grid.m+3;i++) for(int j=-2;j<=fluids_parameters.grid.n+3;j++)
                phi_object(i,j)=-glass.Signed_Distance(xform*fluids_parameters.grid.X(i,j));
            object_levelset.Compute_Cell_Minimum_And_Maximum();}

        std::string filename=STRING_UTILITIES::string_sprintf("Boundary_Velocity_Tests/example_%d.param",test_number);
        if(FILE_UTILITIES::File_Exists(filename)){std::cout << "Reading parameter file '" << filename << "'" << std::endl;parameter_list.Read(filename);}
        THIN_SHELLS_FLUID_COUPLING_UTILITIES<T,RW>::Set_Parameters_From_Parameter_List(*this,parameter_list);

        Set_Parameter(custom_face_to_node,"custom_face_to_node");
    }
    
    ~BOUNDARY_VELOCITY_TESTS() 
    {}

//#####################################################################
// Function Set_Parameter
//#####################################################################
template<class T2> void Set_Parameter(T2& parameter,const std::string& name)
{
    if(parameter_list.Is_Defined(name)){
        parameter=parameter_list.template Get_Parameter<T2>(name);
        std::cout << "[param] set " << name << " = " << parameter << std::endl;}
}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution.phi;

    ARRAY<T,VECTOR<int,2> >::copy(fluids_parameters.grid.min_dx_dy,fluids_parameters.particle_levelset_evolution.phi);

    if(test_number==1){
        CIRCLE<T> circle(VECTOR_2D<T>(.1,.7),.1);
        for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)phi(i,j)=circle.Signed_Distance(grid.X(i,j));}
    else if(test_number==2){
        for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)phi(i,j)=grid.y(j)-(T).400235234;}
    else if(test_number==3){
        for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)phi(i,j)=grid.y(j)-(T).400235234;}
    else if(test_number==4){
        CIRCLE<T> circle(VECTOR_2D<T>((T).5,(T).6),.15);
        for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)
            phi(i,j)=min(circle.Signed_Distance(fluids_parameters.grid.X(i,j)),grid.y(j)-(T).31);}
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    for(int i=0;i<2;i++) if(use_source[i]) Adjust_Phi_With_Source(source[i],world_to_source[i],time);
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects()
{
    if(use_object)
        fluids_parameters.Extrapolate_Phi_Into_Object(fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi,phi_object);
}
//#####################################################################
// Function Adjust_Phi_with_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    if(use_object) fluids_parameters.Adjust_Phi_With_Object(phi_object,V_object,time);
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(VECTOR_2D<T>& X,VECTOR_2D<T>& V,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{          
    if(use_object){
        LEVELSET_2D<T> levelset_object(fluids_parameters.grid,phi_object);
        fluids_parameters.Adjust_Particle_For_Object(levelset_object,V_object,X,V,particle_type,dt,time);
        return true;}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR_2D<T> >& particles,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T time)
{
    if(use_object){
        if(particle_type == PARTICLE_LEVELSET<T,VECTOR_2D<T> >::NEGATIVE || particle_type == PARTICLE_LEVELSET<T,VECTOR_2D<T> >::REMOVED_NEGATIVE){
            for(int k=particles.array_collection->Size();k>=1;k--) if(object_levelset.Lazy_Inside_Extended_Levelset(particles.X(k),-fluids_parameters.grid.dx)) particles.Delete_Particle(k);}
        else for(int k=particles.array_collection->Size();k>=1;k--) if(object_levelset.Lazy_Inside_Extended_Levelset(particles.X(k))) particles.Delete_Particle(k);}
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(use_object) fluids_parameters.Extrapolate_Velocity_Into_Object(phi_object,V_object,3,true,time);
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time)
{
    bool first=true;
    for(int i=0;i<2;i++) if(use_source[i]){Get_Source_Reseed_Mask(source[i],world_to_source[i],cell_centered_mask,first,time);first=false;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    for(int i=0;i<2;i++) if(use_source[i]) Get_Source_Velocities(source[i],world_to_source[i],source_velocity[i],time);
}
//#####################################################################
// Function Average_Face_Velocities_To_Nodes
//#####################################################################
#if 0
void Average_Face_Velocities_To_Nodes()
{
    if(!custom_face_to_node) SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Average_Face_Velocities_To_Nodes();
    std::cout << "custom face->node" << std::endl;
    ARRAY<bool,VECTOR<int,2> >& psi_N_u=fluids_parameters.incompressible.projection.elliptic_solver->psi_N_u;
    ARRAY<bool,VECTOR<int,2> >& psi_N_v=fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v;
    ARRAY<bool,VECTOR<int,2> >& psi_D=fluids_parameters.incompressible.projection.elliptic_solver->psi_D;
    GRID<TV> &grid=fluids_parameters.grid;
    GRID<TV> u_grid=fluids_parameters.incompressible.projection.u_grid,v_grid=fluids_parameters.incompressible.projection.v_grid;
    ARRAY<T,VECTOR<int,2> > &u=fluids_parameters.incompressible.projection.u,&v=fluids_parameters.incompressible.projection.v;
    ARRAY<VECTOR_2D<T> ,VECTOR<int,2> > &V=fluids_parameters.incompressible.V,&V_ghost=fluids_parameters.incompressible.V_ghost;
    ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >::put(V,V_ghost);
    ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >::copy(VECTOR_2D<T>(0,0),V);ARRAY<VECTOR_2D<char> ,VECTOR<int,2> > count(grid);
    
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) if(!psi_N_u(i,j) && !(psi_D(i,j)&&psi_D(i-1,j))){V(i,j).x+=u(i,j);V(i,j+1).x+=u(i,j);count(i,j).x++;count(i,j+1).x++;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) if(!psi_N_v(i,j) && !(psi_D(i,j)&&psi_D(i,j-1))){V(i,j).y+=v(i,j);V(i+1,j).y+=v(i,j);count(i,j).y++;count(i+1,j).y++;}
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){
        if(count(i,j).x>0) V(i,j).x/=(T)count(i,j).x;else V(i,j).x=0;//V_ghost(i,j).x;
        if(count(i,j).y>0) V(i,j).y/=(T)count(i,j).y;else V(i,j).y=0;/*V_ghost(i,j).y;*/}

    fluids_parameters.incompressible.psi_D_old=fluids_parameters.incompressible.projection.elliptic_solver->psi_D;
}
#endif
//#####################################################################
};      
}
#endif    


