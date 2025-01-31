//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET
//##################################################################### 
#include <Core/Math_Tools/pow.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Level_Sets/LEVELSET_COLLIDABLE.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//##################################################################### 
template<class TV> PARTICLE_LEVELSET<TV>::
PARTICLE_LEVELSET(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_ghost_cells_input)
    :only_use_negative_particles(false),use_removed_negative_particles(false),use_removed_positive_particles(false),last_unique_particle_id(0),
    particle_pool(template_particles),reincorporate_removed_particles_everywhere(false),save_removed_particle_times(false),
    delete_positive_particles_crossing_bodies(true),number_of_ghost_cells(number_of_ghost_cells_input),levelset(grid_input,phi_input,collision_body_list_input,number_of_ghost_cells_input),
    deletion_list(levelset.grid.Domain_Indices())
{   
    Set_Outside_Particle_Distance();
    Set_Maximum_Iterations_For_Attraction();
    Bias_Towards_Negative_Particles(false);
    Store_Unique_Particle_Id(false);
    Set_Collision_Distance_Factors();
    Set_Velocity_Interpolation_Collidable();

    Set_Minimum_Particle_Radius((T).1*levelset.grid.dX.Min());Set_Maximum_Particle_Radius((T).5*levelset.grid.dX.Min());
    Set_Band_Width();
    Use_Removed_Positive_Particles(false);
    Use_Removed_Negative_Particles(false);
    Set_Number_Particles_Per_Cell(pow<TV::m>(4));
}   
//#####################################################################
// Destructor
//##################################################################### 
template<class TV> PARTICLE_LEVELSET<TV>::
~PARTICLE_LEVELSET()
{
    removed_negative_particles.Delete_Pointers_And_Clean_Memory();
    removed_positive_particles.Delete_Pointers_And_Clean_Memory();
    positive_particles.Delete_Pointers_And_Clean_Memory();
    negative_particles.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize_Particle_Levelset_Grid_Values
//#####################################################################
template<class TV> void PARTICLE_LEVELSET<TV>::
Initialize_Particle_Levelset_Grid_Values()
{
    positive_particles.Resize(levelset.grid.Block_Indices(number_of_ghost_cells));
    negative_particles.Resize(levelset.grid.Block_Indices(number_of_ghost_cells));
    Use_Removed_Positive_Particles(use_removed_positive_particles);
    Use_Removed_Negative_Particles(use_removed_negative_particles);
    Set_Minimum_Particle_Radius((T).1*levelset.grid.dX.Min());
    Set_Maximum_Particle_Radius((T).5*levelset.grid.dX.Min());
    if(half_band_width && levelset.grid.dX.Min()) Set_Band_Width(half_band_width/((T).5*levelset.grid.dX.Min()));
    else Set_Band_Width();
    levelset.Initialize_Levelset_Grid_Values();
}
//FIX ALL OF THESE TO USE CORRECT PARTICLE TYPE
//#####################################################################
// Function Store_Unique_Particle_Id
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Store_Unique_Particle_Id(const bool store_unique_particle_id_input)
{
    store_unique_particle_id=store_unique_particle_id_input;
    if(store_unique_particle_id){
        template_particles.template Add_Array<int>("id");
        template_removed_particles.template Add_Array<int>("id");}
    else{
        template_particles.Remove_Array("id");
        template_removed_particles.Remove_Array("id");}
}
//#####################################################################
// Function Add_Particle
//##################################################################### 
template<class TV> PARTICLE_LEVELSET_PARTICLES<TV>* PARTICLE_LEVELSET<TV>::
Allocate_Particles(const PARTICLE_LEVELSET_PARTICLES<TV>& clone_particles)
{
    assert(!dynamic_cast<const PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&clone_particles));
    return particle_pool.Allocate_Particle();
}
//#####################################################################
// Function Add_Particle
//##################################################################### 
template<class TV> PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* PARTICLE_LEVELSET<TV>::
Allocate_Particles(const PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& clone_particles)
{
    return clone_particles.Clone();
}
//#####################################################################
// Function Add_Particle
//##################################################################### 
template<class TV> int PARTICLE_LEVELSET<TV>::
Add_Particle(PARTICLE_LEVELSET_PARTICLES<TV>*& particles)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(particles));
    int index;particles=&particle_pool.Add_Particle(*particles,index);
    return index;
}
//#####################################################################
// Function Add_Particle
//##################################################################### 
template<class TV> int PARTICLE_LEVELSET<TV>::
Add_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& particles)
{
    return particles->Add_Element();
}
//#####################################################################
// Function Add_Particle_To_Deletion_List
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Add_Particle_To_Deletion_List(ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >& deletion_list,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&particles));
    deletion_list.Append(PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int>(&particles,index));
    assert(index<particles.Size());
}
//#####################################################################
// Function Add_Particle_To_Deletion_List
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Add_Particle_To_Deletion_List(ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >& deletion_list,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int index)
{
    particles.Add_To_Deletion_List(index);
}
//#####################################################################
// Function Delete_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Delete_Particle_And_Clean_Memory(PARTICLE_LEVELSET_PARTICLES<TV>* head_particles,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&particles));
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(head_particles));
    if(Delete_Particle(particles,index)){
        if(head_particles->next!=0){
            PARTICLE_LEVELSET_PARTICLES<TV> *parent_particles=head_particles;
            while(parent_particles->next->next) parent_particles=parent_particles->next;
            Free_Particle_And_Clear_Pointer(parent_particles->next);}}
}
//#####################################################################
// Function Delete_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Delete_Particle_And_Clean_Memory(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* head_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int index)
{
    Delete_Particle(particles,index);
}
//#####################################################################
// Function Delete_Particle
//##################################################################### 
template<class TV> bool PARTICLE_LEVELSET<TV>::
Delete_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&particles));
    PARTICLE_LEVELSET_PARTICLES<TV> *particles_link=&particles;
    while(particles_link->next){particles_link=particles_link->next;}
    particles.Copy_Element(*particles_link,particles_link->Size()-1,index);
    particles_link->Delete_Element(particles_link->Size()-1);
    return particles_link->Size()==0;
}
//#####################################################################
// Function Delete_Particle
//##################################################################### 
template<class TV> bool PARTICLE_LEVELSET<TV>::
Delete_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int index)
{
    particles.Delete_Element(index);
    return(particles.Size()==0);
}
//#####################################################################
// Function Delete_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Delete_Particles_From_Deletion_List(ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >& deletion_list,PARTICLE_LEVELSET_PARTICLES<TV>& particles,bool verbose)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&particles));
    int level=0;verbose=true;
    PARTICLE_LEVELSET_PARTICLES<TV> *parent_particles=0,*particles_link=&particles,*deletion_particles=0;
    while(particles_link->next){parent_particles=particles_link;particles_link=particles_link->next;level++;}
    deletion_particles=particles_link;
    while(deletion_list.m>0){
        ARRAY<int> deletion_list_local;
        for(int i=deletion_list.m-1;i>=0;i--) if(deletion_particles==deletion_list(i).x){deletion_list_local.Append(deletion_list(i).y);deletion_list.Remove_Index_Lazy(i);}
        deletion_list_local.Sort();
        for(int i=deletion_list_local.m-1;i>=0;i--){int index=deletion_list_local(i);
            deletion_particles->Copy_Element(*particles_link,particles_link->Size()-1,index);
            particles_link->Delete_Element(particles_link->Size()-1);
            if(particles_link->Size()==0 && parent_particles){
                particle_pool.Free_Particle(particles_link);parent_particles->next=0;particles_link=&particles;parent_particles=0;
                while(particles_link->next){parent_particles=particles_link;particles_link=particles_link->next;}}}
        deletion_particles=&particles;level--;
        if(level<0){assert(deletion_list.m==0);return;}
        for(int i=0;i<level;i++){deletion_particles=deletion_particles->next;}}
}
//#####################################################################
// Function Delete_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Delete_Particles_From_Deletion_List(ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >& deletion_list,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles)
{
    particles.Delete_Elements_On_Deletion_List();
}
//#####################################################################
// Function Get_Particle_Link
//##################################################################### 
template<class TV> PARTICLE_LEVELSET_PARTICLES<TV>& PARTICLE_LEVELSET<TV>::
Get_Particle_Link(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int absolute_index,int& index_in_link)
{
    PARTICLE_LEVELSET_PARTICLES<TV>* particles_link=&particles;index_in_link=absolute_index;
    while(index_in_link>=particle_pool.number_particles_per_cell){
        assert(particles_link->next);assert(particles_link->Size()==particle_pool.number_particles_per_cell);
        particles_link=particles_link->next;index_in_link-=particle_pool.number_particles_per_cell;}
    assert(index_in_link<particles_link->Size());
    return *particles_link;
}
//#####################################################################
// Function Copy_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Copy_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_PARTICLES<TV>& to_particles,const int from_absolute_index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&from_particles));
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&to_particles));
    int to_index;PARTICLE_LEVELSET_PARTICLES<TV>& to_particles_link=particle_pool.Add_Particle(to_particles,to_index);
    int from_index;PARTICLE_LEVELSET_PARTICLES<TV>& from_particles_link=Get_Particle_Link(from_particles,from_absolute_index,from_index);
    to_particles_link.Copy_Element(from_particles_link,from_index,to_index);
}
//#####################################################################
// Function Copy_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Copy_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& to_particles,const int from_absolute_index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&from_particles));
    int from_index;PARTICLE_LEVELSET_PARTICLES<TV>& from_particles_link=Get_Particle_Link(from_particles,from_absolute_index,from_index);
    to_particles.Append(from_particles_link,from_index);
}
//#####################################################################
// Function Copy_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Copy_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_PARTICLES<TV>& to_particles,const int from_index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&to_particles));
    int to_index;PARTICLE_LEVELSET_PARTICLES<TV>& to_particles_link=particle_pool.Add_Particle(to_particles,to_index);
    to_particles_link.Copy_Element(from_particles,from_index,to_index);
}
//#####################################################################
// Function Copy_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Copy_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& to_particles,const int from_index)
{
    to_particles.Append(from_particles,from_index);
}
//#####################################################################
// Function Move_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Move_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_PARTICLES<TV>& to_particles,const int from_absolute_index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&from_particles));
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&to_particles));
    int to_index;PARTICLE_LEVELSET_PARTICLES<TV>& to_particles_link=particle_pool.Add_Particle(to_particles,to_index);
    int from_index;PARTICLE_LEVELSET_PARTICLES<TV>& from_particles_link=Get_Particle_Link(from_particles,from_absolute_index,from_index);
    to_particles_link.Copy_Element(from_particles_link,from_index,to_index);Delete_Particle(from_particles_link,from_index);
}
//#####################################################################
// Function Move_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Move_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& to_particles,const int from_absolute_index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&from_particles));
    int to_index;to_index=to_particles.Add_Element();
    int from_index;PARTICLE_LEVELSET_PARTICLES<TV>& from_particles_link=Get_Particle_Link(from_particles,from_absolute_index,from_index);
    to_particles.Copy_Element(from_particles_link,from_index,to_index);Delete_Particle(from_particles_link,from_index);
}
//#####################################################################
// Function Move_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Move_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_PARTICLES<TV>& to_particles,const int from_index)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(&to_particles));
    int to_index;PARTICLE_LEVELSET_PARTICLES<TV>& to_particles_link=particle_pool.Add_Particle(to_particles,to_index);
    to_particles_link.Copy_Element(from_particles,from_index,to_index);Delete_Particle(from_particles,from_index);
}
//#####################################################################
// Function Move_Particle
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Move_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& to_particles,const int from_index)
{
    to_particles.Take(from_particles,from_index);
}
//#####################################################################
// Function Compact_Particles
//##################################################################### 
template<class TV> template<class T_PARTICLES> void PARTICLE_LEVELSET<TV>::
Compact_Particles(T_PARTICLES& particles)
{
    particles.Compact();
}
//#####################################################################
// Function Allocate_Particle
//##################################################################### 
template<class TV> static inline PARTICLE_LEVELSET_PARTICLES<TV>*
Allocate_Particle_Helper(PARTICLE_LEVELSET<TV>& particle_levelset,const PARTICLE_LEVELSET_PARTICLES<TV>*)
{
    return particle_levelset.particle_pool.Allocate_Particle();
}
template<class TV> static inline PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*
Allocate_Particle_Helper(PARTICLE_LEVELSET<TV>& particle_levelset,const PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*)
{
    return particle_levelset.template_removed_particles.Clone();
}
template<class TV> template<class T_PARTICLES> T_PARTICLES* PARTICLE_LEVELSET<TV>::
Allocate_Particle()
{
    return Allocate_Particle_Helper(*this,(T_PARTICLES*)0);
}
//#####################################################################
// Function Free_Particle_And_Clear_Pointer
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Free_Particle_And_Clear_Pointer(PARTICLE_LEVELSET_PARTICLES<TV>*& particles)
{
    assert(!dynamic_cast<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>(particles));
    particle_pool.Free_Particle(particles);particles=0;
}
//#####################################################################
// Function Free_Particle_And_Clear_Pointer
//##################################################################### 
template<class TV> void PARTICLE_LEVELSET<TV>::
Free_Particle_And_Clear_Pointer(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& particles)
{
    delete particles;particles=0;
}
//#####################################################################
// Function Adjust_Particle_For_Object
//##################################################################### 
template<class TV> bool PARTICLE_LEVELSET<TV>::
Adjust_Particle_For_Objects(TV& X,TV& V,const T r, const T collision_distance,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) const
{
    TV X_old=X,X_new=X_old+dt*V;
    if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE){
        if(delete_positive_particles_crossing_bodies) return !levelset.collision_body_list->Any_Crossover(X_old,X_new,dt);
        return true;}

    COLLISION_GEOMETRY_ID body_id;int aggegate_id;T phi_old,phi;TV N,V_object;
    if(levelset.collision_body_list->Get_Body_Penetration(X_old,X_new,levelset.small_number,dt,body_id,aggegate_id,phi_old,phi,N,V_object)){
        X_new+=(max((T)0,-phi+collision_distance))*N;
        T relative_normal_velocity=TV::Dot_Product(V-V_object,N);if(relative_normal_velocity<0) V-=relative_normal_velocity*N;}
    if(levelset.collision_body_list->Get_Body_Penetration(X_new,X_new,collision_distance,1,body_id,aggegate_id,phi_old,phi,N,V_object)){
        X_new+=(max((T)0,-phi+collision_distance))*N;
        T relative_normal_velocity=TV::Dot_Product(V-V_object,N);if(relative_normal_velocity<0) V-=relative_normal_velocity*N;}
    X=X_new-dt*V;
    if(levelset.collision_body_list->Any_Crossover(X_old,X_new,dt)) return false;
    return true;
}
//#####################################################################
namespace PhysBAM{
#define INSTANTIATION_HELPER(TV) \
    template class PARTICLE_LEVELSET<TV>; \
    template PARTICLE_LEVELSET_PARTICLES<TV>* PARTICLE_LEVELSET<TV>::Allocate_Particle<PARTICLE_LEVELSET_PARTICLES<TV> >(); \
    template PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* PARTICLE_LEVELSET<TV>::Allocate_Particle<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(); \
    template void PARTICLE_LEVELSET<TV>::Compact_Particles<PARTICLE_LEVELSET_PARTICLES<TV> >(PARTICLE_LEVELSET_PARTICLES<TV>&); \
    template void PARTICLE_LEVELSET<TV>::Compact_Particles<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> >(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>&);
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(P(VECTOR<float,1>));
INSTANTIATION_HELPER(P(VECTOR<float,2>));
INSTANTIATION_HELPER(P(VECTOR<float,3>));
INSTANTIATION_HELPER(P(VECTOR<double,1>));
INSTANTIATION_HELPER(P(VECTOR<double,2>));
INSTANTIATION_HELPER(P(VECTOR<double,3>));
}
