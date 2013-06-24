//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <Tools/Interpolation/INTERPOLATED_COLOR_MAP.h>
#include <Deformables/Collisions_And_Interactions/VOLUME_COLLISIONS.h>
#include <Deformables/Forces/COLLISION_AREA_PENALTY_FORCE.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
namespace PhysBAM{template<class TV,class ATTR> void Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr);}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COLLISION_AREA_PENALTY_FORCE<TV>::
COLLISION_AREA_PENALTY_FORCE(DEFORMABLE_PARTICLES<TV>& particles)
    :BASE(particles),volume_collisions(*new VOLUME_COLLISIONS<TV>),force_coefficient(1e6)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COLLISION_AREA_PENALTY_FORCE<TV>::
~COLLISION_AREA_PENALTY_FORCE()
{
    delete &volume_collisions;
}
//#####################################################################
// Function Add_Mesh
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Mesh(T_OBJECT& ta)
{
    volume_collisions.objects.Append(&ta);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR COLLISION_AREA_PENALTY_FORCE<TV>::
Potential_Energy(const T time) const
{
    return force_coefficient*sqr(volume_collisions.area);
}
//#####################################################################
// Function Use_Rest_State_For_Strain_Rate
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input)
{
}
//#####################################################################
// Function Limit_Time_Step_By_Strain_Rate
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input,const T max_strain_per_time_step_input)
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    volume_collisions.Compute_Collision_Triangles();
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(typename HASHTABLE<int,TV>::ITERATOR it(volume_collisions.gradient);it.Valid();it.Next())
        F(it.Key())-=2*force_coefficient*volume_collisions.area*it.Data();

    T mx=0;
    for(typename HASHTABLE<int,TV>::ITERATOR it(volume_collisions.gradient);it.Valid();it.Next())
        mx=max(mx,it.Data().Magnitude());

    INTERPOLATED_COLOR_MAP<T> color_map;
    color_map.Initialize_Colors(0,mx,false,true,false);

    for(typename HASHTABLE<int,TV>::ITERATOR it(volume_collisions.gradient);it.Valid();it.Next()){
        Add_Debug_Particle(this->particles.X(it.Key()),color_map(it.Data().Magnitude()));
        Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,-2*force_coefficient*volume_collisions.area*it.Data());}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class TV> int COLLISION_AREA_PENALTY_FORCE<TV>::
Velocity_Dependent_Forces_Size() const
{
    PHYSBAM_FATAL_ERROR();
    return 0;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(typename HASHTABLE<VECTOR<int,2>,MATRIX<T,TV::m> >::ITERATOR it(volume_collisions.hessian);it.Valid();it.Next())
        F(it.Key().x)-=2*force_coefficient*volume_collisions.area*(it.Data()*V(it.Key().y));

    T tot=0;
    for(typename HASHTABLE<int,TV>::ITERATOR it(volume_collisions.gradient);it.Valid();it.Next())
        tot+=TV::Dot_Product(it.Data(),V(it.Key()));
    for(typename HASHTABLE<int,TV>::ITERATOR it(volume_collisions.gradient);it.Valid();it.Next())
        F(it.Key())-=2*force_coefficient*tot*it.Data();
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR COLLISION_AREA_PENALTY_FORCE<TV>::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void COLLISION_AREA_PENALTY_FORCE<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{
}
namespace PhysBAM{
template class COLLISION_AREA_PENALTY_FORCE<VECTOR<float,2> >;
template class COLLISION_AREA_PENALTY_FORCE<VECTOR<float,3> >;
template class COLLISION_AREA_PENALTY_FORCE<VECTOR<double,2> >;
template class COLLISION_AREA_PENALTY_FORCE<VECTOR<double,3> >;
}
