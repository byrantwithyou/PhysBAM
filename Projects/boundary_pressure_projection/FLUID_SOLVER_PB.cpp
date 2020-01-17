//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Projection/BOUNDARY_CONDITION_DOUBLE_FINE.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Fluids/Fluids/FLUID_COLLECTION.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "FLUID_BC_PB.h"
#include "FLUID_BOUNDARY_VECTOR_PB.h"
#include "FLUID_REGIONS_PB.h"
#include "FLUID_SOLVER_PB.h"
#include "FLUID_STATE_PB.h"
#include "SOLID_FLUID_INTERFACE.h"
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_SOLVER_PB<TV>::
FLUID_SOLVER_PB()
{
}

//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_SOLVER_PB<TV>::
~FLUID_SOLVER_PB()
{
}

//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Initialize()
{
    driver->Initialize();
    auto old_cb=driver->example.get_unified_boundary_conditions;
}

//#####################################################################
// Function Write
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Write(int frame) const
{
    driver->Write_Output_Files();
}

//#####################################################################
// Function Read
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Read(int frame)
{
}

//#####################################################################
// Function Compute_Dt
//#####################################################################
template<class TV> auto FLUID_SOLVER_PB<TV>::
Compute_Dt(T time) const -> T
{
    return driver->Compute_Fluids_Dt(time);
}

//#####################################################################
// Function Simulate_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Simulate_Time_Step(FLUID_BOUNDARY_VECTOR<TV>* velocity,T time,T dt)
{
    auto* bc_v=dynamic_cast<FLUID_BOUNDARY_VECTOR_PB<TV>*>(velocity);
    auto* bc_fine=driver->example.fluids_parameters.bc_fine;

    ARRAY<char,TV_INT> bc_type,bc_type_current;
    HASHTABLE<TV_INT,PAIR<T,int> > bc_p;
    HASHTABLE<FACE_INDEX<TV::m>,PAIR<T,int> > bc_u;

    T bc_v_value=0;
    for(auto& h:bc_fine->bc_u)
        if(bc_v->V.Get(h.key,bc_v_value))
            h.data={bc_v_value,1};
            
    bc_type.Exchange(bc_fine->bc_type);
    bc_type_current.Exchange(bc_fine->bc_type_current);
    bc_p.Exchange(bc_fine->bc_p);
    bc_u.Exchange(bc_fine->bc_u);
    auto old_cb=driver->example.get_unified_boundary_conditions;
    driver->example.get_unified_boundary_conditions=
        [&,this](BOUNDARY_CONDITION_DOUBLE_FINE<TV>* bc_fine,const T time)
        {
            bc_type.Exchange(bc_fine->bc_type);
            bc_type_current.Exchange(bc_fine->bc_type_current);
            bc_p.Exchange(bc_fine->bc_p);
            bc_u.Exchange(bc_fine->bc_u);
        };
    driver->Advance_One_Time_Step(dt);
    driver->example.get_unified_boundary_conditions=old_cb;
}

//#####################################################################
// Function Predict_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Predict_Time_Step(T time,T dt)
{
    // Predict that the state does not change.
}

//#####################################################################
// Function Before_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Before_Time_Step(T time)
{
    driver->Setup_Fluids(time);
}
//#####################################################################
// Function After_Time_Step
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
After_Time_Step(T time,T dt)
{
}
//#####################################################################
// Function Before_Frame
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Before_Frame(int frame)
{
    driver->Preprocess_Frame(frame+1);
}
//#####################################################################
// Function After_Frame
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
After_Frame(int frame)
{
    driver->Postprocess_Frame(frame+1);
}
//#####################################################################
// Function Make_State
//#####################################################################
template<class TV> FLUID_STATE<TV>* FLUID_SOLVER_PB<TV>::
Make_State() const
{
    return new FLUID_STATE_PB<TV>;
}

//#####################################################################
// Function Make_BC
//#####################################################################
template<class TV> FLUID_BC<TV>* FLUID_SOLVER_PB<TV>::
Make_BC() const
{
    return new FLUID_BC_PB<TV>;
}

//#####################################################################
// Function Save
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Save(FLUID_STATE<TV>* fluid_state) const
{
    FLUID_STATE_PB<TV>& st=dynamic_cast<FLUID_STATE_PB<TV>&>(*fluid_state);
    st.face_velocities=driver->example.fluid_collection.incompressible_fluid_collection.face_velocities;
    st.pressure=driver->example.fluids_parameters.incompressible->projection.p;
}

//#####################################################################
// Function Restore
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Restore(const FLUID_STATE<TV>* fluid_state)
{
    const FLUID_STATE_PB<TV>& st=dynamic_cast<const FLUID_STATE_PB<TV>&>(*fluid_state);
    driver->example.fluid_collection.incompressible_fluid_collection.face_velocities=st.face_velocities;
    driver->example.fluids_parameters.incompressible->projection.p=st.pressure;
}

//#####################################################################
// Function Diff_u
//#####################################################################
template<class TV> auto FLUID_SOLVER_PB<TV>::
Diff_u(const FLUID_STATE<TV>* fluid_state) const -> T
{
    const FLUID_STATE_PB<TV>& st=dynamic_cast<const FLUID_STATE_PB<TV>&>(*fluid_state);
    return (driver->example.fluid_collection.incompressible_fluid_collection.face_velocities.array-st.face_velocities.array).Max_Abs();
}

//#####################################################################
// Function Diff_p
//#####################################################################
template<class TV> auto FLUID_SOLVER_PB<TV>::
Diff_p(const FLUID_STATE<TV>* fluid_state) const -> T
{
    const FLUID_STATE_PB<TV>& st=dynamic_cast<const FLUID_STATE_PB<TV>&>(*fluid_state);
    return (driver->example.fluids_parameters.incompressible->projection.p.array-st.pressure.array).Max_Abs();
}

//#####################################################################
// Function Get_Constraints
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Get_Constraints(const SOLID_FLUID_INTERFACE<TV>* interface,
    ARRAY<FLUID_BOUNDARY_VECTOR<TV>*>& array,ARRAY<T>& rhs,
    FLUID_REGIONS<TV>* regions) const
{
    typedef FLUID_BOUNDARY_VECTOR_PB<TV> V;
    const auto* es=driver->example.fluids_parameters.incompressible->projection.elliptic_solver;
    const auto& c=es->filled_region_colors;
    const auto& fv=driver->example.fluid_collection.incompressible_fluid_collection.face_velocities;
    const GRID<TV>& grid=*driver->example.fluids_parameters.grid;

    int num_regions=0;
    ARRAY<int> color_map(es->number_of_regions,use_init,-1);
    for(int i=0;i<color_map.m;i++)
        if(!es->filled_region_touches_dirichlet(i))
            color_map(i)=num_regions++;

    array.Resize(num_regions);
    for(auto&a:array) a=Make_Boundary_Vector();
    rhs.Resize(num_regions);

    FLUID_BOUNDARY_VECTOR_PB<TV>* dofs=new FLUID_BOUNDARY_VECTOR_PB<TV>;
    interface->Get_Boundary(dofs);

    TV face_sizes=grid.Face_Sizes();
    for(FACE_RANGE_ITERATOR<TV::m> it(grid.Domain_Indices());it.Valid();it.Next())
    {
        TV_INT i0=it.face.First_Cell_Index(),i1=it.face.Second_Cell_Index();
        int c0=c(i0),c1=c(i1),m0=c0>=0?color_map(c0):-1,m1=c1>=0?color_map(c1):-1;
        if(c0==c1) continue;
        T A=face_sizes(it.face.axis);
        if(dofs->V.Contains(it.face))
        {
            if(m0>=0) static_cast<V*>(array(m0))->V.Get_Or_Insert(it.face)+=A;
            if(m1>=0) static_cast<V*>(array(m1))->V.Get_Or_Insert(it.face)-=A;
        }
        else
        {
            T u=A*fv(it.face);
            if(m0>=0) rhs(m0)-=u;
            if(m1>=0) rhs(m1)+=u;
        }
    }

    auto* r=dynamic_cast<FLUID_REGIONS_PB<TV>*>(regions);
    if(r)
    {
        r->num_regions=num_regions;
        r->regions.Resize(c.domain,no_init);
        for(RANGE_ITERATOR<TV::m> it(c.domain);it.Valid();it.Next())
        {
            int col=c(it.index);
            r->regions(it.index)=col>=0?color_map(col):-1;
        }
    }
    
    delete dofs;
}

//#####################################################################
// Function Compute_Region_Mapping
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Compute_Region_Mapping(const FLUID_REGIONS<TV>* prev,
    const FLUID_REGIONS<TV>* next,ARRAY<int>& next_to_prev) const
{
    const auto* p=dynamic_cast<const FLUID_REGIONS_PB<TV>*>(prev);
    const auto* n=dynamic_cast<const FLUID_REGIONS_PB<TV>*>(next);
    next_to_prev.Resize(n->num_regions);
    ARRAY<ARRAY<int> > counts(n->num_regions);
    for(auto&a:counts) a.Resize(p->num_regions);
    for(RANGE_ITERATOR<TV::m> it(n->regions.domain);it.Valid();it.Next())
    {
        int a=n->regions(it.index);
        int b=p->regions(it.index);
        if(a>=0 && b>=0) counts(a)(b)++;
    }
    for(int i=0;i<counts.m;i++)
    {
        int num=counts(i).Sum();
        int j=counts(i).Arg_Max();
        if(counts(i)(j)*2>num) next_to_prev(i)=j;
        else next_to_prev(i)=-1;
    }
}

//#####################################################################
// Function Make_Boundary_Vector
//#####################################################################
template<class TV> FLUID_BOUNDARY_VECTOR<TV>* FLUID_SOLVER_PB<TV>::
Make_Boundary_Vector() const
{
    return new FLUID_BOUNDARY_VECTOR_PB<TV>;
}

//#####################################################################
// Function Make_Regions
//#####################################################################
template<class TV> FLUID_REGIONS<TV>* FLUID_SOLVER_PB<TV>::
Make_Regions() const
{
    return new FLUID_REGIONS_PB<TV>;
}

//#####################################################################
// Function Get_Force
//#####################################################################
template<class TV> void FLUID_SOLVER_PB<TV>::
Get_Force(FLUID_BOUNDARY_VECTOR<TV>* force) const
{
    const auto& p=driver->example.fluids_parameters.incompressible->projection.p;
    const auto& psi_D=driver->example.fluids_parameters.incompressible->projection.elliptic_solver->psi_D;
    const GRID<TV>& grid=*driver->example.fluids_parameters.grid;
    auto& V=dynamic_cast<FLUID_BOUNDARY_VECTOR_PB<TV>*>(force)->V;

    TV face_sizes=grid.Face_Sizes();
    for(auto&h:V)
    {
        TV_INT i0=h.key.First_Cell_Index(),i1=h.key.Second_Cell_Index();
        T x=0;
        if(!psi_D(i0)) x+=p(i0);
        if(!psi_D(i1)) x-=p(i1);
        h.data=x*face_sizes(h.key.axis);
    }
}

template class FLUID_SOLVER_PB<VECTOR<float,2> >;
template class FLUID_SOLVER_PB<VECTOR<float,3> >;
template class FLUID_SOLVER_PB<VECTOR<double,2> >;
template class FLUID_SOLVER_PB<VECTOR<double,3> >;
}
