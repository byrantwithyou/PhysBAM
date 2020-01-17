//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Vectors/TWIST.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Projection/BOUNDARY_CONDITION_DOUBLE_FINE.h>
#include <Geometry/Topology_Based_Geometry/SIMPLEX_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include "FLUID_BOUNDARY_VECTOR_PB.h"
#include "FLUID_SOLVER_PB.h"
#include "SOLID_BOUNDARY_VECTOR_PB.h"
#include "SOLID_FLUID_INTERFACE_PB.h"
#include "SOLID_SOLVER_PB.h"
namespace PhysBAM{

//#####################################################################
// Function Compute_BC
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Compute_BC(const FLUID_SOLVER<TV>* fluid_solver,SOLID_BC<TV>* solid_bc,T time,T dt) const
{
    LOG::printf("DO THIS: %s\n",__PRETTY_FUNCTION__);
    // TODO go away?
}

//#####################################################################
// Function Compute_BC
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Compute_BC(const SOLID_SOLVER<TV>* solid_solver,FLUID_BC<TV>* fluid_bc,T time,T dt) const
{
    LOG::printf("DO THIS: %s\n",__PRETTY_FUNCTION__);
    // TODO go away?
}

//#####################################################################
// Function Interpolate_Velocity
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Interpolate_Velocity(FLUID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v)
{
    auto* v_pb=dynamic_cast<const SOLID_BOUNDARY_VECTOR_PB<TV>*>(v);
    auto* u_pb=dynamic_cast<FLUID_BOUNDARY_VECTOR_PB<TV>*>(u);
    u_pb->V.Remove_All();
    for(const auto& e:deformable_weights) u_pb->V.Get_Or_Insert(e.face)+=e.w*v_pb->V.Get(e.p)(e.face.axis);
    for(const auto& e:rigid_weights)
    {
        TWIST<TV> t=v_pb->twist.Get(e.p);
        u_pb->V.Get_Or_Insert(e.face)+=e.w*t.linear(e.face.axis)+e.aw.Dot(t.angular);
    }
}

//#####################################################################
// Function Distribute_Force
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Distribute_Force(SOLID_BOUNDARY_VECTOR<TV>* v, const FLUID_BOUNDARY_VECTOR<TV>* u)
{
    auto* v_pb=dynamic_cast<SOLID_BOUNDARY_VECTOR_PB<TV>*>(v);
    auto* u_pb=dynamic_cast<const FLUID_BOUNDARY_VECTOR_PB<TV>*>(u);
    v_pb->V.Remove_All();
    v_pb->twist.Remove_All();
    for(const auto& e:deformable_weights) v_pb->V.Get_Or_Insert(e.p)(e.face.axis)+=e.w*u_pb->V.Get(e.face);
    for(const auto& e:rigid_weights)
    {
        TWIST<TV>& t=v_pb->twist.Get_Or_Insert(e.p);
        T a=u_pb->V.Get(e.face);
        t.linear(e.face.axis)+=e.w*a;
        t.angular+=e.aw*a;
    }
}

//#####################################################################
// Function Get_Boundary
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Get_Boundary(SOLID_BOUNDARY_VECTOR<TV>* v) const
{
    auto* v_pb=dynamic_cast<SOLID_BOUNDARY_VECTOR_PB<TV>*>(v);
    v_pb->V.Remove_All();
    v_pb->twist.Remove_All();
    for(const auto& e:deformable_weights) v_pb->V.Set(e.p,TV());
    for(const auto& e:rigid_weights) v_pb->twist.Set(e.p,TWIST<TV>());
}

//#####################################################################
// Function Get_Boundary
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Get_Boundary(FLUID_BOUNDARY_VECTOR<TV>* v) const
{
    auto* v_pb=dynamic_cast<FLUID_BOUNDARY_VECTOR_PB<TV>*>(v);
    v_pb->V.Remove_All();
    for(const auto& e:deformable_weights) v_pb->V.Set(e.face,0);
    for(const auto& e:rigid_weights) v_pb->V.Set(e.face,0);
}

//#####################################################################
// Function Compute_Coupling_Weights
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Compute_Coupling_Weights(const SOLID_SOLVER<TV>* solid_solver,const FLUID_SOLVER<TV>* fluid_solver,T time,T dt)
{
    auto* ss=dynamic_cast<const SOLID_SOLVER_PB<TV>*>(solid_solver);
    auto* fs=dynamic_cast<const FLUID_SOLVER_PB<TV>*>(fluid_solver);

    auto* bc_fine=fs->driver->example.fluids_parameters.bc_fine;
    bc_fine->Reset(0);

    if(fs->driver->example.get_unified_boundary_conditions)
        fs->driver->example.get_unified_boundary_conditions(bc_fine,time);

    auto& sbc=ss->driver->example.solid_body_collection;
    auto& rbc=sbc.rigid_body_collection;
    auto& dbc=sbc.deformable_body_collection;

    const GRID<TV>& grid=bc_fine->mac_grid;
    VECTOR<ARRAY<bool,TV_INT>,TV::m> psi_D;
    VECTOR<ARRAY<bool,FACE_INDEX<TV::m> >,TV::m> psi_N;
    VECTOR<ARRAY<T,TV_INT>,TV::m> psi_u;
    for(int i=0;i<TV::m;i++)
    {
        GRID<TV> face_grid=grid.Get_Face_Grid(i);
        psi_D(i).Resize(face_grid.Domain_Indices(bc_fine->ghost));
        psi_N(i).Resize(face_grid.Domain_Indices(bc_fine->ghost));
        psi_u(i).Resize(face_grid.Domain_Indices(bc_fine->ghost));
        bc_fine->Get_Viscosity_Boundary_Conditions(psi_D(i),psi_N(i),psi_u(i),i);
    }

    auto coupled_face=[&psi_D,bc_fine](const FACE_INDEX<TV::m>& face)
        {
            TV_INT fine_node_index=(2*face.index+1).Add_Axis(face.axis,-1);
            if(bc_fine->bc_type(fine_node_index.Add_Axis(face.axis,-1))>=0) return true;
            if(bc_fine->bc_type(fine_node_index.Add_Axis(face.axis,1))>=0) return true;
            for(int a=0;a<TV::m;a++)
                for(int s=-1;s<2;s+=2)
                    if(!psi_D(face.axis)(face.index.Add_Axis(a,s)))
                        return true;
            return false;
        };
    
    for(auto* rb:rbc.rigid_body_particles.rigid_body)
    {
        if(!rb) continue;
        auto frame=rb->Frame();
        auto cb=[this,rb,frame,coupled_face](const auto& data)
            {
                RIGID_ENTRY re=
                {
                    rb->particle_index,
                    data.face,
                    1,
                    (data.X-frame.t).Cross(TV::Axis_Vector(data.face.axis))
                };
                rigid_weights.Append(re);
                return 0;
            };
        if(rb->simplicial_object)
        {
            SIMPLEX_OBJECT<TV,TV::m-1> copy;
            copy.mesh.elements=rb->simplicial_object->mesh.elements;
            copy.particles.Add_Elements(rb->simplicial_object->particles.X.m);
            auto A=rb->simplicial_object->particles.X;
            auto B=copy.particles.X;
            for(int i=0;i<A.m;i++) B(i)=frame*A(i);
            bc_fine->Set(copy,bc_fine->bc_noslip,cb,rb->thin_shell);
        }
        else if(rb->implicit_object)
            bc_fine->Set(rb->implicit_object,bc_fine->bc_noslip,cb,rb->thin_shell);
    }

    T kernel_radius=grid.dX.Max()*2;
    auto kernel=[kernel_radius](const TV& dX)
    {
        T d=dX.Magnitude();
        if(d>kernel_radius) return (T)0;
        return 1-d/kernel_radius;
    };

    int s=0;
    for(auto* st:dbc.structures)
    {
        SIMPLEX_OBJECT<TV,TV::m-1>* surface=0;
        if(auto* o=dynamic_cast<SIMPLEX_OBJECT<TV,TV::m-1>*>(st)) surface=o;
        else if(auto* o=dynamic_cast<SIMPLEX_OBJECT<TV,TV::m>*>(st))
            surface=&o->Get_Boundary_Object();
        if(!surface) continue;
        bool thin=s<structure_is_thin.m && structure_is_thin(s);
        HASHTABLE<FACE_INDEX<TV::m>,T> hash;
        bc_fine->Set(*surface,bc_fine->bc_noslip,
            [this,&hash](const auto& data)
            {
                hash.Set(data.face,0);
                return 0;
            },thin);
        
        ARRAY<int> particle_list;
        Get_Unique(particle_list,surface->mesh.elements.Flattened());
        int first_weight=deformable_weights.m;
        for(int p:particle_list)
        {
            TV X=surface->particles.X(p);
            auto box=RANGE<TV>(X).Thickened(kernel_radius);
            box=box.Intersect(grid.domain);
            auto grid_box=grid.Clamp_To_Cell(box,0);
            for(FACE_RANGE_ITERATOR<TV::m> it(grid_box);it.Valid();it.Next())
                if(T* total=hash.Get_Pointer(it.face))
                    if(T w=kernel(grid.Face(it.face)-X))
                    {
                        deformable_weights.Append({p,it.face,w});
                        *total+=w;
                    }
        }
        for(int i=first_weight;i<deformable_weights.m;i++)
            deformable_weights(i).w/=hash.Get(deformable_weights(i).face);
        
        s++;
    }

    HASHTABLE<FACE_INDEX<TV::m>,T> face_weight;
    for(int i=0;i<rigid_weights.m;i++)
    {
        if(!coupled_face(rigid_weights(i).face))
            rigid_weights.Remove_Index_Lazy(i--);
        else face_weight.Get_Or_Insert(rigid_weights(i).face)+=1;
    }

    for(int i=0;i<deformable_weights.m;i++)
    {
        if(!coupled_face(deformable_weights(i).face))
            deformable_weights.Remove_Index_Lazy(i--);
        else face_weight.Get_Or_Insert(deformable_weights(i).face)+=1;
    }

    for(auto& rw:rigid_weights)
    {
        T tot=face_weight.Get(rw.face);
        rw.w/=tot;
        rw.aw/=tot;
    }

    for(auto& dw:deformable_weights)
        dw.w/=face_weight.Get(dw.face);
}

template class SOLID_FLUID_INTERFACE_PB<VECTOR<float,2> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<float,3> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<double,2> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<double,3> >;
}
