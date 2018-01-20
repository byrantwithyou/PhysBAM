//#####################################################################
// Copyright 2018, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Forces_And_Torques/RIGID_PENALTY_WITH_FRICTION.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Forces/SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Collisions/PENALTY_FORCE_COLLECTION.h>
#include <Solids/Forces_And_Torques/RIGID_DEFORMABLE_PENALTY_WITH_FRICTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
namespace PhysBAM{
//#####################################################################
// Function Update_Collision_Detection_Structures
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_Collision_Detection_Structures()
{
    cell_particles.Remove_All();
    cell_vertices.Remove_All();
    rasterized_data.Remove_All();
    int clamp_ghost=1000000;
    
    const DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    if(di_penalty || rd_penalty)
        for(int k=0;k<simulated_particles.m;k++){
            int p=simulated_particles(k);
            TV_INT index=grid.Cell(particles.X(p));
            if(domain_of_interest.Lazy_Inside_Half_Open(index))
                cell_particles.Insert(index,p);}

    const RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    if(rd_penalty || rr_penalty){
        T padding=grid.dX.Max()*.5;
        for(int b=0;b<rigid_body_particles.frame.m;b++){
            RIGID_BODY<TV>& rigid_body=*rigid_body_particles.rigid_body(b);
            rigid_body.Update_Bounding_Box_From_Implicit_Geometry();
            RANGE<TV> box=rigid_body.Axis_Aligned_Bounding_Box().Thickened(padding);
            RANGE<TV_INT> grid_range=grid.Clamp_To_Cell(box,clamp_ghost+1).Intersect(domain_of_interest);
            for(RANGE_ITERATOR<TV::m> it(grid_range);it.Valid();it.Next()){
                TV X=grid.Center(it.index);
                T phi=rigid_body.Implicit_Geometry_Extended_Value(X);
                if(phi<padding)
                    rasterized_data.Insert(it.index,{b,phi});}}}

    if(rr_penalty)
        for(int b=0;b<rigid_body_particles.frame.m;b++){
            RIGID_BODY<TV>& rigid_body=*rigid_body_particles.rigid_body(b);
            for(int v=0;v<rigid_body.simplicial_object->particles.number;v++){
                TV X=rigid_body.Frame()*rigid_body.simplicial_object->particles.X(v);
                TV_INT index=grid.Cell(X);
                if(domain_of_interest.Lazy_Inside_Half_Open(index))
                    cell_vertices.Insert(index,{b,v});}}
}
//#####################################################################
// Function Get_DI_Collision_Candidates
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Get_DI_Collision_Candidates()
{
    for(const auto& k:cell_objects.hash){
        auto D=cell_particles.Get(k.key);
        if(!D.m) continue;
        auto O=cell_objects.Get(k.key);
        for(int i=0;i<D.m;i++)
            for(int j=0;j<O.m;j++)
                di_penalty->Add_Pair(D(i),O(j));}
}
//#####################################################################
// Function Get_RD_Collision_Candidates
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Get_RD_Collision_Candidates()
{
    for(const auto& k:rasterized_data.hash){
        auto D=cell_particles.Get(k.key);
        if(!D.m) continue;
        auto R=rasterized_data.Get(k.key);
        for(int i=0;i<D.m;i++)
            for(int j=0;j<R.m;j++)
                rd_penalty->Add_Pair(D(i),R(j).id);}
}
//#####################################################################
// Function Get_RR_Collision_Candidates
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Get_RR_Collision_Candidates()
{
    for(const auto& k:rasterized_data.hash){
        auto V=cell_vertices.Get(k.key);
        if(!V.m) continue;
        auto R=rasterized_data.Get(k.key);
        for(int i=0;i<V.m;i++)
            for(int j=0;j<R.m;j++)
                if(V(i).x!=R(j).id)
                    rr_penalty->Add_Pair(V(i).x,V(i).y,R(j).id);}
}
//#####################################################################
// Function Get_DD_Collision_Candidates
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Get_DD_Collision_Candidates()
{
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    ARRAY<TV>& Xn=deformable_body_collection.triangle_repulsions_and_collisions_geometry.X_self_collision_free;
    deformable_body_collection.triangle_collisions.Update_Swept_Hierachies_And_Compute_Pairs(
        particles.X,Xn,recently_modified,const_repulsion_thickness);

    LOG::printf("num candidate pairs: %i\n",deformable_body_collection.triangle_collisions.point_face_pairs_internal.m);
    for(int i=0;i<deformable_body_collection.triangle_collisions.point_face_pairs_internal.m;i++){ // p f
        VECTOR<int,TV::m+1> pf=deformable_body_collection.triangle_collisions.point_face_pairs_internal(i);
        int p=pf(0);
        TV_INT f=pf.Remove_Index(0);

        // Detection is expensive; make sure we are not already known.
        auto s_e=dd_penalty->object_from_element.Get(f.Sorted());
        if(dd_penalty->hash.Contains({p,s_e.x})) continue;
        const auto& ts=*dd_penalty->surfaces(s_e.x);
        PHYSBAM_ASSERT(f==ts.mesh.elements(s_e.y));

        // Particle must have exited.
        if(ts.Get_Element(s_e.y).Signed_Distance(particles.X(p))>0) continue;

        // Do the expensive check.
        T_FACE face(Xn.Subset(f));
        T collision_time=0;
        TV normal;
        VECTOR<T,TV::m+1> weights;
        VECTOR<TV,TV::m> V_f(particles.X.Subset(f)-Xn.Subset(f));
        bool in=face.Point_Face_Collision(Xn(p),particles.X(p)-Xn(p),V_f,1,
            const_repulsion_thickness,collision_time,normal,weights,false);
        if(!in) continue;
        
        dd_penalty->Add_Pair(p,s_e.x,weights.Remove_Index(0),s_e.y);}
}
//#####################################################################
// Function Save_State
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Save_State()
{
    const DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    if(dd_penalty){
        solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.X_self_collision_free=particles.X;
        repulsion_thickness.Resize(particles.number,true,true,const_repulsion_thickness);
        recently_modified.Resize(particles.number);
        recently_modified.Fill(true);}
}
//#####################################################################
// Function Init
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Init(T stiffness,T friction,TRIANGLE_COLLISION_PARAMETERS<TV>* param,
    bool use_di,bool use_dd,bool use_rd,bool use_rr)
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    if(use_rr && !rr_penalty && solid_body_collection.rigid_body_collection.rigid_body_particles.number>0){
        rr_penalty=new RIGID_PENALTY_WITH_FRICTION<TV>(
            solid_body_collection.rigid_body_collection,move_rb_diff,
            stiffness,friction);
        rr_penalty->get_candidates=[this](){Get_RR_Collision_Candidates();};
        solid_body_collection.Add_Force(rr_penalty);}

    if(use_rd && !rd_penalty && solid_body_collection.rigid_body_collection.rigid_body_particles.number>0){
        rd_penalty=new RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>(
            solid_body_collection.deformable_body_collection.particles,
            solid_body_collection.rigid_body_collection,move_rb_diff,
            stiffness,friction);
        rd_penalty->get_candidates=[this](){Get_RD_Collision_Candidates();};
        solid_body_collection.Add_Force(rd_penalty);}

    if(use_di && !di_penalty){
        di_penalty=new IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>(
            solid_body_collection.deformable_body_collection.particles,
            stiffness,friction);
        di_penalty->get_candidates=[this](){Get_DI_Collision_Candidates();};
        solid_body_collection.deformable_body_collection.Add_Force(di_penalty);}

    if(use_dd){
        dd_penalty=new SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>(
            particles,stiffness,friction);
        dd_penalty->get_candidates=[this](){Get_DD_Collision_Candidates();};
        solid_body_collection.Add_Force(dd_penalty);
        deformable_body_collection.triangle_collisions.compute_edge_edge_collisions=false;
        for(int i=0;i<deformable_body_collection.structures.m;i++){
            auto* s=deformable_body_collection.structures(i);
            if(auto* p=dynamic_cast<T_SURFACE*>(s))
                dd_penalty->Add_Surface(*p);
            else if(auto* p=dynamic_cast<T_OBJECT*>(s))
                dd_penalty->Add_Surface(p->Get_Boundary_Object());}
        deformable_body_collection.triangle_repulsions_and_collisions_geometry.Build_Collision_Geometry();
        deformable_body_collection.triangle_repulsions_and_collisions_geometry.X_self_collision_free=particles.X;
        repulsion_thickness.Resize(particles.number,true,true,const_repulsion_thickness);
        recently_modified.Resize(particles.number,true,true,true);

        PHYSBAM_ASSERT(param);
        deformable_body_collection.triangle_repulsions_and_collisions_geometry.Initialize(*param);
        deformable_body_collection.triangle_repulsions.Initialize(*param);
        deformable_body_collection.triangle_collisions.Initialize(*param);}
}
//#####################################################################
// Function Rasterize_Implicit_Object
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Rasterize_Implicit_Object(IMPLICIT_OBJECT<TV>* io)
{
    T padding=grid.dX.Magnitude()/2;
    int o=di_penalty->ios.Append(io);
    for(RANGE_ITERATOR<TV::m> it(domain_of_interest);it.Valid();it.Next())
        if(io->Extended_Phi(grid.Center(it.index))<padding)
            cell_objects.Insert(it.index,o);
}
//#####################################################################
// Function Update_Attachments_And_Prune_Pairs
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_Attachments_And_Prune_Pairs()
{
    if(rd_penalty) rd_penalty->Update_Attachments_And_Prune_Pairs();
    if(di_penalty) di_penalty->Update_Attachments_And_Prune_Pairs();
    if(dd_penalty) dd_penalty->Update_Attachments_And_Prune_Pairs();
    if(rr_penalty) rr_penalty->Update_Attachments_And_Prune_Pairs();
}
template class PENALTY_FORCE_COLLECTION<VECTOR<float,2> >;
template class PENALTY_FORCE_COLLECTION<VECTOR<float,3> >;
template class PENALTY_FORCE_COLLECTION<VECTOR<double,2> >;
template class PENALTY_FORCE_COLLECTION<VECTOR<double,3> >;
}
