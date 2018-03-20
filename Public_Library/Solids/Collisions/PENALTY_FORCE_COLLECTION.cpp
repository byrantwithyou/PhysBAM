//#####################################################################
// Copyright 2018, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Spatial_Acceleration/CONTINUOUS_COLLISION_DETECTION.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Forces_And_Torques/RIGID_PENALTY_WITH_FRICTION.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Forces/SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Collisions/PENALTY_FORCE_COLLECTION.h>
#include <Solids/Forces_And_Torques/RIGID_DEFORMABLE_PENALTY_WITH_FRICTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PENALTY_FORCE_COLLECTION<TV>::
~PENALTY_FORCE_COLLECTION()
{
    delete ccd_i;
    delete ccd_d;
    delete ccd_r;
}
//#####################################################################
// Function Update_Collision_Detection_Structures
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_Collision_Detection_Structures()
{
    bool new_grid=Update_Grid() || restarted;
    Update_Cell_Particles(new_grid);
    Update_Rasterized_Data(new_grid);
    Update_Cell_Vertices(new_grid);
    Update_Cell_Objects(new_grid);
    Update_CCD_Positions();
}
//#####################################################################
// Function Reset_Hash_Table
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Reset_Hash_Table()
{
    cell_particles.hash.Clean_Memory();
    cell_objects.hash.Clean_Memory();
    cell_vertices.hash.Clean_Memory();
    rasterized_data.hash.Clean_Memory();
}
//#####################################################################
// Function Update_Cell_Vertices
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_Cell_Vertices(bool new_grid)
{
    if(!rr_penalty || use_rr_ccd) return;
    cell_vertices.Remove_All();
    const RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    for(int b=0;b<rigid_body_particles.frame.m;b++){
        if(exclude_rigid_body_simplices.Contains(b)) continue;
        RIGID_BODY<TV>& rigid_body=*rigid_body_particles.rigid_body(b);
        for(int v=0;v<rigid_body.simplicial_object->particles.number;v++){
            TV X=rigid_body.Frame()*rigid_body.simplicial_object->particles.X(v);
            TV_INT index=grid.Cell(X);
            if(grid.Domain_Indices().Lazy_Inside_Half_Open(index))
                cell_vertices.Insert(index,{b,v});}}
}
//#####################################################################
// Function Update_Rasterized_Data
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_Rasterized_Data(bool new_grid)
{
    if((!rd_penalty || use_rd_ccd) && (!rr_penalty || use_rr_ccd)) return;
    rasterized_data.Remove_All();
    int clamp_ghost=1000000;
    const RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    T padding=grid.dX.Magnitude()*.5;
    for(int b=0;b<rigid_body_particles.frame.m;b++){
        RIGID_BODY<TV>& rigid_body=*rigid_body_particles.rigid_body(b);
        rigid_body.Update_Bounding_Box_From_Implicit_Geometry();
        RANGE<TV> box=rigid_body.Axis_Aligned_Bounding_Box().Thickened(padding);
        RANGE<TV_INT> grid_range=grid.Clamp_To_Cell(box.Intersect(grid.domain),clamp_ghost+1);
        for(RANGE_ITERATOR<TV::m> it(grid_range);it.Valid();it.Next()){
            TV X=grid.Center(it.index);
            T phi=rigid_body.Implicit_Geometry_Extended_Value(X);
            if(phi<padding)
                rasterized_data.Insert(it.index,{b,phi});}}
}
//#####################################################################
// Function Update_Cell_Particles
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_Cell_Particles(bool new_grid)
{
    if((!di_penalty || use_di_ccd) && (!rd_penalty || use_rd_ccd)) return;
    const DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    cell_particles.Remove_All();
    for(int k=0;k<simulated_particles.m;k++){
        int p=simulated_particles(k);
        TV_INT index=grid.Cell(particles.X(p));
        if(grid.Domain_Indices().Lazy_Inside_Half_Open(index))
            cell_particles.Insert(index,p);}
}
//#####################################################################
// Function Update_Cell_Objects
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_Cell_Objects(bool new_grid)
{
    if(!di_penalty || use_di_ccd) return;
    if(!new_grid) return;
    cell_objects.Remove_All();
    T padding=grid.dX.Magnitude()/2;
    int clamp_ghost=1000000;
    for(int o=0;o<di_penalty->ios.m;o++){
        IMPLICIT_OBJECT<TV>* io=di_penalty->ios(o);
        io->Update_Box();
        RANGE<TV> box=io->Box().Thickened(padding);
        RANGE<TV_INT> grid_range=grid.Clamp_To_Cell(box.Intersect(grid.domain),clamp_ghost+1);
        for(RANGE_ITERATOR<TV::m> it(grid_range);it.Valid();it.Next())
            if(io->Extended_Phi(grid.Center(it.index))<padding)
                cell_objects.Insert(it.index,o);}
}
//#####################################################################
// Function Get_DI_Collision_Candidates
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Get_DI_Collision_Candidates()
{
    if(use_di_ccd){
        ARRAY<CCD_PAIR<TV::m> > point_face;
        ccd_d->Compute_Pairs_PF(point_face,*ccd_i,const_repulsion_thickness);
        for(const auto& a:point_face)
            di_penalty->Add_Pair(a.f(0),a.s1);}
    else
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
    if(use_rd_ccd){
        ARRAY<CCD_PAIR<TV::m> > point_face;
        ccd_d->Compute_Pairs_PF(point_face,*ccd_r,const_repulsion_thickness);
        for(const auto& a:point_face)
            rd_penalty->Add_Pair(a.f(0),a.s1);}
    else
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
    if(use_rr_ccd){
        ARRAY<CCD_PAIR<TV::m> > point_face;
        ccd_r->Compute_Pairs_PF(point_face,const_repulsion_thickness);
        for(const auto& a:point_face)
            rr_penalty->Add_Pair(a.s0,a.f(0),a.s1);}
   else
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
    if(!state_saved) return;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    ARRAY<CCD_PAIR<TV::m> > point_face;
    ccd_d->Compute_Pairs_PF(point_face,const_repulsion_thickness);
    
    for(int i=0;i<point_face.m;i++){ // p f
        VECTOR<int,TV::m+1> pf=point_face(i).f;
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
        T_FACE face(X0.Subset(f));
        T collision_time=0;
        TV normal;
        VECTOR<T,TV::m+1> weights;
        VECTOR<TV,TV::m> V_f(particles.X.Subset(f)-X0.Subset(f));
        bool in=face.Point_Face_Collision(X0(p),particles.X(p)-X0(p),V_f,1,
            const_repulsion_thickness,collision_time,normal,weights,false);
        if(!in) continue;

        TV w=weights.Remove_Index(0);
        if(w.Min()<0) continue;
        dd_penalty->Add_Pair(p,s_e.x,w,s_e.y);}
}
//#####################################################################
// Function Save_State
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Save_State()
{
    const DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    Update_CCD_Topology();
    if(dd_penalty || ccd_r){
        state_saved=true;
        X0=particles.X;
        F0=solid_body_collection.rigid_body_collection.rigid_body_particles.frame;}
}
//#####################################################################
// Function Init
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Init(bool use_di,bool use_dd,bool use_rd,bool use_rr)
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    if(use_rr && !rr_penalty && solid_body_collection.rigid_body_collection.rigid_body_particles.number>0){
        rr_penalty=new RIGID_PENALTY_WITH_FRICTION<TV>(
            solid_body_collection.rigid_body_collection,move_rb_diff);
        rr_penalty->get_candidates=[this](){Get_RR_Collision_Candidates();};
        solid_body_collection.Add_Force(rr_penalty);}

    if(use_rd && !rd_penalty && solid_body_collection.rigid_body_collection.rigid_body_particles.number>0){
        rd_penalty=new RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>(
            solid_body_collection.deformable_body_collection.particles,
            solid_body_collection.rigid_body_collection,move_rb_diff);
        rd_penalty->get_candidates=[this](){Get_RD_Collision_Candidates();};
        solid_body_collection.Add_Force(rd_penalty);}

    if(use_di && !di_penalty){
        di_penalty=new IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>(
            solid_body_collection.deformable_body_collection.particles);
        di_penalty->get_candidates=[this](){Get_DI_Collision_Candidates();};
        solid_body_collection.deformable_body_collection.Add_Force(di_penalty);}

    if(use_dd){
        dd_penalty=new SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>(particles);
        dd_penalty->get_candidates=[this](){Get_DD_Collision_Candidates();};
        solid_body_collection.Add_Force(dd_penalty);
        for(int i=0;i<deformable_body_collection.structures.m;i++){
            auto* s=deformable_body_collection.structures(i);
            if(auto* p=dynamic_cast<T_SURFACE*>(s))
                dd_penalty->Add_Surface(*p);
            else if(auto* p=dynamic_cast<T_OBJECT*>(s))
                dd_penalty->Add_Surface(p->Get_Boundary_Object());}}

    Update_CCD_Topology();
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
//#####################################################################
// Function Update_Grid
//#####################################################################
template<class TV> bool PENALTY_FORCE_COLLECTION<TV>::
Update_Grid()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    const RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    RANGE<TV> bounding_box=RANGE<TV>::Bounding_Box(particles.X);
    T max_volume=0;
    
    for(int b=0;b<rigid_body_particles.frame.m;b++){
        if(exclude_rigid_body_simplices.Contains(b)) continue;
        RIGID_BODY<TV>& rigid_body=*rigid_body_particles.rigid_body(b);
        rigid_body.Update_Bounding_Box();
        RANGE<TV> box=rigid_body.Axis_Aligned_Bounding_Box();
        bounding_box.Enlarge_To_Include_Box(box);
        max_volume=std::max(max_volume,bounding_box.Size());}

    // Existing box is good enough.
    if(grid.domain.Contains(bounding_box)) return false;

    T max_dx=bounding_box.Edge_Lengths().Max()/max_resolution;
    max_dx=std::max(max_dx,pow<1,TV::m>(max_volume/max_cells_per_object));
    max_dx=std::max(max_dx,pow<1,TV::m>(bounding_box.Size()/max_cells));
    grid=GRID<TV>::Create_Grid_Given_Cell_Size(bounding_box,max_dx,true,0);
    return true;
}
//#####################################################################
// Function Update_CCD_Topology
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_CCD_Topology()
{
    if(use_di_ccd)
        if(!ccd_i || ccd_i_stale){
            if(!ccd_i) ccd_i=new CONTINUOUS_COLLISION_DETECTION<TV>;
            else ccd_i->Clean_Memory();
            for(int i=0;i<di_penalty->meshes.m;i++)
                ccd_i->Add_Structure(di_penalty->meshes(i),4); // face only
            ccd_i_stale=false;}

    if(use_di_ccd || use_rd_ccd || dd_penalty)
        if(!ccd_d || ccd_d_stale){
            const auto& st=solid_body_collection.deformable_body_collection.structures;
            const DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
            if(!ccd_d) ccd_d=new CONTINUOUS_COLLISION_DETECTION<TV>;
            else ccd_d->Clean_Memory();
            for(int i=0;i<st.m;i++)
                ccd_d->Add_Structure(st(i),4); // face only
            ccd_d->Add_Particles(particles.number);
            ccd_d_stale=false;}

    if(use_rr_ccd || use_rd_ccd)
        if(!ccd_r || ccd_r_stale){
            const auto& rbp=solid_body_collection.rigid_body_collection.rigid_body_particles;
            if(!ccd_r) ccd_r=new CONTINUOUS_COLLISION_DETECTION<TV>;
            else ccd_r->Clean_Memory();
            for(int b=0;b<rbp.frame.m;b++){
                if(exclude_rigid_body_simplices.Contains(b)) continue;
                // face+vertex only
                ccd_r->Add_Structure(rbp.rigid_body(b)->simplicial_object,5);}
            ccd_r_stale=false;}
}
//#####################################################################
// Function Update_CCD_Positions
//#####################################################################
template<class TV> void PENALTY_FORCE_COLLECTION<TV>::
Update_CCD_Positions()
{
    if(ccd_i)
        for(int i=0;i<di_penalty->meshes.m;i++)
            ccd_i->Update_Positions(di_penalty->meshes(i)->particles.X,i);

    if(ccd_d)
        ccd_d->Update_Positions(X0,solid_body_collection.deformable_body_collection.particles.X);

    if(ccd_r){
        const auto& rbp=solid_body_collection.rigid_body_collection.rigid_body_particles;
        for(int b=0;b<rbp.frame.m;b++)
            ccd_r->Update_Positions(rbp.rigid_body(b)->simplicial_object->particles.X,
                F0(b),rbp.frame(b),b);}
}
template class PENALTY_FORCE_COLLECTION<VECTOR<float,2> >;
template class PENALTY_FORCE_COLLECTION<VECTOR<float,3> >;
template class PENALTY_FORCE_COLLECTION<VECTOR<double,2> >;
template class PENALTY_FORCE_COLLECTION<VECTOR<double,3> >;
}
