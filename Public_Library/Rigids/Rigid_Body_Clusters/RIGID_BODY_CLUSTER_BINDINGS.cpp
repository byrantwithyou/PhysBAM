//#####################################################################
// Copyright 2008-2009, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_CLUSTER_BINDINGS
//##################################################################### 
#include <Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <Tools/Data_Structures/TRIPLE.h>
#include <Tools/Log/DEBUG_PRINT.h>
#include <Tools/Matrices/MATRIX.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <Rigids/Joints/JOINT.h>
#include <Rigids/Joints/JOINT_MESH.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
using namespace PhysBAM;
//#####################################################################
// Function RIGID_BODY_CLUSTER_BINDINGS
//#####################################################################
template<class TV> RIGID_BODY_CLUSTER_BINDINGS<TV>::
RIGID_BODY_CLUSTER_BINDINGS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body_input)
    :callbacks(0),rigid_body_collection(rigid_body_collection_input),articulated_rigid_body(articulated_rigid_body_input),
    collide_constituent_bodies(false),saved(false)
{}
//#####################################################################
// Function Update_Joint_Structures
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Update_All_Joint_Structures()
{
    for(T_REVERSE_BINDING_ITERATOR i(reverse_bindings);i.Valid();i.Next())
        Update_Joint_Structures(i.Key());
}
//#####################################################################
// Function Update_Joint_Structures
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Update_Joint_Structures(const int parent) 
{
    CLUSTER& cluster=*reverse_bindings.Get(parent);
    HASHTABLE<int> hashtable;
    articulated_rigid_body.joint_mesh.undirected_graph.Ensure_Number_Nodes(rigid_body_collection.rigid_body_particles.Size());
    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++)
        hashtable.Set(cluster.children(i));
    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++){
        int id=cluster.children(i);
        const ARRAY<JOINT_ID>& list=articulated_rigid_body.joint_mesh.undirected_graph.Adjacent_Edges(id);
        for(int j=0;j<list.m;j++){
            PAIR<int,int> pair=articulated_rigid_body.joint_mesh.undirected_graph.Edges(list(j));
            int oid=pair.x+(Value(pair.y)-Value(id));
            if(!hashtable.Contains(oid)){
                BOUNDARY_JOINT boundary_joint;
                boundary_joint.joint_id=list(j);
                boundary_joint.boundary_body=id;
                boundary_joint.is_parent=(pair.x==id);
                if(boundary_joint.is_parent) boundary_joint.original_frame=articulated_rigid_body.joint_mesh(list(j))->F_pj();
                else boundary_joint.original_frame=articulated_rigid_body.joint_mesh(list(j))->F_cj();
                cluster.boundary_joints.Append(boundary_joint);}
            else if(id<oid){
                INTERNAL_JOINT internal_joint;
                internal_joint.joint=articulated_rigid_body.joint_mesh(list(j));
                internal_joint.parent=pair.x;
                internal_joint.child=pair.y;
                cluster.internal_joints.Append(internal_joint);}}}
}
//#####################################################################
// Function Append_Binding
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Append_Binding(const int particle,const PAIR<int,RIGID_CLUSTER_CONSTITUENT_ID> data)
{
    int size=binding_index(particle).m+1;
    binding_index(particle).Resize(size,true,true);
    binding_index(particle).Last()=binding_index(particle)(0);
    binding_index(particle)(0)=data;
}
//#####################################################################
// Function Bind_Particle
//#####################################################################
template<class TV> int RIGID_BODY_CLUSTER_BINDINGS<TV>::
Add_Binding(const ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID>& child_particles)
{
    RIGID_BODY<TV>* parent_body=new RIGID_BODY<TV>(rigid_body_collection,true);
    int parent=parent_body->particle_index;
    CLUSTER& cluster=*reverse_bindings.Get_Or_Insert(parent,new CLUSTER());
    cluster.parent=parent;
    cluster.stored_active=false;
    binding_index.Resize(rigid_body_collection.rigid_body_particles.Size());
    bool has_static=false,has_kinematic=false;
    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<child_particles.Size();i++){
        cluster.children.Append(child_particles(i));
        cluster.child_to_parent.Append(FRAME<TV>());
        RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(child_particles(i));
        if(body.Has_Infinite_Inertia()){
            if(body.is_static){assert(!has_kinematic);has_static=true;}
            else if(rigid_body_collection.rigid_body_particles.kinematic(child_particles(i))){PHYSBAM_ASSERT(!has_kinematic && !has_static);has_kinematic=true;}
            cluster.infinite_body=child_particles(i);}
        Append_Binding(child_particles(i),PAIR<int,RIGID_CLUSTER_CONSTITUENT_ID>(parent,i));}
    binding_index(parent).Append(PAIR<int,RIGID_CLUSTER_CONSTITUENT_ID>(parent,RIGID_CLUSTER_CONSTITUENT_ID(0)));
    Update_Joint_Structures(parent);
    Build_Aggregate_Geometry(parent);
    Set_Binding_Active(parent,true);
    if(callbacks) callbacks->Create_Cluster(parent);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(parent_body);
    return parent;
}
//#####################################################################
// Function Delete_Binding
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Delete_Binding(const int parent_particle)
{
    Set_Binding_Active(parent_particle,false);
    CLUSTER& cluster=*reverse_bindings.Get(parent_particle);
    rigid_body_collection.Reactivate_Body(cluster.parent);
    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++){
        int child=cluster.children(i);
        for(int j=0;j<binding_index(child).m;j++) if(binding_index(child)(j).x==parent_particle){binding_index(child).Remove_Index_Lazy(j);break;}}
    binding_index(parent_particle).Remove_All();
    reverse_bindings.Delete(parent_particle);
    rigid_body_collection.rigid_body_particles.Remove_Body(parent_particle);
    rigid_body_collection.Destroy_Unreferenced_Geometry();
    if(callbacks) callbacks->Destroy_Cluster(parent_particle);
    delete &cluster;
    // TODO: delete geometry
}
//#####################################################################
// Function Make_Active_Parent
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Make_Active_Parent(const int parent_particle,ARRAY<PAIR<int,RIGID_CLUSTER_CONSTITUENT_ID> >& child_list)
{
    for(int i=0;i<child_list.m;i++) if(child_list(i).x==parent_particle){exchange(child_list(0),child_list(i));return;}
}
//#####################################################################
// Function Set_Binding_Active
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Set_Binding_Active(const int parent_particle,const bool active)
{
    CLUSTER& cluster=*reverse_bindings.Get(parent_particle);
    if(cluster.active==active) return;
    if(active){
        if(!rigid_body_collection.Is_Active(cluster.parent)){
            rigid_body_collection.Reactivate_Body(cluster.parent);
            rigid_body_collection.Rigid_Body(cluster.parent).bounding_box_up_to_date=false;}
        Distribute_Mass_To_Parent(parent_particle);
        Update_Aggregate_Geometry(parent_particle);
        for(int i=0;i<cluster.internal_joints.m;i++)
            articulated_rigid_body.joint_mesh.Deactivate_Articulation(cluster.internal_joints(i).joint->id_number);
        for(int i=0;i<cluster.boundary_joints.m;i++){
            if(cluster.boundary_joints(i).is_parent){articulated_rigid_body.Substitute_Joint_Parent_Body(cluster.boundary_joints(i).joint_id,parent_particle);}
            else{articulated_rigid_body.Substitute_Joint_Child_Body(cluster.boundary_joints(i).joint_id,parent_particle);}}}
    else{
        Clamp_Particles_To_Embedded_Positions(parent_particle);
        Clamp_Particles_To_Embedded_Velocities(parent_particle);
        for(int i=0;i<cluster.internal_joints.m;i++)
            articulated_rigid_body.joint_mesh.Reactivate_Articulation(cluster.internal_joints(i).parent,cluster.internal_joints(i).child,cluster.internal_joints(i).joint);
        for(int i=0;i<cluster.boundary_joints.m;i++){
            if(cluster.boundary_joints(i).is_parent) articulated_rigid_body.Substitute_Joint_Parent_Body(cluster.boundary_joints(i).joint_id,cluster.boundary_joints(i).boundary_body,cluster.boundary_joints(i).original_frame);
            else articulated_rigid_body.Substitute_Joint_Child_Body(cluster.boundary_joints(i).joint_id,cluster.boundary_joints(i).boundary_body,cluster.boundary_joints(i).original_frame);}
        rigid_body_collection.Deactivate_Body(cluster.parent);}

    cluster.active=active;
    if(!active){
        rigid_body_collection.rigid_body_particles.frame(parent_particle).t=TV::All_Ones_Vector()*1e10;
        rigid_body_collection.rigid_body_particles.twist(parent_particle).linear=TV::All_Ones_Vector()*1e10;
        rigid_body_collection.rigid_body_particles.twist(parent_particle).angular=T_SPIN::All_Ones_Vector()*1e10;}
    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++) if(active) Make_Active_Parent(parent_particle,binding_index(cluster.children(i)));
}
//#####################################################################
// Function Deactivate_And_Return_Clusters
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Deactivate_And_Return_Clusters(ARRAY<int>& active_bindings)
{
    for(typename HASHTABLE<int,CLUSTER*>::ITERATOR i(reverse_bindings);i.Valid();i.Next()){
        if(i.Data()->active){active_bindings.Append(i.Key());Set_Binding_Active(i.Key(),false);}}
}
//#####################################################################
// Function Reactivate_Bindings
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Reactivate_Bindings(const ARRAY<int>& active_bindings)
{
    for(int i=0;i<active_bindings.m;i++) Set_Binding_Active(active_bindings(i),true);
}
//#####################################################################
// Function Save_Bindings_State
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Save_Bindings_State()
{
    for(typename HASHTABLE<int,CLUSTER*>::ITERATOR i(reverse_bindings);i.Valid();i.Next()){
        if(!i.Data()->active) continue;
        RIGID_BODY<TV>& parent_body=rigid_body_collection.Rigid_Body(i.Key());
        i.Data()->saved_mass=parent_body.Mass();
        i.Data()->saved_inertia_tensor=parent_body.Inertia_Tensor();
        i.Data()->saved_child_to_parent=i.Data()->child_to_parent;}
}
//#####################################################################
// Function Restore_Bindings_State
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Restore_Bindings_State()
{
    for(typename HASHTABLE<int,CLUSTER*>::ITERATOR i(reverse_bindings);i.Valid();i.Next()){
        if(!i.Data()->active) continue;
        RIGID_BODY<TV>& parent_body=rigid_body_collection.Rigid_Body(i.Key());
        parent_body.Mass()=i.Data()->saved_mass;
        parent_body.Inertia_Tensor()=i.Data()->saved_inertia_tensor;
        i.Data()->child_to_parent=i.Data()->saved_child_to_parent;}
}
//#####################################################################]
// Get_Child_To_Parent_Frame
//#####################################################################
template<class TV> FRAME<TV> RIGID_BODY_CLUSTER_BINDINGS<TV>::
Get_Child_To_Parent_Frame(int child_particle_index) const
{    
    assert(binding_index(child_particle_index)(0).y>=RIGID_CLUSTER_CONSTITUENT_ID());
    assert(Get_Parent_Index(child_particle_index)>=0);
    return reverse_bindings.Get(Get_Parent_Index(child_particle_index))->
        child_to_parent(binding_index(child_particle_index)(0).y);
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Positions
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clamp_Particles_To_Embedded_Positions() const
{
    for(T_REVERSE_BINDING_ITERATOR i(reverse_bindings);i.Valid();i.Next())
        Clamp_Particles_To_Embedded_Positions(i.Key());
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clamp_Particles_To_Embedded_Velocities() const
{
    for(T_REVERSE_BINDING_ITERATOR iterator(reverse_bindings);iterator.Valid();iterator.Next())
        Clamp_Particles_To_Embedded_Velocities(iterator.Key());
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clamp_Particles_To_Embedded_Velocities(ARRAY_VIEW<TWIST<TV> > twist) const
{
    for(T_REVERSE_BINDING_ITERATOR iterator(reverse_bindings);iterator.Valid();iterator.Next())
        Clamp_Particles_To_Embedded_Velocities(iterator.Key(),twist);
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TWIST<TV> > wrench_full) const
{
    for(T_REVERSE_BINDING_ITERATOR iterator(reverse_bindings);iterator.Valid();iterator.Next()){
        int parent=iterator.Key();
        const CLUSTER& cluster=*iterator.Data();
        if(cluster.active && !cluster.infinite_body) for(RIGID_CLUSTER_CONSTITUENT_ID j(0);j<cluster.children.Size();j++){
            const int child=cluster.children(j);
            if(rigid_body_collection.Rigid_Body(child).Has_Infinite_Inertia()) continue;
            wrench_full(parent).linear+=wrench_full(child).linear;
            wrench_full(parent).angular+=wrench_full(child).angular
                +TV::Cross_Product(rigid_body_collection.rigid_body_particles.frame(child).t-rigid_body_collection.rigid_body_particles.frame(parent).t,wrench_full(child).linear);}}
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Positions
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clamp_Particles_To_Embedded_Positions(const int parent) const
{
    const CLUSTER& cluster=*reverse_bindings.Get(parent);
    if(cluster.active) for(RIGID_CLUSTER_CONSTITUENT_ID j(0);j<cluster.children.Size();j++){
        const int child=cluster.children(j);
        if(rigid_body_collection.Rigid_Body(child).Has_Infinite_Inertia()) continue;
        rigid_body_collection.Rigid_Body(child).Frame()=rigid_body_collection.Rigid_Body(parent).Frame()*cluster.child_to_parent(j);}
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clamp_Particles_To_Embedded_Velocities(const int parent) const
{
    const CLUSTER& cluster=*reverse_bindings.Get(parent);
    if(cluster.active) for(RIGID_CLUSTER_CONSTITUENT_ID j(0);j<cluster.children.Size();j++){
        const int child=cluster.children(j);
        if(rigid_body_collection.Rigid_Body(child).Has_Infinite_Inertia()) continue;
        rigid_body_collection.rigid_body_particles.twist(child).linear=rigid_body_collection.Rigid_Body(parent).Pointwise_Object_Velocity(rigid_body_collection.rigid_body_particles.frame(child).t);
        rigid_body_collection.rigid_body_particles.twist(child).angular=rigid_body_collection.rigid_body_particles.twist(parent).angular;
        rigid_body_collection.Rigid_Body(child).Update_Angular_Momentum();}
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clamp_Particles_To_Embedded_Velocities(const int parent,ARRAY_VIEW<TWIST<TV> >& twist) const
{
    const CLUSTER& cluster=*reverse_bindings.Get(parent);
    if(cluster.active) for(RIGID_CLUSTER_CONSTITUENT_ID j(0);j<cluster.children.Size();j++){
        const int child=cluster.children(j);
        if(rigid_body_collection.Rigid_Body(child).Has_Infinite_Inertia()) continue;
        twist(child).linear=rigid_body_collection.Rigid_Body(parent).Pointwise_Object_Velocity(
            twist(parent),rigid_body_collection.rigid_body_particles.frame(parent).t,rigid_body_collection.rigid_body_particles.frame(child).t);
        twist(child).angular=twist(parent).angular;}
}
//#####################################################################
// Function Clear_Hard_Bound_Particles
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clear_Hard_Bound_Particles(ARRAY<bool>& particle_is_simulated) const
{
    for(T_REVERSE_BINDING_ITERATOR i(reverse_bindings);i.Valid();i.Next()){
        const CLUSTER& bindings=*i.Data();
        if(bindings.active) for(RIGID_CLUSTER_CONSTITUENT_ID j(0);j<bindings.children.Size();j++) particle_is_simulated(bindings.children(j))=false;}
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Distribute_Mass_To_All_Active_Parents()
{
    for(T_REVERSE_BINDING_ITERATOR i(reverse_bindings);i.Valid();i.Next())
        if(i.Data()->active) Clamp_Particles_To_Embedded_Velocities(i.Key());
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
namespace{
    void Distribute_Mass_To_Parents_Helper(RIGID_BODY_COLLECTION<VECTOR<float,1> >& rigid_body_collection,RIGID_BODY_CLUSTER_BINDINGS<VECTOR<float,1> >::CLUSTER& cluster,
        const int parent){}
    void Distribute_Mass_To_Parents_Helper(RIGID_BODY_COLLECTION<VECTOR<double,1> >& rigid_body_collection,RIGID_BODY_CLUSTER_BINDINGS<VECTOR<double,1> >::CLUSTER& cluster,
        const int parent){}
    template<class TV> void Distribute_Mass_To_Parents_Helper(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER& cluster,
        const int parent)
    {
        typedef typename TV::SCALAR T;
        typedef typename TV::SPIN T_SPIN;
        RIGID_BODY<TV>& parent_body=rigid_body_collection.Rigid_Body(parent);
        // compute aggregate inertia tensor and angular momentum
        SYMMETRIC_MATRIX<T,TV::SPIN::m> inertia_tensor;
        parent_body.Angular_Momentum()=T_SPIN();
        for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++){
            int child=cluster.children(i);RIGID_BODY<TV>& child_body=rigid_body_collection.Rigid_Body(child);
            TV s=(child_body.Frame().t-parent_body.Frame().t);
            inertia_tensor+=child_body.World_Space_Inertia_Tensor()+MATRIX_POLICY<TV>::CROSS_PRODUCT_MATRIX::Cross_Product_Matrix(child_body.Mass()*s).Times_Cross_Product_Matrix_Transpose_With_Symmetric_Result(s);
            parent_body.Angular_Momentum()+=TV::Cross_Product(child_body.Frame().t-parent_body.Frame().t,child_body.Mass()*child_body.Twist().linear)
                +child_body.Angular_Momentum();}
        // diagonalize inertia tensor and update angular velocity
        parent_body.Diagonalize_Inertia_Tensor(inertia_tensor);
    }
}
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Distribute_Mass_To_Parent(const int parent)
{
    RIGID_BODY_PARTICLES<TV>& rbp=rigid_body_collection.rigid_body_particles;
    CLUSTER& cluster=*reverse_bindings.Get(parent);
    if(cluster.infinite_body){
        RIGID_BODY<TV> *parent_body=&rigid_body_collection.Rigid_Body(parent),*child_body=&rigid_body_collection.Rigid_Body(cluster.infinite_body);
        parent_body->is_static=child_body->is_static;rbp.kinematic(parent)=rigid_body_collection.rigid_body_particles.kinematic(cluster.infinite_body);
        parent_body->Frame()=child_body->Frame();parent_body->Twist()=child_body->Twist();
        parent_body->Mass()=child_body->Mass();parent_body->Angular_Momentum()=child_body->Angular_Momentum();}
    else{
        rbp.frame(parent)=FRAME<TV>();
        rbp.twist(parent)=TWIST<TV>();
        rbp.mass(parent)=T();
        rbp.inertia_tensor(parent)=DIAGONAL_MATRIX<T,TV::SPIN::m>();
        // Find center of mass
        for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++){
            int child=cluster.children(i);
            rbp.mass(parent)+=rbp.mass(child);
            rbp.frame(parent).t+=rbp.frame(child).t*rbp.mass(child);
            rbp.twist(parent).linear+=rbp.twist(child).linear*rbp.mass(child);}
        rbp.frame(parent).t/=rbp.mass(parent);
        rbp.twist(parent).linear/=rbp.mass(parent);
        Distribute_Mass_To_Parents_Helper(rigid_body_collection,cluster,parent);}
    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++){
        int child=cluster.children(i);RIGID_BODY<TV>& child_body=rigid_body_collection.Rigid_Body(child);
        cluster.child_to_parent(i)=rigid_body_collection.Rigid_Body(parent).Frame().Inverse()*child_body.Frame();}
}
//#####################################################################
// Function Build_Aggregate_Geometry
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Build_Aggregate_Geometry(const int parent)
{
    CLUSTER& cluster=*reverse_bindings.Get(parent);
    ARRAY<T_SIMPLICIAL_OBJECT*> objects;
    ARRAY<FRAME<TV> > relative_frames;

    for(int j=0;j<rigid_body_collection.rigid_body_particles.structure_ids(parent).m;j++) if(rigid_body_collection.rigid_body_particles.structure_ids(parent)(j)>=0)
        rigid_body_collection.structure_list.Remove_Element(rigid_body_collection.rigid_body_particles.structure_ids(parent)(j));

    ARRAY<IMPLICIT_OBJECT<TV>*>* implicits=new ARRAY<IMPLICIT_OBJECT<TV>*>;
    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++){
        int child=cluster.children(i);RIGID_BODY<TV>& child_body=rigid_body_collection.Rigid_Body(child);
        if(child_body.simplicial_object){
            objects.Append(child_body.simplicial_object);
            relative_frames.Append(cluster.child_to_parent(i));}
        if(child_body.implicit_object)
            implicits->Append(new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(child_body.implicit_object->object_space_implicit_object,false,&cluster.child_to_parent(i)));}
    assert(implicits->m>0 || objects.m>0);
    if(objects.m>0) rigid_body_collection.Rigid_Body(parent).Add_Structure(*T_SIMPLICIAL_OBJECT::Union_Mesh_Objects_Relatively(objects,relative_frames));
    MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* multibody_levelset=new MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>(implicits);
    multibody_levelset->need_destroy_data=true;
    rigid_body_collection.Rigid_Body(parent).Add_Structure(*multibody_levelset);
    if(objects.m>0) rigid_body_collection.Rigid_Body(parent).simplicial_object->Update_Bounding_Box();
    rigid_body_collection.Rigid_Body(parent).Update_Bounding_Box();
    rigid_body_collection.Destroy_Unreferenced_Geometry();
}
//#####################################################################
// Function Update_Aggregate_Geometry
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Update_Aggregate_Geometry(const int parent)
{
    CLUSTER& cluster=*reverse_bindings.Get(parent);
    ARRAY<T_SIMPLICIAL_OBJECT*> objects;
    ARRAY<FRAME<TV> > relative_frames;

    for(RIGID_CLUSTER_CONSTITUENT_ID i(0);i<cluster.children.Size();i++){
        int child=cluster.children(i);RIGID_BODY<TV>& child_body=rigid_body_collection.Rigid_Body(child);
        if(child_body.simplicial_object){
            objects.Append(child_body.simplicial_object);
            relative_frames.Append(cluster.child_to_parent(i));}}
    if (objects.m>0){
        T_SIMPLICIAL_OBJECT::Union_Mesh_Objects_Relatively(rigid_body_collection.Rigid_Body(parent).simplicial_object,objects,relative_frames);
        rigid_body_collection.Rigid_Body(parent).simplicial_object->Clean_Memory();
        rigid_body_collection.Rigid_Body(parent).simplicial_object->Update_Bounding_Box();
        rigid_body_collection.Rigid_Body(parent).Update_Bounding_Box();}
}
//#####################################################################
// Function Binding
//#####################################################################
template<class TV> int RIGID_BODY_CLUSTER_BINDINGS<TV>::
Binding(const int child_particle,FRAME<TV>& frame) const
{
    if(child_particle>binding_index.m || binding_index(child_particle).m==0) return -1;
    int parent=binding_index(child_particle)(0).x;
    const RIGID_CLUSTER_CONSTITUENT_ID instance=binding_index(child_particle)(0).y;
    if(parent<0) return -1;
    if(!instance) return parent;
    const CLUSTER& binding=*reverse_bindings.Get(parent);
    assert(binding.children(instance)==child_particle);
    if(!binding.active) return -1;
    frame=binding.child_to_parent(instance);
    return parent;
}
//#####################################################################
// Function Clamp_Particle_To_Embedded_Position
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clamp_Particle_To_Embedded_Position(const int child) const
{
    if(rigid_body_collection.Rigid_Body(child).Has_Infinite_Inertia()) return;
    FRAME<TV> child_to_parent;
    int parent=Binding(child,child_to_parent);
    rigid_body_collection.Rigid_Body(child).Frame()=rigid_body_collection.Rigid_Body(parent).Frame()*child_to_parent;
}
//#####################################################################
// Function Clamp_Particle_To_Embedded_Velocity
//#####################################################################
template<class TV> void RIGID_BODY_CLUSTER_BINDINGS<TV>::
Clamp_Particle_To_Embedded_Velocity(const int child) const
{
    if(rigid_body_collection.Rigid_Body(child).Has_Infinite_Inertia()) return;
    FRAME<TV> child_to_parent;
    int parent=Binding(child,child_to_parent);
    rigid_body_collection.rigid_body_particles.twist(child).linear=rigid_body_collection.Rigid_Body(parent).Pointwise_Object_Velocity(rigid_body_collection.rigid_body_particles.frame(child).t);
    rigid_body_collection.rigid_body_particles.twist(child).angular=rigid_body_collection.rigid_body_particles.twist(parent).angular;
    rigid_body_collection.Rigid_Body(child).Update_Angular_Momentum();
}
//#####################################################################
// Function Valid_Cluster_Collision_Helper
//#####################################################################
template<class TV> bool RIGID_BODY_CLUSTER_BINDINGS<TV>::
Valid_Cluster_Collision_Helper(const int rigid_particle_index,const int parent)
{
    if(parent<0) return true;
    bool active_cluster=reverse_bindings.Get(parent)->active,is_parent=parent==rigid_particle_index,is_child=parent>0 && parent!=rigid_particle_index;
    assert(!active_cluster || is_parent || is_child);
    assert(!is_parent || !is_child);
    if(!active_cluster) return is_child;
    if(is_parent) return !collide_constituent_bodies;
    return collide_constituent_bodies;
}
//#####################################################################
// Function Clamp_Particle_To_Embedded_Position
//#####################################################################
template<class TV> bool RIGID_BODY_CLUSTER_BINDINGS<TV>::
Valid_Cluster_Collision(const int rigid_body_1,const int rigid_body_2)
{
    if(rigid_body_1<0 || rigid_body_2<0) return false;
    int parent_1=((rigid_body_1<binding_index.m) && binding_index(rigid_body_1).m>0)?binding_index(rigid_body_1)(0).x:-1,
        parent_2=((rigid_body_2<binding_index.m) && binding_index(rigid_body_2).m>0)?binding_index(rigid_body_2)(0).x:-1;
    if(parent_1<0 && parent_2<0) return true; // neither bound
    if(parent_1==parent_2) return false;
    bool decision=Valid_Cluster_Collision_Helper(rigid_body_1,parent_1)&&Valid_Cluster_Collision_Helper(rigid_body_2,parent_2);
    //PHYSBAM_DEBUG_PRINT("query",rigid_body_1,rigid_body_2,parent_1,parent_2,collide_constituent_bodies,decision);
    return decision;
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER::
CLUSTER()
    :active(false),stored_active(false),infinite_body(0),parent(0)
{}
//#####################################################################
namespace PhysBAM{
template class RIGID_BODY_CLUSTER_BINDINGS<VECTOR<float,1> >;
template class RIGID_BODY_CLUSTER_BINDINGS<VECTOR<float,2> >;
template class RIGID_BODY_CLUSTER_BINDINGS<VECTOR<float,3> >;
template class RIGID_BODY_CLUSTER_BINDINGS<VECTOR<double,1> >;
template class RIGID_BODY_CLUSTER_BINDINGS<VECTOR<double,2> >;
template class RIGID_BODY_CLUSTER_BINDINGS<VECTOR<double,3> >;
}

