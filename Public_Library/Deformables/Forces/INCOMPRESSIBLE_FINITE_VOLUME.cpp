//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/SPARSE_UNION_FIND.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Math_Tools/Robust_Arithmetic.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <Tools/Parallel_Computation/MPI_WORLD.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology/TETRAHEDRON_MESH.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_PARTICLE_STATE.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
namespace PhysBAM{
//#####################################################################
static const int self_collision_subcycles=4;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
INCOMPRESSIBLE_FINITE_VOLUME(STRAIN_MEASURE<TV,d>& strain_measure)
    :DEFORMABLES_FORCES<TV>(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(strain_measure.particles)),strain_measure(strain_measure),
    disable_projection(false),minimum_volume_recovery_time_scale(0),max_cg_iterations(20),mpi_solids(0),merge_at_boundary(false),use_neumann(false),
    use_self_moving_projection(true),use_rigid_clamp_projection(true),use_diagonal_preconditioner(false),repulsions(0)
{
    T_MESH& mesh=strain_measure.mesh;
    if(!strain_measure.mesh.boundary_mesh) strain_measure.mesh.Initialize_Boundary_Mesh();
    if(!mesh.node_on_boundary) mesh.Initialize_Node_On_Boundary();
    if(!mesh.boundary_nodes) mesh.Initialize_Boundary_Nodes();

    LOG::cout<<"self_collision_subcycles "<<self_collision_subcycles<<std::endl;
    T_BOUNDARY_MESH& boundary_mesh=*strain_measure.mesh.boundary_mesh;
    if(!boundary_mesh.incident_elements) boundary_mesh.Initialize_Incident_Elements();
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();
    boundary_to_element.Resize(boundary_mesh.elements.m);
    for(int t=0;t<mesh.elements.m;t++){VECTOR<int,d+1>& element=mesh.elements(t);
        if(VECTOR<bool,d+1>(mesh.node_on_boundary->Subset(element)).Number_True()>=element.m-1) // using Number_True directly on the subset hits a compiler bug in gcc 4.1.1
            for(int i=0;i<element.m;i++){
                int b=boundary_mesh.Simplex(element.Remove_Index(i));
                if(b>=0) boundary_to_element(b).Set(t,i);}}

    node_regions.Resize(particles.Size());
    for(int p=0;p<particles.Size();p++){ARRAY<int>& incident=(*mesh.incident_elements)(p);
        for(int j=0;j<incident.m;j++) node_regions(p).Append_Unique_Elements(strain_measure.mesh.elements(incident(j)));}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
~INCOMPRESSIBLE_FINITE_VOLUME()
{
    cg_vectors.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int d> INCOMPRESSIBLE_FINITE_VOLUME<TV,d>* INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Create(T_OBJECT& object,const bool verbose)
{
    STRAIN_MEASURE<TV,d>* strain_measure=new STRAIN_MEASURE<TV,d>(object);
    if(verbose) strain_measure->Print_Altitude_Statistics();
    return new INCOMPRESSIBLE_FINITE_VOLUME(*strain_measure);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    // TODO: Do not compute force_particles_of_fragment here and in base class.
    Update_Force_Elements(force_elements,strain_measure.mesh.elements,particle_is_simulated);
    Update_Force_Particles(force_dynamic_particles,strain_measure.mesh.elements.Flattened(),particle_is_simulated,true);
    Update_Force_Elements(force_boundary_elements,strain_measure.mesh.boundary_mesh->elements,particle_is_simulated);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    BASE::Update_Position_Based_State(time,is_position_update,update_hessian);
    if(MPI_WORLD::Initialized() && !mpi_solids) PHYSBAM_FATAL_ERROR();
    ARRAY<T> element_volumes(strain_measure.mesh.elements.m,no_init); // TODO: this is not efficient.

    Bs_per_node.Resize(strain_measure.mesh.elements.m,no_init);
    for(int t:force_elements){
        T_MATRIX Ds=strain_measure.Ds(particles.X,t);
        element_volumes(t)=(T)1/factorial(d)*Ds.Parallelepiped_Measure();
        Bs_per_node(t)=(T)1/factorial(d+1)*Ds.Cofactor_Matrix();}

    boundary_normals.Resize(strain_measure.mesh.boundary_mesh->elements.m,no_init);

    for(int b:force_boundary_elements){
        int interior,i;boundary_to_element(b).Get(interior,i);
        if(i>1) boundary_normals(b)=-Bs_per_node(interior).Column(i-1);
        else boundary_normals(b)=Bs_per_node(interior)*VECTOR<T,d>::All_Ones_Vector();}

    volumes_full.Resize(particles.Size(),no_init);
    volumes_full.Subset(strain_measure.mesh.elements.Flattened()).Fill((T)0);;
    for(int t:force_elements){
        volumes_full.Subset(strain_measure.mesh.elements(t))+=element_volumes(t);}

    total_volume=0;
    for(int p:force_dynamic_particles){
        volumes_full(p)*=(T)1/(d+1);
        total_volume+=volumes_full(p);}

    if(mpi_solids) total_volume=mpi_solids->Reduce_Add_Global(total_volume);

    if(!rest_volumes_full.m){
        total_rest_volume=total_volume;
        rest_volumes_full.Resize(particles.Size());
        for(int p:force_dynamic_particles){
            rest_volumes_full(p)=volumes_full(p);}}

        // TODO(jontg): Make it work bleh.
    //if(!mpi_solids) LOG::cout<<"boundary volume: "<<rest_volumes_full.Subset(*strain_measure.mesh.boundary_nodes).Sum()<<std::endl;

    Update_Preconditioner();
}
//#####################################################################
// Class POISSON_SYSTEM
//#####################################################################
namespace{
template<class TV,int d>
class POISSON_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
protected:
    const INCOMPRESSIBLE_FINITE_VOLUME<TV,d>& fvm;
    const ARRAY<int>& dynamic_particles;
public:
    typedef INDIRECT_ARRAY<ARRAY<T> > VECTOR_T;
    typedef KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<T> > > KRYLOV_VECTOR_T;

    POISSON_SYSTEM(const INCOMPRESSIBLE_FINITE_VOLUME<TV,d>& fvm)
        :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),
        fvm(fvm),dynamic_particles(fvm.force_dynamic_particles)
    {}

    void Multiply(const KRYLOV_VECTOR_BASE<T>& bp,KRYLOV_VECTOR_BASE<T>& bresult,bool transpose=false) const override
    {const KRYLOV_VECTOR_T& p=debug_cast<const KRYLOV_VECTOR_T&>(bp);KRYLOV_VECTOR_T& result=debug_cast<KRYLOV_VECTOR_T&>(bresult);
    const KRYLOV_VECTOR_T* use_p=&p;
    INDIRECT_ARRAY<const ARRAY<T> > diagonal_preconditioner(fvm.diagonal_preconditioner_full,dynamic_particles);
    if(fvm.use_diagonal_preconditioner){use_p=&result;result.v=p.v;result.v*=diagonal_preconditioner;}
    fvm.Gradient(use_p->v.array,fvm.gradient_full);
    fvm.gradient_full.Subset(dynamic_particles)*=fvm.particles.one_over_mass.Subset(dynamic_particles);
    fvm.Project_Vector_Field(fvm.gradient_full);
    fvm.Negative_Divergence(fvm.gradient_full,result.v.array);
    if(fvm.use_diagonal_preconditioner) result.v*=diagonal_preconditioner;}

    void Project(KRYLOV_VECTOR_BASE<T>& p) const override
    {}

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& p) const override
    {Project(p);}

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bx,const KRYLOV_VECTOR_BASE<T>& by) const override
    {const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx),&y=debug_cast<const KRYLOV_VECTOR_T&>(by);
    assert(x.v.Size()==dynamic_particles.Size());
    T inner_product=x.v.Dot(y.v);
    if(fvm.mpi_solids) inner_product=fvm.mpi_solids->Reduce_Add(inner_product);
    return inner_product;}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bx) const override
    {const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);
    assert(x.v.Size()==dynamic_particles.Size());
    T convergence_norm=x.v.Maximum_Magnitude();
    if(fvm.mpi_solids) convergence_norm=fvm.mpi_solids->Reduce_Max(convergence_norm);
    return convergence_norm;}

    T Magnitude(const KRYLOV_VECTOR_BASE<T>& x) const
    {return sqrt((T)Inner_Product(x,x));}
};
}
//#####################################################################
// Function Make_Incompressible
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Make_Incompressible(const T dt,const bool correct_volume)
{
    LOG::SCOPE scope("PROJECTION","projection, dt = %g",dt);
    if(correct_volume) LOG::cout<<"Correcting Volume"<<std::endl;
    else LOG::cout<<"Correcting Divergence"<<std::endl;

    POISSON_SYSTEM<TV,d> system(*this);
    T max_error=0;
    
    for(int p:force_dynamic_particles){
        if(rest_volumes_full(p)){T error=(volumes_full(p)-rest_volumes_full(p))/rest_volumes_full(p);max_error=max(max_error,abs(error));}}
    if(mpi_solids) max_error=mpi_solids->Reduce_Max(max_error);
    LOG::cout<<"max error = "<<max_error<<", total volume = "<<total_volume<<", error = "<<Robust_Divide(total_volume-total_rest_volume,total_rest_volume)<<std::endl;

    if(disable_projection) return;

    if(boundary_pressures.m){
        if(TV::m!=d) PHYSBAM_FATAL_ERROR();
        gradient_full.Resize(particles.Size(),no_init);
        gradient_full.Subset(force_dynamic_particles).Fill(TV());
        T_BOUNDARY_MESH& boundary_mesh=*strain_measure.mesh.boundary_mesh;
        for(int t:force_boundary_elements){
            VECTOR<int,d>& element=boundary_mesh.elements(t);
            INDIRECT_ARRAY<ARRAY<T>,VECTOR<int,d>&> boundary_pressures_subset=boundary_pressures.Subset(element);
            T p_sum=boundary_pressures_subset.Sum();
            for(int i=0;i<d;i++) gradient_full(element[i])-=(p_sum+boundary_pressures(element[i]))*boundary_normals(t);}
        for(int p:force_dynamic_particles){
            particles.V(p)+=particles.one_over_mass(p)*gradient_full(p);}}

    Negative_Divergence(particles.V,divergence_full);
    KRYLOV_VECTOR_T divergence(divergence_full.Subset(force_dynamic_particles));

    if(correct_volume){
        T one_over_dt=1/dt,maximum_volume_recovery_fraction=dt/max(minimum_volume_recovery_time_scale,(T)1e-10);
        for(int p:force_dynamic_particles){
            T volume_error=rest_volumes_full(p)-volumes_full(p);
            volume_error=sign(volume_error)*min(abs(volume_error),maximum_volume_recovery_fraction*rest_volumes_full(p));
            divergence_full(p)+=one_over_dt*volume_error;}}

    system.Project(divergence);
    LOG::cout<<"divergence magnitude = "<<system.Magnitude(divergence)<<std::endl;

    pressure_full.Resize(particles.Size(),no_init);
    KRYLOV_VECTOR_T pressure(pressure_full.Subset(force_dynamic_particles));
    pressure.v.Fill((T)0);

    {CONJUGATE_RESIDUAL<T> cr;
    INDIRECT_ARRAY<ARRAY<T> > diagonal_preconditioner(diagonal_preconditioner_full,force_dynamic_particles);
    if(use_diagonal_preconditioner) divergence.v*=diagonal_preconditioner;
    T tolerance=max((T).01*system.Convergence_Norm(divergence),(T)1e-10);
    bool converged=cr.Solve(system,pressure,divergence,cg_vectors,tolerance,0,max_cg_iterations);
    if(use_diagonal_preconditioner) pressure.v*=diagonal_preconditioner;
    LOG::Stat("divergence magnitude",sqrt(cr.residual_magnitude_squared));
    LOG::Stat("nullspace measure",cr.nullspace_measure);
    if(!converged) LOG::cout<<"CONJUGATE_RESIDUAL FAILED - GIVING UP"<<std::endl;}

    Gradient(pressure_full,gradient_full);
    gradient_full.Subset(force_dynamic_particles)*=particles.one_over_mass.Subset(force_dynamic_particles);
    Project_Vector_Field(gradient_full);
    for(int i=0;i<force_dynamic_particles.m;i++){int p=force_dynamic_particles(i);particles.V(p)-=gradient_full(p);}
}
//#####################################################################
// Function Test_System
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Test_System()
{
    if(disable_projection) return;
    RANDOM_NUMBERS<T> random;random.Set_Seed(1823);
    Update_Position_Based_State(0,true,true);

    pressure_full.Resize(particles.Size(),no_init);
    divergence_full.Resize(particles.Size(),no_init);
    const ARRAY<int> &fragment_dynamic_particles=force_dynamic_particles,
        &fragment_particles=force_dynamic_particles;
    INDIRECT_ARRAY<ARRAY<T> > volumes(volumes_full,fragment_dynamic_particles);
    KRYLOV_VECTOR_T pressure(pressure_full.Subset(fragment_dynamic_particles)),divergence(divergence_full.Subset(fragment_dynamic_particles));
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > X(particles.X.Subset(fragment_dynamic_particles)),V(particles.V.Subset(fragment_dynamic_particles));
    INDIRECT_ARRAY<ARRAY<TV> > gradient(gradient_full.Subset(fragment_dynamic_particles));
    ARRAY<T> old_volumes(volumes);
    particles.V.Subset(fragment_particles).Fill(TV());
    pressure_full.Subset(fragment_particles).Fill(T());

    T dt=(T)1e-3;
    POISSON_SYSTEM<TV,d> system(*this);

    for(int iteration=0;iteration<10;iteration++){
        for(int p=0;p<fragment_dynamic_particles.m;p++){
            random.Set_Seed(fragment_dynamic_particles(p)*123+20+23*iteration);
            T scale=random.Get_Uniform_Number(-(T)10,(T)10);TV direction=random.template Get_Direction<TV>();
            V(p)=exp(scale)*direction;}
        for(int i=0;i<fragment_dynamic_particles.m;i++){int p=fragment_dynamic_particles(i);
            random.Set_Seed(p*123+21+23*iteration);
            pressure_full(p)=exp(random.Get_Uniform_Number(-(T)10,(T)10));}
        system.Project(pressure);

        T convergence_norm=V.Maximum_Magnitude();
        if(mpi_solids) convergence_norm=mpi_solids->Reduce_Max(convergence_norm);
        LOG::cout<<"|V|_inf = "<<convergence_norm<<std::endl;
        V/=convergence_norm;

        T magnitude_squared=V.Dot(V);
        if(mpi_solids) magnitude_squared=mpi_solids->Reduce_Add(magnitude_squared);
        T magnitude=sqrt(magnitude_squared);
        LOG::cout<<"|V|_before = "<<magnitude<<std::endl;
        pressure.v/=system.Convergence_Norm(pressure);
        LOG::cout<<"|V| = "<<magnitude<<", |p| = "<<system.Magnitude(pressure)<<std::endl;

        Project_Vector_Field(particles.V);
        Negative_Divergence(particles.V,divergence_full);system.Project(divergence);
        Gradient(pressure_full,gradient_full);
        Project_Vector_Field(gradient_full);

        T gradient_dot_V=gradient.Dot(V),divergence_dot_pressure=(T)system.Inner_Product(divergence,pressure);
        if(mpi_solids) gradient_dot_V=mpi_solids->Reduce_Add(gradient_dot_V);
        LOG::cout<<"gradient_dot_V = "<<gradient_dot_V<<", divergence_dot_pressure = "<<divergence_dot_pressure
            <<", relative error: "<<abs(Robust_Divide(gradient_dot_V,divergence_dot_pressure)-1)<<std::endl;

        ARRAY<TV> X_old(X);
        X+=dt*V;
        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(particles.X);
        Update_Position_Based_State(0,true,true);

        ARRAY<T> dv_full(volumes_full.m);
        KRYLOV_VECTOR_T dv(dv_full.Subset(fragment_dynamic_particles));
        dv.v=volumes-old_volumes;
        divergence*=dt;

        LOG::cout<<"|dv| = "<<system.Magnitude(dv)<<std::endl;
        LOG::cout<<"dv_dot_divergence / |divergence|^2 = "<<system.Inner_Product(dv,divergence)/system.Inner_Product(divergence,divergence)<<std::endl;

        // restore old positions
        X=X_old;
        Update_Position_Based_State(0,true,true);}

    exit(1);
}
//#####################################################################
// Function Gradient
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Gradient(ARRAY_VIEW<const T> p,ARRAY<TV>& gradient) const
{
    ARRAY_VIEW<T> modifiable_p(const_cast<T*>(p.Get_Array_Pointer()),p.m);
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(modifiable_p);
    gradient.Resize(particles.Size(),no_init);
    gradient.Subset(force_dynamic_particles).Fill(TV());
    for(int t:force_elements){
        VECTOR<int,d+1>& element=strain_measure.mesh.elements(t);
        INDIRECT_ARRAY<ARRAY_VIEW<const T>,VECTOR<int,d+1>&> p_subset=p.Subset(element);
        strain_measure.Distribute_Force(gradient,element,Bs_per_node(t)*-p_subset.Sum());}
}
//#####################################################################
// Function Negative_Divergence
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Negative_Divergence(ARRAY_VIEW<const TV> V,ARRAY<T>& divergence) const
{
    ARRAY_VIEW<TV> modifiable_V(const_cast<TV*>(V.Get_Array_Pointer()),V.m);
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(modifiable_V);
    divergence.Resize(particles.Size(),no_init);
    INDIRECT_ARRAY<ARRAY<T>,ARRAY<int>&> divergence_subset=divergence.Subset(force_dynamic_particles);
    divergence_subset.Fill(T());
    for(int t:force_elements){
        VECTOR<int,d+1>& element=strain_measure.mesh.elements(t);
        T divergence_per_node=T_MATRIX::Inner_Product(Bs_per_node(t),strain_measure.Ds(V,t));
        divergence.Subset(element)-=divergence_per_node;}
}
//#####################################################################
// Function Diagonal_Elements
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Diagonal_Elements(ARRAY<T>& D) const
{
    ARRAY<TV> forces(particles.Size());
    ARRAY<ARRAY<int> >& incident_elements=*strain_measure.mesh.incident_elements;
    for(int p:force_dynamic_particles){
        D(p)=T();
        forces.Subset(node_regions(p)).Fill(TV());
        for(int j=0;j<incident_elements(p).m;j++){int t=incident_elements(p)(j);
            strain_measure.Distribute_Force(forces,strain_measure.mesh.elements(t),-Bs_per_node(t));}
        forces.Subset(node_regions(p))*=particles.one_over_mass.Subset(node_regions(p));
        for(int j=0;j<incident_elements(p).m;j++){int t=incident_elements(p)(j);
            D(p)-=T_MATRIX::Inner_Product(Bs_per_node(t),strain_measure.Ds(forces,t));}}
}
//#####################################################################
// Function Diagonal_Elements
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Update_Preconditioner()
{
    if(!use_diagonal_preconditioner) return;
    diagonal_preconditioner_full.Resize(particles.Size(),no_init);
    Diagonal_Elements(diagonal_preconditioner_full);
    for(int p:force_dynamic_particles){
        diagonal_preconditioner_full(p)=1/sqrt(diagonal_preconditioner_full(p));}
}
//#####################################################################
// Function Project_Vector_Field
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Project_Vector_Field(ARRAY_VIEW<TV> field) const
{
    if(!use_neumann) return;
    if(self_collision_subcycles==1) PHYSBAM_FATAL_ERROR();
    Project_All_Isolated_Clamping_Constraints(field,projection_data);
    Project_All_Clamping_Constraints(field,projection_data);
    if(repulsions && use_self_moving_projection && d==3 && (projection_data.point_face_pairs.m || projection_data.edge_edge_pairs.m))
        for(int i=0;i<self_collision_subcycles;i++){
            TRIANGLE_REPULSIONS<TV>::Project_All_Moving_Constraints(projection_data.point_face_precomputed,projection_data.edge_edge_precomputed,field);
            Project_All_Clamping_Constraints(field,projection_data);}
}
//#####################################################################
// Function Project_All_Clamping_Constraints
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Project_All_Clamping_Constraints(ARRAY_VIEW<TV> field,const PROJECTION_DATA& data) const
{
    for(int i=0;i<data.neumann_boundary_nodes.m;i++){int p=data.neumann_boundary_nodes(i);
        field(p)-=TV::Dot_Product(field(p),data.neumann_boundary_normals(p))*data.neumann_boundary_normals(p);}
    for(int i=0;i<data.fixed_nodes.m;i++) field(data.fixed_nodes(i))=TV(); // TODO: generalize constructs
}
//#####################################################################
// Function Project_All_Isolated_Clamping_Constraints
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Project_All_Isolated_Clamping_Constraints(ARRAY_VIEW<TV> field,const PROJECTION_DATA& data) const
{
    for(int i=0;i<data.neumann_boundary_nodes_isolated.m;i++){int p=data.neumann_boundary_nodes_isolated(i);
        field(p)-=TV::Dot_Product(field(p),data.neumann_boundary_normals(p))*data.neumann_boundary_normals(p);}
}
//#####################################################################
// Function Set_Neumann_Boundary_Conditions
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Set_Neumann_Boundary_Conditions(const ARRAY<COLLISION_PARTICLE_STATE<TV> >* particle_states,TRIANGLE_REPULSIONS<TV>* repulsions_input)
{
    if(repulsions_input) repulsions=repulsions_input;

    neumann_boundary_count.Resize(particles.Size());
    // TODO: Unnecessarily expensive
    neumann_boundary_count.Fill(0);
    // TODO: Assuming we are created new simplifies this a little.
    projection_data.point_face_pairs.Remove_All();
    projection_data.edge_edge_pairs.Remove_All();
    projection_data.neumann_boundary_nodes.Remove_All();
    // TODO: Unnecessarily expensive
    projection_data.neumann_boundary_normals.Resize(particles.Size());projection_data.neumann_boundary_normals.Fill(TV());

    if(particle_states && use_rigid_clamp_projection)
        for(int p:force_dynamic_particles){
            const COLLISION_PARTICLE_STATE<TV>& collision=(*particle_states)(p);
            if(collision.enforce){
                projection_data.neumann_boundary_normals(p)=collision.normal;
                if(!neumann_boundary_count(p))
                    projection_data.neumann_boundary_nodes.Append(p);
                neumann_boundary_count(p)=1;}}

    if(repulsions && use_self_moving_projection){
        repulsions->Set_Collision_Pairs(projection_data.point_face_precomputed,projection_data.edge_edge_precomputed,projection_data.point_face_pairs,projection_data.edge_edge_pairs,(T)2); // TODO: Fix me.
        for(int i=0;i<projection_data.point_face_pairs.m;i++) neumann_boundary_count.Subset(projection_data.point_face_pairs(i).nodes)+=1;
        for(int i=0;i<projection_data.edge_edge_pairs.m;i++) neumann_boundary_count.Subset(projection_data.edge_edge_pairs(i).nodes)+=1;}

    int j=0;
    for(int i=0;i<projection_data.neumann_boundary_nodes.m;i++){int p=projection_data.neumann_boundary_nodes(i);
        if(neumann_boundary_count(p)==1) projection_data.neumann_boundary_nodes_isolated.Append(p);
        else projection_data.neumann_boundary_nodes(j++)=p;}
    projection_data.neumann_boundary_nodes.m=j;

    LOG::cout<<"moving constraints: "<<projection_data.point_face_pairs.m<<std::endl;
}
//#####################################################################
// Function Max_Relative_Velocity_Error
//#####################################################################
template<class TV,int d> typename TV::SCALAR INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Max_Relative_Velocity_Error()
{
    Update_Position_Based_State(0,true,true);
    T max_error=0;
    
    for(int i=0;i<force_dynamic_particles.m;i++){int p=force_dynamic_particles(i);
        if(rest_volumes_full(p)){T error=(volumes_full(p)-rest_volumes_full(p))/rest_volumes_full(p);max_error=max(max_error,abs(error));}}
    if(mpi_solids) max_error=mpi_solids->Reduce_Max_Global(max_error);
    return max_error;
}
//#####################################################################
// Function Save_Volumes
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Save_Volumes()
{
    Update_Position_Based_State(0,true,true);
    saved_volumes_full=volumes_full;
}
//#####################################################################
// Function Check_Improvement
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Check_Improvement()
{
    Update_Position_Based_State(0,true,true);
    T ave_improve=0,max_improve=-FLT_MAX,min_improve=FLT_MAX,old_total_volume_accumulated=0,new_total_volume_accumulated=0;
    int better=0,worse=0,same=0;
    const INDIRECT_ARRAY<const ARRAY<T> > volumes(volumes_full.Subset(force_dynamic_particles));
    const INDIRECT_ARRAY<const ARRAY<T> > saved_volumes(saved_volumes_full.Subset(force_dynamic_particles));
    const INDIRECT_ARRAY<const ARRAY<T> > rest_volumes(rest_volumes_full.Subset(force_dynamic_particles));
    T new_total_volume=volumes.Sum();
    if(mpi_solids) new_total_volume=mpi_solids->Reduce_Add_Global(new_total_volume);
    new_total_volume_accumulated+=new_total_volume;
    T old_total_volume=saved_volumes.Sum();
    if(mpi_solids) old_total_volume=mpi_solids->Reduce_Add_Global(old_total_volume);
    old_total_volume_accumulated+=old_total_volume;
    if(volumes.Size()!=rest_volumes.Size() || volumes.Size()!=saved_volumes.Size()){LOG::cout<<"Volume arrays have different sizes!"<<std::endl;PHYSBAM_FATAL_ERROR();}
    for(int i=0;i<volumes.Size();i++) if(rest_volumes(i)){
        T old_v=abs((saved_volumes(i)-rest_volumes(i))/rest_volumes(i)),new_v=abs((volumes(i)-rest_volumes(i))/rest_volumes(i)),diff=new_v-old_v;
        ave_improve+=diff;
        max_improve=max(max_improve,diff);
        min_improve=min(min_improve,diff);
        if(abs(new_v-old_v)<1e-2) same++;
        else if(new_v<old_v) better++;
        else worse++;}

     if(mpi_solids) max_improve=mpi_solids->Reduce_Max_Global(max_improve);
     if(mpi_solids) min_improve=mpi_solids->Reduce_Min_Global(min_improve);
     if(mpi_solids) ave_improve=mpi_solids->Reduce_Add_Global(ave_improve);
     if(mpi_solids) better=mpi_solids->Reduce_Add_Global(better);
     if(mpi_solids) same=mpi_solids->Reduce_Add_Global(same);
     if(mpi_solids) worse=mpi_solids->Reduce_Add_Global(worse);
    ave_improve/=particles.Size();
    LOG::cout<<"INCOMP STAT total volume error: "<<abs(old_total_volume_accumulated-total_rest_volume)<<" => "<<abs(new_total_volume_accumulated-total_rest_volume)<<std::endl;
    LOG::cout<<"INCOMP STAT better: "<<better<<"    same: "<<same<<"    worse: "<<worse<<std::endl;
    LOG::cout<<"INCOMP STAT improvements -   min: "<<min_improve<<"    ave: "<<ave_improve<<"    max: "<<max_improve<<std::endl;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV,int d> void INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    strain_measure.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<float,2>,2>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<float,3>,2>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<float,3>,3>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<double,2>,2>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<double,3>,2>;
template class INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<double,3>,3>;
}
