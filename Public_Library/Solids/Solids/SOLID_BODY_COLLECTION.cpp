//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Levine, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_BODY_COLLECTION
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Slender_Rods/RIGID_SLENDER_ROD_COLLECTION.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_BODY_COLLECTION<TV>::
SOLID_BODY_COLLECTION(DEFORMABLE_PARTICLES<TV>* particles)
    :collision_body_list(*new COLLISION_BODY_COLLECTION<TV>),
    deformable_body_collection(*new DEFORMABLE_BODY_COLLECTION<TV>(particles,&collision_body_list)),
    rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(&collision_body_list)),
    rigid_slender_rod_collection(*new RIGID_SLENDER_ROD_COLLECTION<TV>),
    print_energy(false),simulate(true),iterations_used_diagnostic(0)
{
    Print_Diagnostics();
    Print_Residuals(false);
    Set_CFL_Number();
    Set_Implicit_Damping();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_BODY_COLLECTION<TV>::
~SOLID_BODY_COLLECTION()
{
    solids_forces.Delete_Pointers_And_Clean_Memory();
    delete &deformable_body_collection;
    delete &rigid_body_collection;
    delete &rigid_slender_rod_collection;
    delete &collision_body_list;
}
//#####################################################################
// Function Delete_Forces
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Delete_Forces()
{
    solids_forces.Clean_Memory();
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Update_Simulated_Particles()
{
    rigid_body_collection.Update_Simulated_Particles();
    deformable_body_collection.Update_Simulated_Particles();
    int particles_number=deformable_body_collection.particles.Size();
    int rigid_particles_number=rigid_body_collection.rigid_body_particles.Size();

    ARRAY<bool> particle_is_simulated(particles_number);
    particle_is_simulated.Subset(deformable_body_collection.simulated_particles).Fill(true);

    ARRAY<bool> rigid_particle_is_simulated(rigid_particles_number);
    rigid_particle_is_simulated.Subset(rigid_body_collection.simulated_rigid_body_particles).Fill(true);
    for(int i=0;i<solids_forces.m;i++) solids_forces(i)->Update_Mpi(particle_is_simulated,rigid_particle_is_simulated,deformable_body_collection.mpi_solids);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    for(int k=0;k<solids_forces.m;k++) solids_forces(k)->Update_Position_Based_State(time);
    for(int k=0;k<rigid_body_collection.rigids_forces.m;k++) rigid_body_collection.rigids_forces(k)->Update_Position_Based_State(time);
    for(int k=0;k<deformable_body_collection.deformables_forces.m;k++) deformable_body_collection.deformables_forces(k)->Update_Position_Based_State(time,is_position_update,update_hessian);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    assert(F.V.array.Size()==deformable_body_collection.particles.Size());
    for(int k=0;k<solids_forces.m;k++)
        solids_forces(k)->Add_Velocity_Independent_Forces(F,time);
    for(int k=0;k<rigid_body_collection.rigids_forces.m;k++)
        rigid_body_collection.rigids_forces(k)->Add_Velocity_Independent_Forces(F.rigid_V.array,time);
    for(int k=0;k<deformable_body_collection.deformables_forces.m;k++)
        deformable_body_collection.deformables_forces(k)->Add_Velocity_Independent_Forces(F.V.array,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
// can depend on position too
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time) const
{
    assert(F.V.array.Size()==deformable_body_collection.particles.Size());
    for(int k=0;k<solids_forces.m;k++)
        solids_forces(k)->Add_Velocity_Dependent_Forces(V,F,time);
    for(int k=0;k<rigid_body_collection.rigids_forces.m;k++)
        rigid_body_collection.rigids_forces(k)->Add_Velocity_Dependent_Forces(V.rigid_V.array,F.rigid_V.array,time);
    for(int k=0;k<deformable_body_collection.deformables_forces.m;k++)
        deformable_body_collection.deformables_forces(k)->Add_Velocity_Dependent_Forces(V.V.array,F.V.array,time);
}
//#####################################################################
// Function Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,
    GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose) const
{
    assert(V.V.array.Size()==deformable_body_collection.particles.Size() && F.V.array.Size()==deformable_body_collection.particles.Size());
    assert(F.rigid_V.array.Size()==rigid_body_collection.rigid_body_particles.Size());
    for(int k=0;k<solids_forces.m;k++)
        solids_forces(k)->Add_Implicit_Velocity_Independent_Forces(V,F,time,transpose);
    for(int k=0;k<rigid_body_collection.rigids_forces.m;k++)
        rigid_body_collection.rigids_forces(k)->Add_Implicit_Velocity_Independent_Forces(V.rigid_V.array,F.rigid_V.array,time,transpose);
    for(int k=0;k<deformable_body_collection.deformables_forces.m;k++)
        deformable_body_collection.deformables_forces(k)->Add_Implicit_Velocity_Independent_Forces(V.V.array,F.V.array,time,transpose);
}
//#####################################################################
// Function Force_Differential
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Force_Differential(const GENERALIZED_VELOCITY<TV>& dX,GENERALIZED_VELOCITY<TV>& dF,const T time) const
{
    assert(dX.V.array.Size()==deformable_body_collection.particles.Size() && dF.V.array.Size()==deformable_body_collection.particles.Size());
    dF.V.array.Subset(deformable_body_collection.simulated_particles).Fill(TV());
    for(int k=0;k<deformable_body_collection.deformables_forces.m;k++)
        deformable_body_collection.deformables_forces(k)->Add_Implicit_Velocity_Independent_Forces(dX.V.array,dF.V.array,time);
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
    for(int k=0;k<deformable_body_collection.deformables_forces.m;k++)
        deformable_body_collection.deformables_forces(k)->Enforce_Definiteness(enforce_definiteness_input);
}
//#####################################################################
// Function Update_CFL
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Update_CFL()
{
    bool cfl_valid=true;
    if(solids_forces.m || rigid_body_collection.rigids_forces.m || deformable_body_collection.deformables_forces.m){
        for(int i=0;i<solids_forces.m;i++){if(!solids_forces(i)->CFL_Valid()){cfl_valid=false;break;}}
        for(int i=0;i<rigid_body_collection.rigids_forces.m;i++){if(!rigid_body_collection.rigids_forces(i)->CFL_Valid()){cfl_valid=false;break;}}
        for(int i=0;i<deformable_body_collection.deformables_forces.m;i++){if(!deformable_body_collection.deformables_forces(i)->CFL_Valid()){cfl_valid=false;break;}}}
    else cfl_valid=false;
    if(!cfl_valid){
        frequency.Resize(deformable_body_collection.particles.Size(),no_init);
        frequency.Subset(deformable_body_collection.simulated_particles).Fill(T_FREQUENCY_DEFORMABLE());

        rigid_frequency.Resize(rigid_body_collection.rigid_body_particles.Size(),no_init);
        rigid_frequency.Subset(rigid_body_collection.simulated_rigid_body_particles).Fill(T_FREQUENCY_RIGID());

        for(int i=0;i<solids_forces.m;i++){solids_forces(i)->Initialize_CFL(frequency,rigid_frequency);solids_forces(i)->Validate_CFL();}
        for(int i=0;i<rigid_body_collection.rigids_forces.m;i++){rigid_body_collection.rigids_forces(i)->Initialize_CFL(rigid_frequency);rigid_body_collection.rigids_forces(i)->Validate_CFL();}
        for(int i=0;i<deformable_body_collection.deformables_forces.m;i++){deformable_body_collection.deformables_forces(i)->Initialize_CFL(frequency);deformable_body_collection.deformables_forces(i)->Validate_CFL();}
        cfl_elastic=FLT_MAX;cfl_damping=FLT_MAX;
        for(int i=0;i<deformable_body_collection.simulated_particles.m;i++){int p=deformable_body_collection.simulated_particles(i);
            cfl_elastic=min(cfl_elastic,Robust_Inverse(sqrt(frequency(p).elastic_squared)));
            cfl_damping=min(cfl_damping,Robust_Inverse(frequency(p).damping));}
        for(int i=0;i<rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=rigid_body_collection.simulated_rigid_body_particles(i);
            cfl_elastic=min(cfl_elastic,Robust_Inverse(sqrt(rigid_frequency(p).elastic_squared)));
            cfl_damping=min(cfl_damping,Robust_Inverse(rigid_frequency(p).damping));}}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR SOLID_BODY_COLLECTION<TV>::
CFL(const bool verbose)
{
    T dt_elastic_and_damping=CFL_Elastic_And_Damping(),dt_strain_rate=CFL_Strain_Rate();
    if(verbose){
        LOG::cout<<"dt_elastic_and_damping = "<<dt_elastic_and_damping<<std::endl;
        LOG::cout<<"dt_strain_rate = "<<dt_strain_rate<<std::endl;
        LOG::cout<<"min = "<<min(dt_elastic_and_damping,dt_strain_rate)<<std::endl;}
    return min(dt_elastic_and_damping,dt_strain_rate);
}
//#####################################################################
// Function CFL_Elastic_And_Damping
//#####################################################################
template<class TV> typename TV::SCALAR SOLID_BODY_COLLECTION<TV>::
CFL_Elastic_And_Damping()
{
    T dt_elastic=CFL_Elastic();
    T dt_damping=FLT_MAX;if(!implicit_damping) dt_damping=CFL_Damping();
    T one_over_dt_full=1/dt_elastic+1/dt_damping;
    return Robust_Divide((T)1,one_over_dt_full);
}
//#####################################################################
// Function CFL_Elastic
//#####################################################################
template<class TV> typename TV::SCALAR SOLID_BODY_COLLECTION<TV>::
CFL_Elastic()
{
    Update_CFL();
    return cfl_elastic;
}
//#####################################################################
// Function CFL_Damping
//#####################################################################
template<class TV> typename TV::SCALAR SOLID_BODY_COLLECTION<TV>::
CFL_Damping()
{
    Update_CFL();
    return cfl_damping;
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR SOLID_BODY_COLLECTION<TV>::
CFL_Strain_Rate()
{
    T dt_strain=FLT_MAX;
    for(int k=0;k<solids_forces.m;k++) if(solids_forces(k)->limit_time_step_by_strain_rate) dt_strain=min(dt_strain,solids_forces(k)->CFL_Strain_Rate()); // otherwise not included
    for(int k=0;k<deformable_body_collection.deformables_forces.m;k++) if(deformable_body_collection.deformables_forces(k)->limit_time_step_by_strain_rate) dt_strain=min(dt_strain,deformable_body_collection.deformables_forces(k)->CFL_Strain_Rate()); // otherwise not included
    return dt_strain;
}
//#####################################################################
// Function Compute_Linear_Momentum 
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Compute_Linear_Momentum(TV& linear_momentum) const
{
    linear_momentum=TV();
    for(int i=0;i<deformable_body_collection.dynamic_particles.m;i++){int p=deformable_body_collection.dynamic_particles(i);
        linear_momentum+=deformable_body_collection.particles.mass(p)*deformable_body_collection.particles.V(p);}
    for(int i=0;i<rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        linear_momentum+=rigid_body_collection.Rigid_Body(p).Mass()*rigid_body_collection.Rigid_Body(p).Twist().linear;}
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const
{
    potential_energy=0;
    kinetic_energy=0;
    for(int i=0;i<solids_forces.m;i++) potential_energy+=solids_forces(i)->Potential_Energy(time);
    for(int i=0;i<rigid_body_collection.rigids_forces.m;i++) potential_energy+=rigid_body_collection.rigids_forces(i)->Potential_Energy(time);
    for(int i=0;i<deformable_body_collection.deformables_forces.m;i++) potential_energy+=deformable_body_collection.deformables_forces(i)->Potential_Energy(time);
    for(int i=0;i<deformable_body_collection.dynamic_particles.m;i++){int p=deformable_body_collection.dynamic_particles(i);
        kinetic_energy+=(T).5*deformable_body_collection.particles.mass(p)*TV::Dot_Product(deformable_body_collection.particles.V(p),deformable_body_collection.particles.V(p));}
    for(int i=0;i<rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        kinetic_energy+=rigid_body_collection.Rigid_Body(p).Kinetic_Energy();}
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> TV SOLID_BODY_COLLECTION<TV>::
Compute_Momentum() const
{
    TV momentum=deformable_body_collection.particles.V.Subset(deformable_body_collection.dynamic_particles).
        Weighted_Sum(deformable_body_collection.particles.mass.Subset(deformable_body_collection.dynamic_particles));
    for(int i=0;i<rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=rigid_body_collection.dynamic_rigid_body_particles(i);
        momentum+=rigid_body_collection.Rigid_Body(p).Momentum();}
    return momentum;
}
//#####################################################################
// Function Print_Energy
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Print_Energy(const T time,const int step) const
{
    if(print_energy){
        T potential_energy=0,kinetic_energy=0;
        Compute_Energy(time,kinetic_energy,potential_energy);
        LOG::cout<<"total energy = "<<(potential_energy+kinetic_energy)<<"    (KE = "<<kinetic_energy<<"   PE = "<<potential_energy<<")  Step "<<step<<std::endl;
        LOG::cout<<"total momentum = "<<Compute_Momentum()<<"  Step "<<step<<std::endl;}
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Read(const VIEWER_DIR& viewer_dir,const bool static_variables_every_frame,
    ARRAY<int>* needs_init,ARRAY<int>* needs_destroy)
{
    deformable_body_collection.Read(viewer_dir,static_variables_every_frame);
    rigid_body_collection.Read(viewer_dir,needs_init,needs_destroy);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir,const bool static_variables_every_frame,const bool write_from_every_process) const
{
    deformable_body_collection.Write(stream_type,viewer_dir,static_variables_every_frame,write_from_every_process);
    ARRAY<FORCE_DATA<TV> > spring_data_list;
    for(int i=0;i<solids_forces.m;i++) solids_forces(i)->Add_Force_Data(spring_data_list);
    for(int i=0;i<deformable_body_collection.deformables_forces.m;i++) deformable_body_collection.deformables_forces(i)->Add_Force_Data(spring_data_list);
    if(spring_data_list.m!=0) Write_To_File(stream_type,viewer_dir.current_directory+"/force_data",spring_data_list);
    rigid_body_collection.Write(stream_type,viewer_dir);
    //deformable_body_collection.triangle_repulsions.Output_Interaction_Pairs(stream_type,viewer_dir.current_directory+"/interaction_pairs");
}
//#####################################################################
// Function Add_All_Forces
//#####################################################################
template<class TV> void SOLID_BODY_COLLECTION<TV>::
Add_All_Forces(GENERALIZED_VELOCITY<TV>& F,const T time,const bool damping_only)
{
    if(!damping_only) Add_Velocity_Independent_Forces(F,time);
    GENERALIZED_VELOCITY<TV> V(*this);
    Add_Velocity_Dependent_Forces(V,F,time);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int SOLID_BODY_COLLECTION<TV>::
Add_Force(SOLIDS_FORCES<TV>* force)
{
    solids_forces.Append(force);
    force->Set_CFL_Number(cfl_number);
    return solids_forces.m-1;
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int SOLID_BODY_COLLECTION<TV>::
Add_Force(DEFORMABLES_FORCES<TV>* force)
{
    deformable_body_collection.deformables_forces.Append(force);
    force->Set_CFL_Number(cfl_number);
    return deformable_body_collection.deformables_forces.m-1;
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int SOLID_BODY_COLLECTION<TV>::
Add_Force(RIGIDS_FORCES<TV>* force)
{
    rigid_body_collection.rigids_forces.Append(force);
    force->Set_CFL_Number(cfl_number);
    return rigid_body_collection.rigids_forces.m-1;
}
//#####################################################################
// Function New_Generalized_Velocity
//#####################################################################
template<class TV> GENERALIZED_VELOCITY<TV>& SOLID_BODY_COLLECTION<TV>::
New_Generalized_Velocity() const
{
    GENERALIZED_VELOCITY<TV> gv(*this);
    return static_cast<GENERALIZED_VELOCITY<TV>&>(*gv.Clone_Default());
}
//#####################################################################
namespace PhysBAM{
template class SOLID_BODY_COLLECTION<VECTOR<float,1> >;
template class SOLID_BODY_COLLECTION<VECTOR<float,2> >;
template class SOLID_BODY_COLLECTION<VECTOR<float,3> >;
template class SOLID_BODY_COLLECTION<VECTOR<double,1> >;
template class SOLID_BODY_COLLECTION<VECTOR<double,2> >;
template class SOLID_BODY_COLLECTION<VECTOR<double,3> >;
}

