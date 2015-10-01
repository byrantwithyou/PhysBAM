//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/KD_TREE.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/PINNING_FORCE.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Hybrid_Methods/Forces/MPM_OLDROYD_FINITE_ELEMENTS.h>
#include <Hybrid_Methods/Forces/MPM_VISCOSITY.h>
#include <Hybrid_Methods/Forces/OLDROYD_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Forces/VOLUME_PRESERVING_OB_NEO_HOOKEAN.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include "STANDARD_TESTS_KKT_2D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> STANDARD_TESTS_KKT<VECTOR<T,2> >::
STANDARD_TESTS_KKT(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :STANDARD_TESTS_KKT_BASE<TV>(stream_type_input,parse_args),
    use_surface_tension(false),Nsurface(0)
{
    parse_args.Parse();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> STANDARD_TESTS_KKT<VECTOR<T,2> >::
~STANDARD_TESTS_KKT()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,2> >::
Write_Output_Files(const int frame)
{
    BASE::Write_Output_Files(frame);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,2> >::
Read_Output_Files(const int frame)
{
    BASE::Read_Output_Files(frame);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,2> >::
Initialize()
{
    switch(test_number)
    {
        case 1:{ // circle free fall
            grid.Initialize(TV_INT()+resolution*2-1,RANGE<TV>::Unit_Box(),true);
            SPHERE<TV> sphere(TV(.5,.5),.3);
            T density=2*scale_mass;
            Seed_Particles_Helper(sphere,[=](const TV& X){return TV();},
                [=](const TV&){return MATRIX<T,2>();},
                density,particles_per_cell);
            T total_mass=particles.mass.Sum();
            TV total_momentum=particles.V.Weighted_Sum(particles.mass);
            TV dV=total_momentum/total_mass;
            particles.V-=dV;
            Add_Gravity(TV(0,-0.8));
            Add_Neo_Hookean(1e3*scale_E,0.3);
        } break;
        case 2:{ // full box
            grid.Initialize(TV_INT()+resolution*2-1,RANGE<TV>::Unit_Box(),true);
            RANGE<TV> box(TV(),TV::All_Ones_Vector());
            T density=2*scale_mass;
            Seed_Particles_Helper(box,0,0,density,particles_per_cell);
            Add_Gravity(TV(0,-1.8));
        } break;
        default: PHYSBAM_FATAL_ERROR("test number not implemented");
    }
    // initialize coarse grid
    coarse_grid.Initialize(TV_INT()+resolution,grid.domain.Thickened(grid.dX/2),true);
}
//#####################################################################
// Function Begin_Frame
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,2> >::
Begin_Frame(const int frame)
{
}
//#####################################################################
// Function End_Frame
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,2> >::
End_Frame(const int frame)
{
}
//#####################################################################
// Function Begin_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,2> >::
Begin_Time_Step(const T time)
{
    if(test_number==12){
        if(time>=10/24.0){
            lagrangian_forces.Delete_Pointers_And_Clean_Memory();
            this->deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();
            this->output_structures_each_frame=true;
            Add_Walls(-1,COLLISION_TYPE::separate,1.9,.1+(T)(time-10/24.0)*0.08,true);
            Add_Gravity(TV(0,-9.8));}}
    
    if(use_surface_tension){

        bool use_bruteforce=true;
        bool use_kdtree=false;

        // Remove old surface particles
        int N_non_surface=particles.number-Nsurface;
        for(int k=N_non_surface;k<particles.number;k++){
            particles.Add_To_Deletion_List(k);
            int m=steal(k-N_non_surface);
            TV old_momentum=particles.mass(m)*particles.V(m)+particles.mass(k)*particles.V(k);
            particles.mass(m)+=particles.mass(k);
            particles.volume(m)+=particles.volume(k);
            particles.V(m)=old_momentum/particles.mass(m);
            if(use_affine) particles.B(m)=(particles.B(m)+particles.B(k))*0.5;}
        LOG::cout<<"deleting "<<Nsurface<<" particles..."<<std::endl;
        particles.Delete_Elements_On_Deletion_List();
        lagrangian_forces.Delete_Pointers_And_Clean_Memory();
        this->deformable_body_collection.structures.Delete_Pointers_And_Clean_Memory();

        // Dirty hack of rasterizing mass
        this->simulated_particles.Remove_All();
        for(int p=0;p<this->particles.number;p++)
            if(this->particles.valid(p))
                this->simulated_particles.Append(p);
        this->particle_is_simulated.Remove_All();
        this->particle_is_simulated.Resize(this->particles.X.m);
        this->particle_is_simulated.Subset(this->simulated_particles).Fill(true);
        this->weights->Update(this->particles.X);
        this->gather_scatter.Prepare_Scatter(this->particles);
        MPM_PARTICLES<TV>& my_particles=this->particles;
#pragma omp parallel for
        for(int i=0;i<this->mass.array.m;i++)
            this->mass.array(i)=0;
        this->gather_scatter.template Scatter<int>(
            [this,&my_particles](int p,const PARTICLE_GRID_ITERATOR<TV>& it,int data)
            {
                T w=it.Weight();
                TV_INT index=it.Index();
                this->mass(index)+=w*my_particles.mass(p);
            },false);

        // Marching cube
        SEGMENTED_CURVE_2D<T>* surface=SEGMENTED_CURVE_2D<T>::Create();
        MARCHING_CUBES<TV>::Create_Surface(*surface,grid,mass,particles.mass(0)*0.4);

        // Seed surface particles
        int Nold=particles.number;
        Nsurface=surface->particles.number;
        LOG::cout<<"adding "<<Nsurface<<" particles..."<<std::endl;
        SEGMENTED_CURVE_2D<T>& new_sc=Seed_Lagrangian_Particles(*surface,0,0,(T)0.001,true);
        ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
        for(int i=Nold;i<particles.X.m;i++) (*color_attribute)(i)=VECTOR<T,3>(0,1,0);

        // Build K-d tree for non-surface particles
        LOG::cout<<"building kdtree..."<<std::endl;
        KD_TREE<TV> kdtree;
        ARRAY<TV> nodes(Nold);
        if(use_kdtree){
            for(int p=0;p<Nold;p++) nodes(p)=particles.X(p);
            kdtree.Create_Left_Balanced_KD_Tree(nodes);}

        // Assign physical quantities
        steal.Clean_Memory();
        for(int k=Nold;k<particles.number;k++){
            T dist2=FLT_MAX;
            int m=-1;

            // Find closest interior particle using brute force
            if(use_bruteforce){
                for(int q=0;q<Nold;q++){
                    T dd=(particles.X(q)-particles.X(k)).Magnitude_Squared();
                    if(dd<dist2){dist2=dd;m=q;}}}

            // Find closest interior particle using kdtree
            int number_of_points_in_estimate=1;
            ARRAY<int> points_found(number_of_points_in_estimate);
            ARRAY<T> squared_distance_to_points_found(number_of_points_in_estimate);
            if(use_kdtree){
                int number_of_points_found;T max_squared_distance_to_points_found;
                kdtree.Locate_Nearest_Neighbors(particles.X(k),FLT_MAX,points_found,
                    squared_distance_to_points_found,number_of_points_found,max_squared_distance_to_points_found,nodes);
                PHYSBAM_ASSERT(number_of_points_found==number_of_points_in_estimate);}

            // Debug kdtree
            if(use_bruteforce && use_kdtree && m!=points_found(0)){
                LOG::cout<<"Disagree!"<<std::endl; PHYSBAM_FATAL_ERROR();}

            if(use_kdtree) m=points_found(0);

            T split_mass=particles.mass(m)*.5;
            T split_volume=particles.volume(m)*.5;
            TV com=particles.X(m);
            particles.mass(m)=split_mass;
            particles.volume(m)=split_volume;
            particles.mass(k)=split_mass;
            particles.volume(k)=split_volume;
            particles.V(k)=particles.V(m);
            if(use_affine) particles.B(k)=particles.B(m);
            steal.Append(m);}

        // Add surface tension force
        SURFACE_TENSION_FORCE<TV>* stf=new SURFACE_TENSION_FORCE<TV>(new_sc,(T)1e-2);
        Add_Force(*stf);
    }

}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class T> void STANDARD_TESTS_KKT<VECTOR<T,2> >::
End_Time_Step(const T time)
{
}
template class STANDARD_TESTS_KKT<VECTOR<float,2> >;
template class STANDARD_TESTS_KKT<VECTOR<double,2> >;
}
