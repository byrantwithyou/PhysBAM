//#####################################################################
// Copyright 2004-2008, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINTS_POOL
//#####################################################################
#ifndef __POINTS_POOL__
#define __POINTS_POOL__

#include <Core/Data_Structures/STACK.h>
#include <Core/Log/LOG.h>
#include <Tools/Clone/CLONE_ARRAY.h>
namespace PhysBAM{

template<class T_PARTICLES>
class POINTS_POOL:public NONCOPYABLE
{
private:
    typedef typename T_PARTICLES::VECTOR_T TV;

    const T_PARTICLES& template_particles;
    ARRAY<CLONE_ARRAY<T_PARTICLES>*> allocated_batches; // pointer to the first particle in the batch (of size allocation_batch_size)
    STACK<T_PARTICLES*> free_pool;
    int allocation_batch_size;
public:
    int number_particles_per_cell;

    POINTS_POOL(const T_PARTICLES& template_particles,const int allocation_batch_size_input=100000)
        :template_particles(template_particles),allocation_batch_size(allocation_batch_size_input),number_particles_per_cell(0)
    {
    }

    ~POINTS_POOL()
    {allocated_batches.Delete_Pointers_And_Clean_Memory();}

    void Set_Number_Particles_Per_Cell(const int number_particles_per_cell_input)
    {number_particles_per_cell=number_particles_per_cell_input;}

    T_PARTICLES* Allocate_Particle()
    {if(free_pool.Empty()){
        Allocate_New_Batch();
    }
    T_PARTICLES* particles=free_pool.Pop();
    return particles;}

    void Free_Particle(T_PARTICLES* particle_class)
    {if(particle_class){
        Free_Particle(particle_class->next);particle_class->next=0;particle_class->Delete_All_Elements();
        free_pool.Push(particle_class);}}

    T_PARTICLES& Add_Particle(T_PARTICLES& particles,int& index)
    {T_PARTICLES* particles_link=&particles;
    while(particles_link->Size()==number_particles_per_cell){ // find the right link (allocate if necessary)
        if(!particles_link->next) particles_link->next=Allocate_Particle();
        particles_link=particles_link->next;}
    index=particles_link->Add_Element();
    return *particles_link;}

private:
    void Allocate_New_Batch()
    {PHYSBAM_ASSERT(!template_particles.Size()); // make sure our clones are empty
    LOG::cout<<"Particle Pool: Allocating another "<<allocation_batch_size<<" particles"<<std::endl;
    CLONE_ARRAY<T_PARTICLES>* batch=new CLONE_ARRAY<T_PARTICLES>(template_particles,allocation_batch_size);
    for(int i=0;i<allocation_batch_size;i++) (*batch)(i).Preallocate(number_particles_per_cell);
    allocated_batches.Append(batch);
    free_pool.Increase_Size(allocation_batch_size);
    for(int i=0;i<allocation_batch_size;i++) free_pool.Push(&(*batch)(i));
    }

//#####################################################################
};
}
#endif
