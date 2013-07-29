//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Frank Losasso, Duc Nguyen, Nick Rasmussen, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_UNIFORM
//#####################################################################
#ifndef __PARTICLE_LEVELSET_UNIFORM__
#define __PARTICLE_LEVELSET_UNIFORM__

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Incompressible/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD.h>
#ifdef USE_PTHREADS
#include <Tools/Parallel_Computation/PTHREAD.h>
#endif
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET.h>
namespace PhysBAM{

template<class TV> class FLUID_MASS_CONSERVATION_MODIFY_PHI;
template<class T> class MPI_UNIFORM_GRID;

template<class TV>
class PARTICLE_LEVELSET_UNIFORM:public PARTICLE_LEVELSET<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename GRID<TV>::BLOCK T_BLOCK;
    typedef typename ARRAY<T,TV_INT>::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename ARRAY<T,TV_INT>::template REBIND<char>::TYPE T_ARRAYS_CHAR;
    typedef typename ARRAY<T,TV_INT>::template REBIND<ARRAY<bool> >::TYPE T_ARRAYS_ARRAY_BOOL;
    typedef typename ARRAY<T,TV_INT>::template REBIND<PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename ARRAY<T,TV_INT>::template REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
    typedef typename ARRAY<T,TV_INT>::template REBIND<ARRAY<TV>*>::TYPE T_ARRAYS_ARRAY_TV;
    typedef LINEAR_INTERPOLATION_MAC_HELPER<TV> T_LINEAR_INTERPOLATION_MAC_HELPER;
    typedef typename LINEAR_INTERPOLATION_UNIFORM<TV,T>::template REBIND<TV>::TYPE T_LINEAR_INTERPOLATION_VECTOR;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> T_FACE_LOOKUP_COLLIDABLE;
public:
    typedef PARTICLE_LEVELSET<TV> BASE;
    using BASE::half_band_width;using BASE::use_removed_negative_particles;using BASE::use_removed_positive_particles;using BASE::store_unique_particle_id;using BASE::last_unique_particle_id;
    using BASE::number_particles_per_cell;using BASE::maximum_particle_radius;using BASE::minimum_particle_radius;using BASE::random;using BASE::outside_particle_distance_multiplier;
    using BASE::maximum_iterations_for_attraction;using BASE::bias_towards_negative_particles;using BASE::Set_Number_Particles_Per_Cell;
    using BASE::min_collision_distance_factor;using BASE::max_minus_min_collision_distance_factor_over_max_short;
    using BASE::reincorporate_removed_particles_everywhere;using BASE::save_removed_particle_times;using BASE::removed_particle_times;using BASE::only_use_negative_particles;
    using BASE::velocity_interpolation_collidable_contour_value;using BASE::template_particles;using BASE::template_removed_particles;
    using BASE::levelset;using BASE::negative_particles;using BASE::positive_particles;using BASE::removed_negative_particles;using BASE::removed_positive_particles;
    using BASE::escaped_negative_particles;using BASE::escaped_positive_particles;using BASE::Set_Band_Width;using BASE::cfl_number;using BASE::number_of_ghost_cells;using BASE::particle_pool;
    using BASE::Adjust_Particle_For_Objects;using BASE::Delete_Particle_And_Clean_Memory;using BASE::Particle_Collision_Distance;using BASE::Allocate_Particles;
    using BASE::Free_Particle_And_Clear_Pointer;using BASE::Add_Particle;using BASE::Delete_Particles_From_Deletion_List;using BASE::Move_Particle;using BASE::Copy_Particle;
    using BASE::Get_Particle_Link;using BASE::Add_Particle_To_Deletion_List;using BASE::deletion_list;

    MPI_UNIFORM_GRID<TV>* mpi_grid;

    T_ARRAYS_ARRAY_TV positive_particle_positions,negative_particle_positions;

    THREAD_QUEUE* thread_queue;
#ifdef USE_PTHREADS
    pthread_mutex_t cell_lock;
    pthread_barrier_t cell_barr;
#endif

    PARTICLE_LEVELSET_UNIFORM(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_ghost_cells_input);
    ~PARTICLE_LEVELSET_UNIFORM();

    void Set_Thread_Queue(THREAD_QUEUE* thread_queue_input)
    {thread_queue=thread_queue_input;}

//#####################################################################
    void Seed_Particles(const T time,const bool verbose=true);
    void Adjust_Particle_Radii();
    void Adjust_Particle_Radii_Threaded(RANGE<TV_INT>& domain);
    void Modify_Levelset_Using_Escaped_Particles(ARRAY<T,FACE_INDEX<TV::m> >* V,ARRAY<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES*>* other_positive_particles=0);
    void Update_Particles_To_Reflect_Mass_Conservation(ARRAY<T,TV_INT>& phi_old,const bool update_particle_cells=true,const bool verbose=false);
    void Update_Particles_To_Reflect_Mass_Conservation(const LEVELSET<TV>& levelset_old,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const bool update_particle_cells,const bool verbose);
    void Euler_Step_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,const T dt,const T time,const bool use_second_order_for_nonremoved_particles=false,
        const bool update_particle_cells_after_euler_step=true,const bool verbose=true,const bool analytic_test=false);
    void Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    void Euler_Step_Removed_Particles(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    template<class T_ARRAYS_PARTICLES> void Euler_Step_Particles_Wrapper(const ARRAY<T,FACE_INDEX<TV::m> >& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool update_particle_cells_after_euler_step=true,const bool assume_particles_in_correct_blocks=true,const bool enforce_domain_boundaries=true);
    template<class T_ARRAYS_PARTICLES> void Euler_Step_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool update_particle_cells_after_euler_step=true,const bool assume_particles_in_correct_blocks=true,const bool enforce_domain_boundaries=true);
    template<class T_ARRAYS_PARTICLES> void Euler_Step_Particles_Threaded(RANGE<TV_INT>& domain,const ARRAY<T,FACE_INDEX<TV::m> >& V,T_ARRAYS_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,
        const bool assume_particles_in_correct_blocks=true,const bool enforce_domain_boundaries=true);
    template<class T_ARRAYS_PARTICLES> void Second_Order_Runge_Kutta_Step_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,T_ARRAYS_PARTICLES& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true);
    template<class T_ARRAYS_PARTICLES> void Second_Order_Runge_Kutta_Step_Particles_Threaded(RANGE<TV_INT>& domain,const ARRAY<T,FACE_INDEX<TV::m> >& V,T_ARRAYS_PARTICLES& particles,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time,const bool verbose=true);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells(T_ARRAYS_PARTICLES& particles);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Part_One_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block);
#ifdef USE_PTHREADS
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Part_Two_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process,ARRAY<pthread_mutex_t>* mutexes,const ARRAY<int,TV_INT>* domain_index);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Part_Three_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process,ARRAY<pthread_mutex_t>* mutexes,const ARRAY<int,TV_INT>* domain_index);
#else
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Part_Two_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process,void* mutexes,const void* domain_index);
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Part_Three_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process,void* mutexes,const void* domain_index);
#endif
    template<class T_ARRAYS_PARTICLES> void Update_Particle_Cells_Find_List(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles,const ARRAY<ARRAY<int>,TV_INT>& number_of_particles_per_block,ARRAY<ARRAY<TRIPLE<TV_INT,typename T_ARRAYS_PARTICLES::ELEMENT,int> >,TV_INT>& list_to_process);
    template<class T_ARRAYS_PARTICLES> void Delete_Marked_Particles(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLES& particles);
    template<class T_ARRAYS_PARTICLES> void Consistency_Check(RANGE<TV_INT> domain,T_ARRAYS_PARTICLES& particles);
    int Reseed_Particles(const T time,ARRAY<bool,TV_INT>* cell_centered_mask=0);
    int Reseed_Delete_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign);
    void Identify_Escaped_Particles(RANGE<TV_INT>& domain,const int sign=0);
    void Delete_All_Particles_In_Cell(const TV_INT& block);
    void Delete_Deep_Escaped_Particles(const T radius_fraction=1.5,const bool need_to_identify_escaped_particles=true,const bool verbose=false);
    void Delete_Particles_Outside_Grid();
    void Delete_Particles_Far_Outside_Grid(RANGE<TV_INT>& domain);
    void Delete_Particles_Near_Outside_Grid(RANGE<TV_INT>& domain,int axis,int side);
    void Delete_Particles_In_Local_Maximum_Phi_Cells(const int sign);
    void Delete_Particles_In_Local_Maximum_Phi_Cells_Threaded(RANGE<TV_INT>& domain,const int sign,const T tolerance);
    void Identify_And_Remove_Escaped_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,const T radius_fraction=1.5,const T time=0,const bool verbose=true);
    void Identify_And_Remove_Escaped_Particles_Threaded(RANGE<TV_INT>& domain,const ARRAY<T,FACE_INDEX<TV::m> >& V,const T radius_fraction=1.5,const T time=0,const bool verbose=true);
    void Identify_And_Remove_Escaped_Positive_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,const T radius_fraction=1.5,const T time=0,const bool verbose=true);
    void Reincorporate_Removed_Particles(const T radius_fraction,const T mass_scaling,ARRAY<T,FACE_INDEX<TV::m> >* V,const bool conserve_momentum_for_removed_negative_particles=true);
    void Reincorporate_Removed_Particles_Threaded(RANGE<TV_INT>& domain,const T radius_fraction,const T mass_scaling,ARRAY<T,FACE_INDEX<TV::m> >* V);
    bool Fix_Momentum_With_Escaped_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,const ARRAY<T,TV_INT>& momentum_lost,const T radius_fraction,const T mass_scaling,const T time,const bool force=true);
    bool Fix_Momentum_With_Escaped_Particles(const TV& location,const ARRAY<T,FACE_INDEX<TV::m> >& V,const T momentum_lost,const T radius_fraction,const T mass_scaling,const T time,const bool force=true);

    void Delete_Particles_Far_From_Interface(const int discrete_band=4);
    void Delete_Particles_Far_From_Interface_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface);
    void Delete_Particles_Far_From_Interface_Part_Two(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface,const int discrete_band=4);
    void Delete_Particles_Far_From_Interface_Part_Three(RANGE<TV_INT>& domain,T_ARRAYS_CHAR& near_interface);
    void Add_Negative_Particle(const TV& particle_location,const TV& particle_velocity,const unsigned short quantized_collision_distance);
    void Exchange_Overlap_Particles();
    template<class T_ARRAYS_PARTICLES> void Modify_Levelset_Using_Escaped_Particles_Threaded(RANGE<TV_INT>& domain,ARRAY<T,TV_INT>& phi,T_ARRAYS_PARTICLES& particles,ARRAY<T,FACE_INDEX<TV::m> >* V,const int sign,const T one_over_radius_multiplier);
protected:
#ifdef USE_PTHREADS
    bool Attract_Individual_Particle_To_Interface_And_Adjust_Radius(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles,
        const T phi_min,const T phi_max,const BLOCK_UNIFORM<TV>& block,const int index,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const bool delete_particles_that_leave_original_cell,
        RANDOM_NUMBERS<T>& local_random,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>& list_to_process,ARRAY<pthread_mutex_t>* mutexes,const ARRAY<int,TV_INT>* domain_index);
#else
    bool Attract_Individual_Particle_To_Interface_And_Adjust_Radius(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles,
        const T phi_min,const T phi_max,const BLOCK_UNIFORM<TV>& block,const int index,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const bool delete_particles_that_leave_original_cell,
        RANDOM_NUMBERS<T>& local_random,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>& list_to_process,void* mutexes,const void* domain_index);
#endif
    void Adjust_Particle_Radii(const BLOCK_UNIFORM<TV>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int sign);
    template<class T_ARRAYS_PARTICLES> void Modify_Levelset_Using_Escaped_Particles(ARRAY<T,TV_INT>& phi,T_ARRAYS_PARTICLES& particles,ARRAY<T,FACE_INDEX<TV::m> >* V,const int sign);
    int Reseed_Add_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& other_particles,const int sign,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,ARRAY<bool,TV_INT>* cell_centered_mask);
    void Reseed_Add_Particles_Threaded_Part_One(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& other_particles,const int sign,ARRAY<bool,TV_INT>* cell_centered_mask,ARRAY<int,TV_INT>& number_of_particles_to_add);
#ifdef USE_PTHREADS
    void Reseed_Add_Particles_Threaded_Part_Two(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const ARRAY<int,TV_INT>& number_of_particles_to_add,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>& list_to_process,ARRAY<pthread_mutex_t>* mutexes,const ARRAY<int,TV_INT>* domain_index);
#else
    void Reseed_Add_Particles_Threaded_Part_Two(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const int sign,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time,const ARRAY<int,TV_INT>& number_of_particles_to_add,ARRAY<ARRAY<TRIPLE<TV_INT,PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT>& list_to_process,void* mutexes,const void* domain_index);
#endif
    //void Copy_From_Move_List_Threaded(RANGE<TV_INT>& domain,T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,ARRAY<ARRAY<TRIPLE<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT,int> >,TV_INT>& move_particles,ARRAY<pthread_mutex_t>* mutexes,const ARRAY<int,TV_INT>* domain_index);
    void Identify_Escaped_Particles(const BLOCK_UNIFORM<TV>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign);
    int Delete_Deep_Escaped_Particles(const BLOCK_UNIFORM<TV>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,ARRAY<bool>& escaped,const int sign,const T radius_fraction,
        const bool need_to_identify_escaped_particles);
    void Delete_Particles_Outside_Grid(const T domain_boundary,const int axis,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int side);
    void Delete_Particles_Outside_Grid(const T domain_boundary,const int axis,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int side);
    template<class T_FACE_LOOKUP_LOOKUP> int Remove_Escaped_Particles(const BLOCK_UNIFORM<TV>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& removed_particles,const T_FACE_LOOKUP_LOOKUP& V,const T radius_fraction,const T time);
    int Remove_Escaped_Particles(const BLOCK_UNIFORM<TV>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const ARRAY<bool>& escaped,const int sign,const T radius_fraction);
    void Reincorporate_Removed_Particles(const BLOCK_UNIFORM<TV>& block,PARTICLE_LEVELSET_PARTICLES<TV>*& particles,const int sign,
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& removed_particles,const T radius_fraction,const T mass_scaling,ARRAY<T,FACE_INDEX<TV::m> >* V);
//#####################################################################
};
}
#endif
