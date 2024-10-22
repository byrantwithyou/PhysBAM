//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET
//##################################################################### 
#ifndef __PARTICLE_LEVELSET__
#define __PARTICLE_LEVELSET__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Particles/POINTS_POOL.h>
#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <Incompressible/Level_Sets/LEVELSET_COLLIDABLE.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <climits>
namespace PhysBAM{

template<class TV> class LEVELSET_COLLIDABLE;

template<class TV>
class PARTICLE_LEVELSET
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;typedef typename GRID<TV>::BLOCK T_BLOCK;
public:
    int number_particles_per_cell;
    bool only_use_negative_particles;
    T half_band_width;
    T minimum_particle_radius,maximum_particle_radius;
    T outside_particle_distance_multiplier; // how far outside a particle needs to be before it counts as outside
    int maximum_iterations_for_attraction;
    RANDOM_NUMBERS<T> random;
    bool use_removed_negative_particles,use_removed_positive_particles;
    bool bias_towards_negative_particles; // in the error correction process
    bool store_unique_particle_id;
    int last_unique_particle_id;
    PARTICLE_LEVELSET_PARTICLES<TV> template_particles; // defines the set of attributes used for all particles
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> template_removed_particles; // defines the set of attributes used for all remove particles
    POINTS_POOL<PARTICLE_LEVELSET_PARTICLES<TV> > particle_pool;
    T min_collision_distance_factor,max_collision_distance_factor,max_minus_min_collision_distance_factor_over_max_short;
    bool reincorporate_removed_particles_everywhere;
    bool save_removed_particle_times;
    ARRAY<PAIR<int,T> > removed_particle_times;
    T velocity_interpolation_collidable_contour_value;
    bool delete_positive_particles_crossing_bodies;
    T cfl_number;
    int number_of_ghost_cells;

    LEVELSET_COLLIDABLE<TV> levelset;
    ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*,TV_INT> positive_particles,negative_particles;
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> removed_negative_particles,removed_positive_particles;
    ARRAY<ARRAY<bool>,TV_INT> escaped_positive_particles,escaped_negative_particles;
    ARRAY<ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >,TV_INT> deletion_list;

    PARTICLE_LEVELSET(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,const int number_of_ghost_cells_input);
    PARTICLE_LEVELSET(const PARTICLE_LEVELSET&) = delete;
    void operator=(const PARTICLE_LEVELSET&) = delete;
    virtual ~PARTICLE_LEVELSET();

    void Set_Band_Width(const T number_of_cells=6)
    {half_band_width=number_of_cells*(T).5*levelset.grid.dX.Min();
    levelset.Set_Band_Width(max((T)8,number_of_cells+2));}

    void Set_Number_Particles_Per_Cell(const int number,const int particle_pool_size=0)
    {number_particles_per_cell=number;
    if(particle_pool_size){assert(particle_pool_size>=number);particle_pool.Set_Number_Particles_Per_Cell(particle_pool_size);}
    else particle_pool.Set_Number_Particles_Per_Cell((int)((T)1.5*number));}

    void Use_Removed_Negative_Particles(const bool use_removed_negative_particles_input=true)
    {
        if(use_removed_negative_particles!=use_removed_negative_particles_input){
            use_removed_negative_particles=use_removed_negative_particles_input;
            removed_negative_particles.Delete_Pointers_And_Clean_Memory();
            if(use_removed_negative_particles) removed_negative_particles.Resize(levelset.grid.Block_Indices(number_of_ghost_cells));}
    }
    
    void Use_Removed_Positive_Particles(const bool use_removed_positive_particles_input=true)
    {
        if(use_removed_positive_particles!=use_removed_positive_particles_input){
            use_removed_positive_particles=use_removed_positive_particles_input;
            removed_positive_particles.Delete_Pointers_And_Clean_Memory();
            if(use_removed_positive_particles) removed_positive_particles.Resize(levelset.grid.Block_Indices(number_of_ghost_cells));}
    }

    void Set_Minimum_Particle_Radius(const T radius)
    {minimum_particle_radius=radius;}

    void Set_Maximum_Particle_Radius(const T radius)
    {maximum_particle_radius=radius;}

    void Set_Outside_Particle_Distance(const T fraction_of_the_radius=1)
    {outside_particle_distance_multiplier=fraction_of_the_radius;}
    
    void Set_Maximum_Iterations_For_Attraction(const int maximum_iterations=15)
    {maximum_iterations_for_attraction=maximum_iterations;}

    void Bias_Towards_Negative_Particles(const bool bias_towards_negative_particles_input=true)
    {bias_towards_negative_particles=bias_towards_negative_particles_input;}

    void Set_Collision_Distance_Factors(const T min_collision_distance_factor_input=(T).1,const T max_collision_distance_factor_input=(T)1)
    {min_collision_distance_factor=min_collision_distance_factor_input;max_collision_distance_factor=max_collision_distance_factor_input;
    max_minus_min_collision_distance_factor_over_max_short=(max_collision_distance_factor-min_collision_distance_factor)/USHRT_MAX;}

    T Particle_Collision_Distance(const unsigned short quantized_collision_distance)
    {return (min_collision_distance_factor+(T)quantized_collision_distance*max_minus_min_collision_distance_factor_over_max_short)*levelset.grid.dX.Min();}

    void Set_Velocity_Interpolation_Collidable(const T contour_value_input=0)
    {velocity_interpolation_collidable_contour_value=contour_value_input;}

    static const std::string& Particle_Type_Name(PARTICLE_LEVELSET_PARTICLE_TYPE particle_type)
    {static std::string particle_type_names[4]={"positive","negative","removed positive","removed negative"};return particle_type_names[(int)particle_type];}

//#####################################################################
    void Store_Unique_Particle_Id(const bool store_unique_particle_id_input=true);
    template<class T_PARTICLES> void Compact_Particles(T_PARTICLES& particles);
    int Add_Particle(PARTICLE_LEVELSET_PARTICLES<TV>*& particles);
    int Add_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& particles);
    void Add_Particle_To_Deletion_List(ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >& deletion_list,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index);
    void Add_Particle_To_Deletion_List(ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >& deletion_list,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int index);
    PARTICLE_LEVELSET_PARTICLES<TV>* Allocate_Particles(const PARTICLE_LEVELSET_PARTICLES<TV>& clone_particles);
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* Allocate_Particles(const PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& clone_particles);
    void Delete_Particle_And_Clean_Memory(PARTICLE_LEVELSET_PARTICLES<TV>* head_particles,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index);
    void Delete_Particle_And_Clean_Memory(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* head_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int index);
private:
    bool Delete_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index);
    bool Delete_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int index);
public:
    void Delete_Particles_From_Deletion_List(ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >& deletion_list,PARTICLE_LEVELSET_PARTICLES<TV>& particles,bool verbose=false);
    void Delete_Particles_From_Deletion_List(ARRAY<PAIR<PARTICLE_LEVELSET_PARTICLES<TV>*,int> >& deletion_list,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles);
    PARTICLE_LEVELSET_PARTICLES<TV>& Get_Particle_Link(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int absolute_index,int& index_in_link);
    void Move_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_PARTICLES<TV>& to_particles,const int from_absolute_index);
    void Move_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& to_particles,const int from_absolute_index);
    void Move_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_PARTICLES<TV>& to_particles,const int from_index);
    void Move_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& to_particles,const int from_index);
    void Copy_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_PARTICLES<TV>& to_particles,const int from_absolute_index);
    void Copy_Particle(PARTICLE_LEVELSET_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& to_particles,const int from_absolute_index);
    void Copy_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_PARTICLES<TV>& to_particles,const int from_index);
    void Copy_Particle(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& from_particles,PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& to_particles,const int from_index);
    template<class T_PARTICLES> T_PARTICLES* Allocate_Particle();
    void Free_Particle_And_Clear_Pointer(PARTICLE_LEVELSET_PARTICLES<TV>*& particles);
    void Free_Particle_And_Clear_Pointer(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*& particles);
    bool Adjust_Particle_For_Objects(TV& X,TV& V,const T r,const T collision_distance,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) const;
    void Initialize_Particle_Levelset_Grid_Values();
//#####################################################################
};   
}
#endif

