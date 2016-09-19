//#####################################################################
// Copyright 2016, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __POUR_SOURCE__
#define __POUR_SOURCE__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <functional>
namespace PhysBAM{

template<class TV> class STANDARD_TESTS_BASE;
template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class POUR_SOURCE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    STANDARD_TESTS_BASE<TV>& example;
    IMPLICIT_OBJECT<TV>& seed_area;
    ARRAY<TV> waiting_particles;
    TV velocity,normal;
    TV gravity,p;
    T max_shift;
    T buffer_capacity;
    T buffer_left;
    T mass;
    T volume;
    T cur_time;
    std::string output_file;
    std::function<TV(const TV&)> func_V;
    std::function<MATRIX<T,TV::m>(const TV&)> func_dV;
    bool show_waiting_particles;
    ARRAY<VECTOR<T,3> > waiting_particle_color;
    int next_color;
    RANDOM_NUMBERS<T> random;

    POUR_SOURCE(STANDARD_TESTS_BASE<TV>& example,IMPLICIT_OBJECT<TV>& seed_area,
        const TV& release_plane_normal,const TV& release_plane_pt,const TV& velocity,
        const TV& gravity,T max_shift,T seeding_buffer,T mass,T volume);
    ~POUR_SOURCE();
    void Write_Output_Files(const int frame);
    void Read_Output_Files(const int frame);
    void Begin_Time_Step(const T time);
    void End_Time_Step(const T time);
    void Refill(bool init);
    void Emit();
    void Advance_To_Time(T time);
};
}

#endif
