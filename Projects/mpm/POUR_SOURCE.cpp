//#####################################################################
// Copyright 2016, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Seeding/POISSON_DISK.h>
#include "POUR_SOURCE.h"
#include "STANDARD_TESTS_BASE.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> POUR_SOURCE<TV>::
POUR_SOURCE(STANDARD_TESTS_BASE<TV>& example,IMPLICIT_OBJECT<TV>& seed_area,
    const TV& release_plane_normal,const TV& release_plane_pt,const TV& velocity,const TV& gravity,
    T max_shift,T seeding_buffer,T mass,T volume)
    :example(example),seed_area(seed_area),velocity(velocity),normal(release_plane_normal),gravity(gravity),
    p(release_plane_pt),max_shift(max_shift),buffer_capacity(seeding_buffer),buffer_left(-1),mass(mass),
    volume(volume),cur_time(0),output_file("spout"),show_waiting_particles(true),next_color(1)
{
    T vn=velocity.Dot(normal),gn=gravity.Dot(normal);
    PHYSBAM_ASSERT(vn>0 && gn>=0);

    func_V=[=](const TV& X)
        {
            T a=sqrt(vn*vn+2*gn*((X-p).Dot(normal)));
            return velocity+(a-vn)/gn*gravity;
        };
    func_dV=[=](const TV& X)
        {
            T a=sqrt(vn*vn+2*gn*((X-p).Dot(normal)));
            return MATRIX<T,TV::m>::Outer_Product(gravity/a,normal);
        };
    Refill(true);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> POUR_SOURCE<TV>::
~POUR_SOURCE()
{
    delete &seed_area;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void POUR_SOURCE<TV>::
Write_Output_Files(const int frame)
{
    std::string file=LOG::sprintf("%s/%d/%s",example.output_directory.c_str(),frame,output_file.c_str());
    FILE_UTILITIES::Write_To_File(example.stream_type,file,waiting_particles,buffer_left,cur_time,waiting_particle_color,next_color,random);
    if(show_waiting_particles){
        T vn=velocity.Dot(normal);
        T speed=velocity.Magnitude();
        for(int i=0;i<waiting_particles.m;i++){
            TV X=waiting_particles(i),V=velocity;
            T dot=(X-p).Dot(normal);
            if(dot>=0){
                T t=-dot/vn;
                V=velocity+gravity*t;
                X+=(T).5*t*t*gravity;}
            Add_Debug_Particle(X,dot<=0?VECTOR<T,3>(.5,.5,.5):waiting_particle_color(i));
            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,V);}}
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void POUR_SOURCE<TV>::
Read_Output_Files(const int frame)
{
    std::string file=LOG::sprintf("%s/%d/%s",example.output_directory.c_str(),frame,output_file.c_str());
    FILE_UTILITIES::Read_From_File(example.stream_type,file,waiting_particles,buffer_left,cur_time,waiting_particle_color,next_color,random);
}
//#####################################################################
// Function Begin_Time_Step
//#####################################################################
template<class TV> void POUR_SOURCE<TV>::
Begin_Time_Step(const T time)
{
    Emit();
}
//#####################################################################
// Function End_Time_Step
//#####################################################################
template<class TV> void POUR_SOURCE<TV>::
End_Time_Step(const T time)
{
    Advance_To_Time(time);
}
//#####################################################################
// Function Refill
//#####################################################################
template<class TV> void POUR_SOURCE<TV>::
Refill(bool init)
{
    POISSON_DISK<TV> poisson_disk(1);
    poisson_disk.Set_Distance_By_Volume(volume);
    poisson_disk.Sample(random,seed_area,waiting_particles);
    buffer_left=buffer_capacity;
    if(!init)
        for(int i=waiting_particles.m-1;i>=waiting_particle_color.m;i--)
            if((waiting_particles(i)-p).Dot(normal)>=0)
                waiting_particles.Remove_Index_Lazy(i);
    waiting_particle_color.Resize(waiting_particles.m,true,true,VECTOR<T,3>(next_color&1,(next_color>>1)&1,next_color>>2));
    next_color=next_color%6+1;
}
//#####################################################################
// Function Emit
//#####################################################################
template<class TV> void POUR_SOURCE<TV>::
Emit()
{
    T vn=velocity.Dot(normal),gn=gravity.Dot(normal);
    T speed=velocity.Magnitude();
    for(int i=waiting_particles.m-1;i>=0;i--){
        TV X=waiting_particles(i);
        T dot=(X-p).Dot(normal);
        if(dot<0) continue;
        waiting_particles.Remove_Index_Lazy(i);
        waiting_particle_color.Remove_Index_Lazy(i);
        T t=dot/vn;
        TV X1=X+(T).5*t*t*gravity;
        T a=sqrt(vn*vn+2*gn*((X1-p).Dot(normal)));
        example.Add_Particle(X1,func_V,func_dV,mass,volume);}
}
//#####################################################################
// Function Advance_To_Time
//#####################################################################
template<class TV> void POUR_SOURCE<TV>::
Advance_To_Time(T time)
{
    if(buffer_left<max_shift) Refill(false);
    T dt=time-cur_time;
    waiting_particles+=dt*velocity;
    buffer_left-=dt*velocity.Dot(normal);
    cur_time=time;
}
template class POUR_SOURCE<VECTOR<float,2> >;
template class POUR_SOURCE<VECTOR<float,3> >;
template class POUR_SOURCE<VECTOR<double,2> >;
template class POUR_SOURCE<VECTOR<double,3> >;
}


