//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Ray_Tracing/Rendering/PHOTON_MAP.h>
#include <Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <Ray_Tracing/Rendering/RENDERING_RAY.h>
#include <Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <Ray_Tracing/Rendering_Lights/RENDERING_RECTANGLE_LIGHT.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_TRIANGULATED_SURFACE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> RENDERING_RECTANGLE_LIGHT<T>::
RENDERING_RECTANGLE_LIGHT(const VECTOR<T,3>& position_input,const VECTOR<T,3>& color_input,
    const T brightness_input,const VECTOR<T,3>& u_direction_input,const VECTOR<T,3>& v_direction_input,
    const int u_samples_input,const int v_samples_input,RENDER_WORLD<T>& world_input,const bool supports_global_photons,
    const bool supports_caustic_photons,const bool supports_volume_photons,const bool use_stratified_sampling_on_photon_emit)
    :RENDERING_LIGHT<T>(position_input,color_input,brightness_input,world_input,supports_global_photons,
        supports_caustic_photons,supports_volume_photons),
    RENDERING_TRIANGULATED_SURFACE<T>(*(TRIANGULATED_SURFACE<T>::Create())),
    u_direction(u_direction_input),v_direction(v_direction_input),
    u_samples(u_samples_input),v_samples(v_samples_input),
    use_stratified_sampling_on_photon_emit(use_stratified_sampling_on_photon_emit)
{
    triangulated_surface.Clean_Memory();
    one_over_u_samples=1/T(u_samples);
    one_over_v_samples=1/T(v_samples);
    normal=VECTOR<T,3>::Cross_Product(u_direction,v_direction).Normalized();
    GEOMETRY_PARTICLES<VECTOR<T,3> >& particles=triangulated_surface.particles;
    particles.Preallocate(4);
    int index_1=particles.Add_Element();int index_2=particles.Add_Element();
    int index_3=particles.Add_Element();int index_4=particles.Add_Element();
    particles.X(index_1)=VECTOR<T,3>(position);
    particles.X(index_2)=VECTOR<T,3>(position+u_direction);
    particles.X(index_3)=VECTOR<T,3>(position+v_direction);
    particles.X(index_4)=VECTOR<T,3>(position+u_direction+v_direction);
    ARRAY<VECTOR<int,3> > triangles(2);triangles(0).Set(index_1,index_2,index_3);triangles(1).Set(index_2,index_4,index_3);
    triangulated_surface.mesh.Initialize_Mesh(particles.Size(),triangles);
    triangulated_surface.Update_Triangle_List();
    triangulated_surface.Initialize_Hierarchy(true);
    triangulated_surface.Update_Bounding_Box();
    triangulated_surface.Use_Vertex_Normals();
    triangulated_surface.Update_Vertex_Normals();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> RENDERING_RECTANGLE_LIGHT<T>::
~RENDERING_RECTANGLE_LIGHT()
{
    delete &triangulated_surface;
}
//#####################################################################
// Function Sample_Points
//#####################################################################
template<class T> void RENDERING_RECTANGLE_LIGHT<T>::
Sample_Points(const VECTOR<T,3>& surface_position,const VECTOR<T,3>& surface_normal,
    ARRAY<RAY<VECTOR<T,3> > >& sample_array)const
{
    sample_array.Resize(u_samples*v_samples);
    for(int u=0;u<u_samples;u++)
        for(int v=0;v<v_samples;v++){
            T u_jitter_offset=sample_points_random.Get_Uniform_Number(T(-.5),T(.5)),v_jitter_offset=sample_points_random.Get_Uniform_Number(T(-.5),T(.5));
            sample_array(u+u_samples*v+1)=
                RAY<VECTOR<T,3> >(SEGMENT_3D<T>(surface_position,position+(u+u_jitter_offset+T(.5))*one_over_u_samples*u_direction+(v+v_jitter_offset+T(.5))*one_over_v_samples*v_direction));}
}
//#####################################################################
// Function Emitted_Light
//#####################################################################
template<class T> VECTOR<T,3> RENDERING_RECTANGLE_LIGHT<T>::
Emitted_Light(const RENDERING_RAY<T>& ray) const
{
    T cosine_of_angle=-VECTOR<T,3>::Dot_Product(ray.ray.direction,normal);
    if(cosine_of_angle<0)return VECTOR<T,3>(0,0,0);
    else return (1/((T)pi*4))*cosine_of_angle*color*brightness/sqr(ray.ray.t_max);
}
//#####################################################################
// Function Emit_Photons
//#####################################################################
template<class T> int RENDERING_RECTANGLE_LIGHT<T>::
Emit_Photons(RENDERING_RAY<T>& parent_ray,PHOTON_MAP<T>& photon_map,const typename PHOTON_MAP<T>::PHOTON_MAP_TYPE type) const 
{
    int number_emitted=0;
    RANDOM_NUMBERS<T>* random=0;
    if(type==PHOTON_MAP<T>::GLOBAL_PHOTON_MAP) random=&global_photon_random;
    else if(type==PHOTON_MAP<T>::CAUSTIC_PHOTON_MAP) random=&caustic_photon_random;
    else if(type==PHOTON_MAP<T>::VOLUME_PHOTON_MAP) random=&volume_photon_random;
    else PHYSBAM_FATAL_ERROR();
    if(use_stratified_sampling_on_photon_emit){
        const int theta_samples=10,phi_samples=10;
        T one_over_theta_samples=T(1)/T(theta_samples),one_over_phi_samples=T(1)/T(phi_samples);
        while(true)
            for(int u=0;u<u_samples;u++)
                for(int v=0;v<v_samples;v++)
                    for(int theta=0;theta<theta_samples;theta++)
                        for(int phi=0;phi<phi_samples;phi++){
                            if(!photon_map.Light_Emission_Quota_Remains())return number_emitted;
                            world.random.Set_Seed(abs((int)random->Get_Uniform_Number((T)0,(T)1)));
                            T u_jitter_offset=world.random.Get_Uniform_Number((T)-.5,(T).5),v_jitter_offset=world.random.Get_Uniform_Number(T(-.5),T(.5));
                            T xi_1_jitter_offset=world.random.Get_Uniform_Number((T)-.5,(T).5),xi_2_jitter_offset=world.random.Get_Uniform_Number(T(-.5),T(.5));
                            T xi_1=(xi_1_jitter_offset+T(.5)+theta)*one_over_theta_samples,xi_2=(xi_2_jitter_offset+T(.5)+phi)*one_over_phi_samples;
                            VECTOR<T,3> point_on_light=position+(u+u_jitter_offset+T(.5))*one_over_u_samples*u_direction+(v+v_jitter_offset+T(.5))*one_over_v_samples*v_direction;
                            VECTOR<T,3> robust_point_on_light=normal*small_number*4+point_on_light;
                            VECTOR<T,3> direction=world.Random_Reflected_Direction_In_Hemisphere(normal,xi_1,xi_2);
                            RENDERING_RAY<T> photon_reflect_ray(RAY<VECTOR<T,3> >(robust_point_on_light,direction),1,parent_ray.current_object);
                            world.Cast_Photon(photon_reflect_ray,parent_ray,color*brightness,type,0,0);
                            number_emitted++;}}
    else while(photon_map.Light_Emission_Quota_Remains()){ 
            world.random.Set_Seed(abs((int)random->Get_Uniform_Number((T)0,(T)1)));
            VECTOR<T,3> point_on_light=normal*small_number*4+position+world.random.Get_Uniform_Number(T(0),T(1))*u_direction+world.random.Get_Uniform_Number(T(0),T(1))*v_direction;
            VECTOR<T,3> direction=world.Random_Reflected_Direction_In_Hemisphere(normal,world.random.Get_Uniform_Number(T(0),T(1)),world.random.Get_Uniform_Number(T(0),T(1)));
            RENDERING_RAY<T> photon_emit_ray(RAY<VECTOR<T,3> >(point_on_light,direction),1,parent_ray.current_object);
            world.Cast_Photon(photon_emit_ray,parent_ray,color*brightness,type,0,0);
            number_emitted++;}
    return number_emitted;
}
template class RENDERING_RECTANGLE_LIGHT<double>;
template class RENDERING_RECTANGLE_LIGHT<float>;
}
