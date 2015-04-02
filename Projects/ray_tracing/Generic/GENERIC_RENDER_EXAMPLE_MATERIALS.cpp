//#####################################################################
// Copyright 2004-2009, Zhaosheng Bao, Eilene Hao, Sergey Koltakov, Nipun Kwatra, Frank Losasso, Andrew Selle, Jerry Talton, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt. 
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <Ray_Tracing/Rendering_Lights/RENDERING_LIGHT.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_IMPLICIT_SURFACE.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_LEVELSET_MULTIPLE_OBJECT.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_ABSORPTION_SPECTRUM_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BLEND_IMPLICIT_SURFACE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BLEND_TRIANGULATED_SURFACE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_IMAGE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_CAMERA_ORIENTED_NORMAL_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_CHECKERBOARD_TEXTURE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_COLOR_BLEND_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_FIXED_BLEND_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_HOMOGENEOUS_VOLUME_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_INFINITE_REFLECTION_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_KAJIYA_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_LAMBERTIAN_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_LIGHT_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_MARBLE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_MASKED_BLEND_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_NOISE_TEXTURE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_NORMAL_ILLUSTRATION_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_PHONG_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_PHOTON_SINK.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_REFLECTION_MAP_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_REFLECTION_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_SHELL_EMISSION_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_SUBDIVISION_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_SUM_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_TEXTURE_IMAGE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_TORRANCE_SPARROW_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_TRANSLUCENCY_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_TRANSPARENT_MATERIAL_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_TRANSPARENT_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_TRIANGULATED_SURFACE_VERTEX_COLOR_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_TWO_SIDE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_UNIFORM_COLOR_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_VOXEL_FIRE_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_VOXEL_SHADER.h>
#include <Ray_Tracing/Rendering_Shaders/RENDERING_WOOD_SHADER.h>
#include <cstdio>
#include "GENERIC_RENDER_EXAMPLE.h"
//#include <Rendering_Shaders/RENDERING_HAIR_SHADER.h>
using namespace PhysBAM;
//#####################################################################
// Material - Prepare materials
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Material(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters)
{
    std::string type=parameters.Get_Parameter("Type",std::string("Null"));
    std::string name=parameters.Get_Parameter("Name",std::string("Name"));
    if(type=="Color"){
        TV color=parameters.Get_Parameter("Color",TV(0.5,0.5,0.5));
        shaders.Set(name,new RENDERING_UNIFORM_COLOR_SHADER<T>(color,world));
        LOG::cout<<"Material '"<<name<<"' Color="<<color<<std::endl;}
    else if(type=="Vertex_Color") {
        std::string inputfile=parameters.Get_Parameter("Input_File", std::string("<unknown>"));
        if(inputfile == "<unknown>"){LOG::cerr<<"You have to supply a file of vertex colors (array of colors)!"<<std::endl;exit(1);}
        RENDERING_TRIANGULATED_SURFACE_VERTEX_COLOR_SHADER<T,RW>* vertex_shader=new RENDERING_TRIANGULATED_SURFACE_VERTEX_COLOR_SHADER<T,RW>(world);
        vertex_shader->Initialize(inputfile);
        shaders.Set(name,vertex_shader);
        LOG::cout<<"Material '"<<name<<"' Source File='"<<inputfile<<"'"<<std::endl;}
    else if(type=="Checker"){
        T scaling=parameters.Get_Parameter("Scaling",(T)1);
        std::string shader1name=parameters.Get_Parameter("Shader1",std::string("<unknown>"));
        std::string shader2name=parameters.Get_Parameter("Shader2",std::string("<unknown>"));
        MATERIAL_SHADER<T> *shader1=shaders.Get(shader1name),*shader2=shaders.Get(shader2name);
        if(shader1==0){LOG::cerr<<"Invalid shader 1 "<<shader1name<<" specified in Checker"<<std::endl;exit(1);}
        if(shader2==0){LOG::cerr<<"Invalid shader 2"<<shader2name<<" specified in Checker"<<std::endl;exit(1);}
        shaders.Set(name,new RENDERING_CHECKERBOARD_TEXTURE_SHADER<T>(*shader1,*shader2,scaling,world));
        LOG::cout<<"Checkerboard '"<<name<<"' Scaling="<<scaling<<std::endl;}
    else if(type=="Two_Side"){
        T scaling=parameters.Get_Parameter("Scaling",(T)1);
        std::string frontname=parameters.Get_Parameter("Front_Shader",std::string("<unknown>"));
        std::string backname=parameters.Get_Parameter("Back_Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T> *front=shaders.Get(frontname),*back=shaders.Get(backname);
        if(front==0){LOG::cerr<<"Invalid Front shader"<<frontname<<" specified in Two_Side"<<std::endl;exit(1);}
        if(back==0){LOG::cerr<<"Invalid Back shader"<<backname<<" specified in Two_Side"<<std::endl;exit(1);}
        shaders.Set(name,new RENDERING_TWO_SIDE_SHADER<T>(*front,*back,world));
        LOG::cout<<"Checkerboard '"<<name<<"' Scaling="<<scaling<<std::endl;}
    else if(type=="Sum_Shader"){
        std::cout<<"Material "<<name<<std::endl;
        RENDERING_SUM_SHADER<T>* shader=new RENDERING_SUM_SHADER<T>(world);
        for(int i=1;;i++){
            std::string shader_name=parameters.Get_Parameter(LOG::sprintf("Shader_%d",i),std::string("<unknown>"));
            std::cout<<"   Shader "<<i<<" is "<<shader_name<<std::endl;
            if(shader_name=="<unknown>") break;
            shader->shaders.Append(shaders.Get(shader_name));}
        shaders.Set(name,shader);}
    else if(type=="Marble"){
        std::cout<<"Material "<<name<<std::endl;
        TV color1=parameters.Get_Parameter("Color1",TV(0,0,0));
        TV color2=parameters.Get_Parameter("Color2",TV(0,0,0));
        int octaves=parameters.Get_Parameter("Octaves",(int)10);
        T gain=parameters.Get_Parameter("Gain",(T).5);
        T lacunarity=parameters.Get_Parameter("Lacunarity",(T)1.99);
        T low=parameters.Get_Parameter("Low",(T)0);
        T high=parameters.Get_Parameter("High",(T)1);
        T vain_value=parameters.Get_Parameter("Vain_Value",(T).4);
        T vain_width=parameters.Get_Parameter("Vain_Width",(T).05);
        T clamp_width=parameters.Get_Parameter("Clamp_Width",(T).4);
        RENDERING_MARBLE_SHADER<T>* shader=new RENDERING_MARBLE_SHADER<T>(color1,color2,octaves,lacunarity,gain,low,high,vain_value,vain_width,clamp_width,world);
        shaders.Set(name,shader);}
    else if(type=="Wood"){
        std::cout<<"Material "<<name<<std::endl;
        TV color1=parameters.Get_Parameter("Color1",TV(0,0,0));
        TV color2=parameters.Get_Parameter("Color2",TV(0,0,0));
        T ring_frequency=parameters.Get_Parameter("Ring_Frequency",(T)8);
        T ring_noise=parameters.Get_Parameter("Ring_Noise",(T)0.02);
        T ring_noise_frequency=parameters.Get_Parameter("Ring_Noise_Frequency",(T)1);
        T trunk_wobble=parameters.Get_Parameter("Trunk_Wobble",(T).15);
        T trunk_wobble_frequency=parameters.Get_Parameter("Trunk_Wobble_Frequency",(T)1);
        T angular_wobble=parameters.Get_Parameter("Angular_Wobble",(T)1);
        T angular_wobble_frequency=parameters.Get_Parameter("Angular_Wobble_Frequency",(T)1.5);
        T grain_frequency=parameters.Get_Parameter("Grain_Frequency",(T)25);
        T ringy=parameters.Get_Parameter("Ring_Y",(T)1);
        T grainy=parameters.Get_Parameter("Grain_Y",(T)1);
        RENDERING_WOOD_SHADER<T>* shader=new RENDERING_WOOD_SHADER<T>(color1,color2,ring_frequency,ring_noise,ring_noise_frequency,trunk_wobble,trunk_wobble_frequency,angular_wobble,angular_wobble_frequency,
            grain_frequency,grainy,ringy,world);
        shaders.Set(name,shader);}
    else if(type=="Torrance_Sparrow"){
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* color_shader=shaders.Get(shader_name);
        if(color_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in lambertian"<<std::endl;exit(1);}
        T exponent=parameters.Get_Parameter("Exponent",(T)3);
        T index_of_refraction=parameters.Get_Parameter("Index_Of_Refraction",(T).177); // silver defaults
        T absorption_coefficient=parameters.Get_Parameter("Absorption_Coefficient",(T)3.638); // silver defaults
        MATERIAL_SHADER<T>* shader=new RENDERING_TORRANCE_SPARROW_SHADER<T>(*color_shader,exponent,index_of_refraction,absorption_coefficient,world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Torrance_Sparrow shader="<<shader_name<<" exponent="<<exponent<<std::endl;}
    else if(type=="Lambertian"){
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* color_shader=shaders.Get(shader_name);
        if(color_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in lambertian"<<std::endl;exit(1);}
        TV ambient=parameters.Get_Parameter("Ambient",TV(0,0,0));
        T reflectance=parameters.Get_Parameter("Reflectance",(T)0.5);
        bool visualize_photon_map_directly=parameters.Get_Parameter("Visualize_Photon_Map_Directly",false);
        MATERIAL_SHADER<T>* shader=new RENDERING_LAMBERTIAN_SHADER<T>(*color_shader,1,ambient,reflectance,TV(1,1,1),world,visualize_photon_map_directly);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Lambertian shader="<<shader_name<<std::endl;}
    else if(type=="Phong"){
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* child_shader=shaders.Get(shader_name);
        if(child_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in phong"<<std::endl;exit(1);}
        TV ambient=parameters.Get_Parameter("Ambient",TV(0,0,0));
        TV diffuse=parameters.Get_Parameter("Diffuse",TV((T)0.5,(T)0.5,(T)0.5));
        TV specular=parameters.Get_Parameter("Specular",TV((T)1.0,(T)1.0,(T)1.0));
        T specular_exponent=parameters.Get_Parameter("Specular_Exponent",(T)20);
        MATERIAL_SHADER<T>* shader=new RENDERING_PHONG_SHADER<T>(ambient,diffuse,specular,specular_exponent,*child_shader,world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Phong shader="<<shader_name<<std::endl;}
    else if(type=="Blend"){
        std::string child_name_1=parameters.Get_Parameter("Shader1",std::string("<unknown>"));
        std::string child_name_2=parameters.Get_Parameter("Shader2",std::string("<unknown>"));
        bool direct_shading_only=parameters.Get_Parameter("Direct_Shading_Only",false);
        T blending_fraction=parameters.Get_Parameter("Blend_Fraction",(T)0.5);
        MATERIAL_SHADER<T>* child_shader_1=shaders.Get(child_name_1);
        if(child_shader_1==0){LOG::cerr<<"Invalid shader "<<child_name_1<<" specified in blend shader"<<std::endl;exit(1);}
        MATERIAL_SHADER<T>* child_shader_2=shaders.Get(child_name_2);
        if(child_shader_2==0){LOG::cerr<<"Invalid shader "<<child_name_2<<" specified in blend shader"<<std::endl;exit(1);}
        MATERIAL_SHADER<T>* shader=new RENDERING_FIXED_BLEND_SHADER<T>(blending_fraction,*child_shader_1,*child_shader_2,direct_shading_only,world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Blend shaders="<<child_name_1<<", "<<child_name_2<<std::endl;}
    else if(type=="MaskedBlend"){
        std::string mask_name=parameters.Get_Parameter("MaskShader",std::string("<unknown>"));
        std::string child_name_1=parameters.Get_Parameter("Shader1",std::string("<unknown>"));
        std::string child_name_2=parameters.Get_Parameter("Shader2",std::string("<unknown>"));
        std::string channel=parameters.Get_Parameter("Channel",std::string("gray")); // red, green, blue, or gray.
        bool direct_shading_only=parameters.Get_Parameter("Direct_Shading_Only",false);
        MATERIAL_SHADER<T>* mask_shader=shaders.Get(mask_name);
        if(mask_shader==0){LOG::cerr<<"Invalid shader "<<mask_name<<" specified in blend shader"<<std::endl;exit(1);}
        MATERIAL_SHADER<T>* child_shader_1=shaders.Get(child_name_1);
        if(child_shader_1==0){LOG::cerr<<"Invalid shader "<<child_name_1<<" specified in blend shader"<<std::endl;exit(1);}
        MATERIAL_SHADER<T>* child_shader_2=shaders.Get(child_name_2);
        if(child_shader_2==0){LOG::cerr<<"Invalid shader "<<child_name_2<<" specified in blend shader"<<std::endl;exit(1);}
        RENDERING_MASKED_BLEND_SHADER<T>* shader=new RENDERING_MASKED_BLEND_SHADER<T>(*mask_shader,*child_shader_1,*child_shader_2,direct_shading_only,world);
        if(channel=="gray") shader->Use_Gray();
        else if(channel=="red") shader->Use_Red();
        else if(channel=="green") shader->Use_Green();
        else if(channel=="blue") shader->Use_Blue();
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Blend shaders="<<mask_name<<", "<<child_name_1<<", "<<child_name_2<<std::endl;}
    else if(type=="Camera_Oriented_Normal"){
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* child_shader=shaders.Get(shader_name);
        if(child_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in camera oriented normal shader"<<std::endl;exit(1);}
        MATERIAL_SHADER<T>* shader=new RENDERING_CAMERA_ORIENTED_NORMAL_SHADER<T>(*child_shader,world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Camera-oriented normal shader="<<shader_name<<std::endl;}
    else if (type=="Kajiya"){
        TV diffuse=parameters.Get_Parameter("Diffuse",TV((T)0.5,(T)0.5,(T)0.5));
        TV specular=parameters.Get_Parameter("Specular",TV((T)1.0,(T)1.0,(T)1.0));
        T diffuse_coefficient=parameters.Get_Parameter("Diffuse_Coefficient",(T)10);
        T specular_coefficient=parameters.Get_Parameter("Specular_Coefficient",(T)10);
        T specular_exponent=parameters.Get_Parameter("Specular_Exponent",(T)20);
        MATERIAL_SHADER<T>* shader=new RENDERING_KAJIYA_SHADER<T>(diffuse,diffuse_coefficient,specular,specular_coefficient,specular_exponent,world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Kajiya Diffuse="<<diffuse<<" Specular="<<specular<<std::endl;}
    else if(type=="Transparent"){
        T reflectivity=parameters.Get_Parameter("Reflectivity",(T)1.0);
        bool fresnel=parameters.Get_Parameter("Fresnel",true);
        std::string shift_direction=parameters.Get_Parameter("Shift_Direction",std::string("<unknown>"));
        MATERIAL_SHADER<T>* shader=new RENDERING_TRANSPARENT_SHADER<T>(reflectivity,fresnel,(shift_direction=="ray")?true:false,world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Transparent shader reflectivity="<<reflectivity<<" fresnel="<<fresnel<<std::endl;}
    else if(type=="Transparent_Blend"){
        std::string shader_name=parameters.Get_Parameter("Blend_Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* child_shader=shaders.Get(shader_name);
        if(!child_shader){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in blend"<<std::endl;exit(1);}
        RENDERING_BLEND_SHADER<T>* blend_shader=dynamic_cast<RENDERING_BLEND_SHADER<T>*>(child_shader);
        if(!blend_shader){LOG::cerr<<shader_name<<" is not a blend shader"<<std::endl;exit(1);}
        T reflectivity0=parameters.Get_Parameter("Reflectivity0",(T)1.0);
        T reflectivity1=parameters.Get_Parameter("Reflectivity1",(T)1.0);
        bool fresnel=parameters.Get_Parameter("Fresnel",true);
        std::string shift_direction=parameters.Get_Parameter("Shift_Direction",std::string("<unknown>"));
        MATERIAL_SHADER<T>* shader=new RENDERING_TRANSPARENT_SHADER<T>(reflectivity0,reflectivity1,*blend_shader,fresnel,(shift_direction=="ray")?true:false,world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Transparent blend shader reflectivity="<<reflectivity0<<" "<<reflectivity1<<" blend_shader="<<shader_name<<" fresnel="<<fresnel<<std::endl;}
    else if(type=="Transparent_Material"){
        T reflectivity=parameters.Get_Parameter("Reflectivity",(T)1.0);
        bool fresnel=parameters.Get_Parameter("Fresnel",true);
        bool use_reflected_ray=parameters.Get_Parameter("Use_Reflected_Ray",true);
        bool use_russian_roulette=parameters.Get_Parameter("Use_Russian_Roulette",false);
        std::string shift_direction=parameters.Get_Parameter("Shift_Direction",std::string("<unknown>"));
        std::string shader_name=parameters.Get_Parameter("Color_Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* color_shader=0;
        if(shader_name!="<unknown>"){
            color_shader=shaders.Get(shader_name);
            if(color_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in Transparent_Material"<<std::endl;exit(1);}}
        shader_name=parameters.Get_Parameter("Surface_Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* surface_shader=0;
        if(shader_name!="<unknown>"){
            surface_shader=shaders.Get(shader_name);
            if(surface_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in Transparent_Material"<<std::endl;exit(1);}}
        MATERIAL_SHADER<T>* shader=new RENDERING_TRANSPARENT_MATERIAL_SHADER<T>(reflectivity,fresnel,(shift_direction=="ray")?true:false,use_reflected_ray,use_russian_roulette,world,color_shader,surface_shader);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Transparent shader reflectivity="<<reflectivity<<" fresnel="<<fresnel<<std::endl;}
    else if(type=="Shell_Emission"){
        T shell_thickness=parameters.Get_Parameter("Shell_Thickness",(T).01);
        T shell_radius_of_curvature=parameters.Get_Parameter("Shell_Radius_Of_Curvature",(T).1);
        T shell_amplification_factor=parameters.Get_Parameter("Shell_Amplification_Factor",(T)10);
        TV shell_emission_color=parameters.Get_Parameter("Shell_Emission_Color",TV(0,0,1));
        MATERIAL_SHADER<T>* shader=new RENDERING_SHELL_EMISSION_SHADER<T>(shell_thickness,shell_radius_of_curvature,shell_amplification_factor,shell_emission_color,world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Shell emission shader thick="<<shell_thickness<<" radius_of_curvature="<<shell_radius_of_curvature<<" shell_amplification="<<shell_amplification_factor<<" color="<<shell_emission_color<<std::endl;}
    else if(type=="Image_Texture"){
        std::string filename=parameters.Get_Parameter("Filename",std::string("<unknown>"));
        if(filename == "<unknown>" || filename.length()<4){LOG::cerr<<"You have to supply an image file for texture-mapping!"<<std::endl;exit(1);}
        bool wrap_s=parameters.Get_Parameter("Wrap_S",(bool)true);bool wrap_t=parameters.Get_Parameter("Wrap_T",(bool)true);
        bool cubic_interpolation=parameters.Get_Parameter("Cubic_Interpolation",(bool)true);
        RENDERING_TEXTURE_IMAGE_SHADER<T>* texture_map_shader=new RENDERING_TEXTURE_IMAGE_SHADER<T>(world,wrap_s,wrap_t,cubic_interpolation);
        texture_map_shader->Initialize(filename);
        shaders.Set(name,texture_map_shader);
        LOG::cout<<"Material '"<<name<<"' Filename="<<filename<<std::endl;}
    else if(type=="Translucency"){
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* color_shader=shaders.Get(shader_name);
        if(color_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in lambertian"<<std::endl;exit(1);}
        T translucency=parameters.Get_Parameter("Translucency",(T)0.5);
        MATERIAL_SHADER<T>* shader=new RENDERING_TRANSLUCENCY_SHADER<T>(*color_shader,world,translucency);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Lambertian shader="<<shader_name<<std::endl;}
    else if(type=="Infinite_Reflection"){
        std::string filename=parameters.Get_Parameter("Filename",std::string("<unknown>"));
        if(filename == "<unknown>"){LOG::cerr<<"You have to supply an image file for texture-mapping!"<<std::endl;exit(1);}
        T scale=parameters.Get_Parameter("Scale",(T)1);
        T rotation=parameters.Get_Parameter("Rotation",(T)0)*(T)pi/(T)180;
        RENDERING_INFINITE_REFLECTION_SHADER<T>* texture_map_shader=new RENDERING_INFINITE_REFLECTION_SHADER<T>(world,rotation);
        texture_map_shader->Initialize(filename,scale);
        shaders.Set(name,texture_map_shader);
        LOG::cout<<"Material '"<<name<<"' Filename="<<filename<<std::endl;}
    else if(type=="Reflection_Map"){
        std::string filename=parameters.Get_Parameter("Filename",std::string("<unknown>"));
        if(filename == "<unknown>"){LOG::cerr<<"You have to supply an image file for reflection-mapping!"<<std::endl;exit(1);}
        RENDERING_REFLECTION_MAP_SHADER<T>* reflection_map_shader=new RENDERING_REFLECTION_MAP_SHADER<T>(world);
        reflection_map_shader->Initialize(filename);
        shaders.Set(name,reflection_map_shader);
        LOG::cout<<"Material '"<<name<<"' Filename="<<filename<<std::endl;}
    else if(type=="Image_Bump_Map"){
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* child_shader=shaders.Get(shader_name);
        if(child_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in bump map"<<std::endl;exit(1);}
        T bump_height=parameters.Get_Parameter("Bump_Height",(T).01);
        std::string filename=parameters.Get_Parameter("Filename",std::string("<unknown>"));
        if(filename == "<unknown>"){LOG::cerr<<"You have to supply an image file for bump-mapping!"<<std::endl;exit(1);}
        RENDERING_BUMP_MAP_IMAGE_SHADER<T>* bump_map_shader=new RENDERING_BUMP_MAP_IMAGE_SHADER<T>(child_shader,world);
        bool load_success=bump_map_shader->Initialize(filename,bump_height);
        if(!load_success){LOG::cerr<<"Failed to load texture!"<<std::endl;exit(1);}
        shaders.Set(name,bump_map_shader);
        LOG::cout<<"Material '"<<name<<"' Filename="<<filename<<" Shader="<<shader_name<<std::endl;}
    else if(type=="Normal_Image_Bump_Map"){
        std::string shader_name=parameters.Get_Parameter("Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* child_shader=shaders.Get(shader_name);
//        if(child_shader==0){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in bump map"<<std::endl;exit(1);}
        std::string filename=parameters.Get_Parameter("Filename",std::string("<unknown>"));
        if(filename == "<unknown>"){LOG::cerr<<"You have to supply an image file for bump-mapping!"<<std::endl;exit(1);}
        RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER<T>* bump_map_shader=new RENDERING_BUMP_MAP_NORMAL_IMAGE_SHADER<T>(child_shader,world);
        bump_map_shader->Initialize(filename);
        shaders.Set(name,bump_map_shader);
        LOG::cout<<"Material '"<<name<<"' Filename="<<filename<<" Shader="<<shader_name<<std::endl;}
    else if(type=="Implicit_Surface_Blend"){
        std::string shader1_name=parameters.Get_Parameter("Shader1",std::string("<unknown>")),shader2_name=parameters.Get_Parameter("Shader2",std::string("<unknown>"));
        MATERIAL_SHADER<T> *child_shader1=shaders.Get(shader1_name),*child_shader2=shaders.Get(shader2_name);
        if(!child_shader1){LOG::cerr<<"Invalid shader "<<shader1_name<<" specified in blend"<<std::endl;exit(1);}
        if(!child_shader2){LOG::cerr<<"Invalid shader "<<shader2_name<<" specified in blend"<<std::endl;exit(1);}
        std::string implicit_surface_type=parameters.Get_Parameter("Implicit_Surface_Type",std::string("<unknown>"));
        std::string raw_grid_filename=parameters.Get_Parameter("Grid_Filename",std::string("<unknown_grid_file>"));
        std::string raw_phi_filename=parameters.Get_Parameter("Phi_Filename",std::string("<unknown_phi_file>"));
        IMPLICIT_OBJECT<TV>* implicit_surface;
            {LOG::cout<<"Unknown implicit surface type "<<implicit_surface_type<<std::endl;exit(1);}
        INTERVAL<T> value_range(parameters.Get_Parameter("Low_Value",(T)0),parameters.Get_Parameter("High_Value",(T)1));
        INTERVAL<T> weight_range(parameters.Get_Parameter("Min_Weight",(T)0),parameters.Get_Parameter("Max_Weight",(T)1));
        shaders.Set(name,new RENDERING_BLEND_IMPLICIT_SURFACE_SHADER<T>(*implicit_surface,value_range,weight_range,*child_shader1,*child_shader2,world));
        LOG::cout<<"Material '"<<name<<"' Range="<<value_range<<" Weight_Range="<<weight_range<<" Shader1="<<shader1_name<<" Shader2="<<shader2_name<<std::endl;}
    else if(type=="Triangulated_Surface_Blend"){
        bool direct_shading_only=parameters.Get_Parameter("Direct_Shading_Only",false);
        std::string shader1_name=parameters.Get_Parameter("Shader1",std::string("<unknown>")),shader2_name=parameters.Get_Parameter("Shader2",std::string("<unknown>"));
        MATERIAL_SHADER<T> *child_shader1=shaders.Get(shader1_name),*child_shader2=shaders.Get(shader2_name);
        if(!child_shader1){LOG::cerr<<"Invalid shader "<<shader1_name<<" specified in blend"<<std::endl;exit(1);}
        if(!child_shader2){LOG::cerr<<"Invalid shader "<<shader2_name<<" specified in blend"<<std::endl;exit(1);}
        std::string raw_field_file=parameters.Get_Parameter("Field",std::string("<unknown>"));
        if(raw_field_file == "<unknown>"){LOG::cerr<<"You have to supply a field file for blending."<<std::endl;exit(1);}
        std::string field_file=Animated_Filename(raw_field_file,frame);
        ARRAY<T>* field=new ARRAY<T>;FILE_UTILITIES::Read_From_File<RW>(field_file,*field);
        T low_value=parameters.Get_Parameter("Low_Value",(T)0),high_value=parameters.Get_Parameter("High_Value",(T)1);
        shaders.Set(name,new RENDERING_BLEND_TRIANGULATED_SURFACE_SHADER<T>(*field,low_value,high_value,*child_shader1,*child_shader2,direct_shading_only,world));
        LOG::cout<<"Material '"<<name<<"' Field="<<field_file<<" Range=["<<low_value<<","<<high_value<<"] Shader1="<<shader1_name<<" Shader2="<<shader2_name<<std::endl;}
    else if(type=="Color_Blend"){
        RENDERING_COLOR_BLEND_SHADER<T> *shader=new RENDERING_COLOR_BLEND_SHADER<T>(world);
        shaders.Set(name,shader);
        LOG::cout<<"Material '"<<name<<"' Color blend shader: ";
        for(int i=1;;i++){
            std::string child_shader_tag=LOG::sprintf("Shader%d",i);
            std::string weight_tag=LOG::sprintf("Weight%d",i);
            if(!parameters.Is_Defined(child_shader_tag) || !parameters.Is_Defined(weight_tag)) break;
            std::string child_shader_name=parameters.Get_Parameter(child_shader_tag,std::string("<unknown>"));
            T weight=parameters.Get_Parameter(weight_tag,(T)1);
            MATERIAL_SHADER<T> *child_shader=shaders.Get(child_shader_name);
            if(!child_shader){LOG::cerr<<"Invalid shader "<<child_shader_name<<" specified in blend"<<std::endl;exit(1);}
            LOG::cout<<"('"<<child_shader_name<<"',"<<weight<<") ";
            shader->Add_Shader(child_shader,weight);}
        LOG::cout<<std::endl;}
    else if(type=="Normal_Illustration"){
        shaders.Set(name,new RENDERING_NORMAL_ILLUSTRATION_SHADER<T>(world));
        LOG::cout<<"Material '"<<name<<"' Normal illustration shader"<<std::endl;}
    else if(type=="Noise_Shader"){
        TV color0=parameters.Get_Parameter("First_Color",TV((T)0,(T)0,(T)0));
        TV color1=parameters.Get_Parameter("Second_Color",TV((T)1,(T)1,(T)1));
        T texture_scaling_factor=parameters.Get_Parameter("Texture_Scaling_Factor",(T)1.0);
        MATERIAL_SHADER<T>* shader=new RENDERING_NOISE_TEXTURE_SHADER<T>(color0,color1,texture_scaling_factor,world);
        shaders.Set(name,shader);
        LOG::cout<<"Noise shader '"<<name<<"' first color="<<color0<<" second_color="<<color1<<" texture scaling factor="<<texture_scaling_factor<<std::endl;}
    else if(type=="Subdivision_Shader"){
        std::string shader_name=parameters.Get_Parameter("Material_Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* child_shader=shaders.Get(shader_name);
        if(!child_shader){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in material"<<std::endl;exit(1);}
        RENDERING_SUBDIVISION_SHADER<T>* shader=new RENDERING_SUBDIVISION_SHADER<T>(*child_shader,world);
        shaders.Set(name,shader);
        LOG::cout<<"Subdivision shader '"<<name<<"' child shader='"<<shader_name<<"'"<<std::endl;}
    else{LOG::cerr<<"ERROR: Unknown shader type '"<<type<<"'."<<std::endl;exit(1);}
}
//#####################################################################
// Volume_Material - Prepare materials
//#####################################################################
template<class T,class RW> void GENERIC_RENDER_EXAMPLE<T,RW>::
Volume_Material(RENDER_WORLD<T>& world,const int frame,PARAMETER_LIST& parameters)
{
    std::string type=parameters.Get_Parameter("Type",std::string("Null"));
    std::string name=parameters.Get_Parameter("Name",std::string("Name"));
    if(type=="Voxel_Fire_Shader"){
        LOG::cout<<"Voxel_Fire_Shader"<<std::endl;
        T absorption=parameters.Get_Parameter("Absorption",(T)0);
        T scattering=parameters.Get_Parameter("Scattering",(T)0);
        T inscattering_amplification=parameters.Get_Parameter("Inscattering_Amplification",(T)0);
        T white_point_temperature=parameters.Get_Parameter("White_Point_Temperature",(T)1500);
        T emission_amplification=parameters.Get_Parameter("Emission_Amplification",(T)0);
        T temperature_scale=parameters.Get_Parameter("Temperature_Scale",(T)1);
        T temperature_offset=parameters.Get_Parameter("Temperature_Offset",(T)0);
        bool clamp_low_temperature=parameters.Get_Parameter("Clamp_Low_Temperature",false);
        T temperature_lowest=parameters.Get_Parameter("Temperature_Lowest",(T)0);
        bool use_lms_scaling=parameters.Get_Parameter("Use_LMS_Scaling",(bool)true);
        T density_scale=parameters.Get_Parameter("Density_Scale",(T)1);
        T density_offset=parameters.Get_Parameter("Density_Offset",(T)0);
        bool clamp_low_density=parameters.Get_Parameter("Clamp_Low_Density",false);
        T density_lowest=parameters.Get_Parameter("Density_Lowest",(T)0);
        std::string levelset_object_name=parameters.Get_Parameter("Blue_Core_Levelset",std::string("<unknown>"));
        RENDERING_OBJECT<T>* levelset_object=0;
        LEVELSET<TV>* blue_core_levelset;
        if(!objects.Get(levelset_object_name,levelset_object)){
            fprintf(stderr,"Invalid blue core levelset object named '%s' specified\n",levelset_object_name.c_str());
            blue_core_levelset=0;}
        else{
            RENDERING_IMPLICIT_SURFACE<T>& surface=dynamic_cast<RENDERING_IMPLICIT_SURFACE<T>&>(*levelset_object);
            LEVELSET_IMPLICIT_OBJECT<TV>& levelset_implicit=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>&>(*surface.implicit_surface);
            blue_core_levelset=&(levelset_implicit.levelset);}
        RENDERING_VOXEL_FIRE_SHADER<T>* fire_shader=new RENDERING_VOXEL_FIRE_SHADER<T>(absorption*TV(1,1,1),scattering*TV(1,1,1),
            inscattering_amplification,emission_amplification,temperature_scale,temperature_offset,clamp_low_temperature,temperature_lowest,
            white_point_temperature,density_scale,density_offset,clamp_low_density,density_lowest,use_lms_scaling,blue_core_levelset,world);
        volume_shaders.Set(name,fire_shader);}
    else if(type=="Voxel_Shader"){
        LOG::cout<<"Voxel_Shader"<<std::endl;
        T absorption=parameters.Get_Parameter("Absorption",(T)0);
        T scattering=parameters.Get_Parameter("Scattering",(T)0);
        T inscattering_amplification=parameters.Get_Parameter("Inscattering_Amplification",(T)0);
        T emission_amplification=parameters.Get_Parameter("Emission_Amplification",(T)0);
        T white_point_temperature=parameters.Get_Parameter("White_Point_Temperature",(T)1500);
        bool use_lms_scaling=parameters.Get_Parameter("Use_LMS_Scaling",(bool)false);
        bool use_blackbody_ramp=parameters.Get_Parameter("Use_Blackbody_Ramp",(bool)true);
        bool use_constant_emission_color=parameters.Get_Parameter("Use_Constant_Emission_Color",(bool)false);
        TV constant_emission_color=parameters.Get_Parameter("Constant_Emission_Color",TV(0.5,0.5,0.5));
        RENDERING_VOXEL_SHADER<T>* voxel_shader=new RENDERING_VOXEL_SHADER<T>(absorption*TV(1,1,1),scattering*TV(1,1,1),inscattering_amplification,emission_amplification,white_point_temperature,use_lms_scaling,use_blackbody_ramp,use_constant_emission_color,constant_emission_color,world);
        std::string empty_levelset_object_name=parameters.Get_Parameter("Empty_Levelset",std::string("none"));
        bool use_empty_levelset_for_light_attenuation=parameters.Get_Parameter("Use_Empty_Levelset_For_Light_Attenuation",true);
        RENDERING_OBJECT<T>* empty_levelset_object=0;
        if(!objects.Get(empty_levelset_object_name,empty_levelset_object)){
            if(empty_levelset_object_name!="none"){LOG::cout<<"Invalid levelset object \""<<empty_levelset_object<<"\"."<<std::endl;exit(1);}}
        else{
            IMPLICIT_OBJECT<TV>* implicit_surface=dynamic_cast<RENDERING_IMPLICIT_SURFACE<T>&>(*empty_levelset_object).implicit_surface;
            implicit_surface->Compute_Cell_Minimum_And_Maximum();
            voxel_shader->Use_Empty_Implicit_Surface(*implicit_surface,use_empty_levelset_for_light_attenuation);}

        std::string multiple_levelset_object_name=parameters.Get_Parameter("Multiple_Levelset",std::string("none"));
        int non_empty_region = parameters.Get_Parameter("Non_Empty_Region", 0);
        RENDERING_OBJECT<T>* multiple_levelset_object=0;
        if(!objects.Get(multiple_levelset_object_name,multiple_levelset_object)){
            if(multiple_levelset_object_name!="none"){LOG::cout<<"Invalid multiple levelset object \""<<empty_levelset_object<<"\"."<<std::endl;exit(1);}}
        else{
            RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<TV> >& multiple_levelset=
                dynamic_cast<RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<TV> >&>(*multiple_levelset_object);
            voxel_shader->Use_Levelset_Multiple_Region(&multiple_levelset,non_empty_region);}

        T absorption_shadow=parameters.Get_Parameter("Absorption_Shadow",absorption);
        voxel_shader->Set_Absorption_Shadow(absorption_shadow*TV(1,1,1));
        std::string temperature_remap_filename=parameters.Get_Parameter("Temperature_Remap",std::string("none"));
        if(temperature_remap_filename!="none"){
            ARRAY<T,VECTOR<int,1> > temperature_remap;FILE_UTILITIES::Read_From_File<RW>(temperature_remap_filename,temperature_remap);voxel_shader->Use_Temperature_Remap(temperature_remap);}
        volume_shaders.Set(name,voxel_shader);
        LOG::cout<<"Volume shader '"<<name<<"' Absorption="<<absorption<<"' Absorption_Shadow="<<absorption_shadow<<" Scattering="<<scattering<<" Inscattering_Amplification="<<inscattering_amplification<<" Emission_Amplification="<<emission_amplification
            <<" Use_LMS_Scaling="<<use_lms_scaling<<" Empty_Levelset="<<empty_levelset_object_name<<std::endl;}
    else if(type=="Homogeneous_Volume_Shader"){
        LOG::cout<<"Homogeneous_Volume_Shader"<<std::endl;
        T absorption=parameters.Get_Parameter("Absorption",(T)0);
        T scattering=parameters.Get_Parameter("Scattering",(T)0);
        T volumetric_step=parameters.Get_Parameter("Volume_Step",(T)0.01);
        volume_shaders.Set(name,new RENDERING_HOMOGENEOUS_VOLUME_SHADER<T>(absorption,scattering,volumetric_step,world));}
    else if(type=="Attenuation_Shader"){
        T absorption_coefficient=parameters.Get_Parameter("Absorption",(T)0.5);
        T absorption_clamp=parameters.Get_Parameter("Absorption_Clamp",(T)0);
        TV absorption_spectrum=parameters.Get_Parameter("Spectrum",TV((T)0.5,(T)0.5,(T)0.5));
        volume_shaders.Set(name,new RENDERING_ABSORPTION_SPECTRUM_SHADER<T>(absorption_coefficient,absorption_spectrum,world,absorption_clamp));
        LOG::cout<<"Material '"<<name<<"' Absorption shader absorption="<<absorption_coefficient<<" spectrum="<<absorption_spectrum<<" clamp="<<absorption_clamp<<std::endl;}
    else if(type=="Attenuation_Blend_Shader"){
        std::string shader_name=parameters.Get_Parameter("Blend_Shader",std::string("<unknown>"));
        MATERIAL_SHADER<T>* child_shader=shaders.Get(shader_name);
        if(!child_shader){LOG::cerr<<"Invalid shader "<<shader_name<<" specified in blend"<<std::endl;exit(1);}
        RENDERING_BLEND_SHADER<T>* blend_shader=dynamic_cast<RENDERING_BLEND_SHADER<T>*>(child_shader);
        if(!blend_shader){LOG::cerr<<shader_name<<" is not a blend shader"<<std::endl;exit(1);}
        T absorption_coefficient1=parameters.Get_Parameter("Absorption1",(T)0.5);
        T absorption_coefficient2=parameters.Get_Parameter("Absorption2",(T)0.5);
        TV absorption_spectrum1=parameters.Get_Parameter("Spectrum1",TV((T)0.5,(T)0.5,(T)0.5));
        TV absorption_spectrum2=parameters.Get_Parameter("Spectrum2",TV((T)0.5,(T)0.5,(T)0.5));
        volume_shaders.Set(name,new RENDERING_ABSORPTION_SPECTRUM_SHADER<T>(absorption_coefficient1,absorption_coefficient2,absorption_spectrum1,absorption_spectrum2,*blend_shader,world));
        LOG::cout<<"Material '"<<name<<"' Absorption blend shader absorption="<<absorption_coefficient1<<" "<<absorption_coefficient2<<" spectrum="<<absorption_spectrum1<<" "
            <<absorption_spectrum2<<" blend shader="<<shader_name<<std::endl;}
}
//#####################################################################
template void GENERIC_RENDER_EXAMPLE<float,float>::Material(RENDER_WORLD<float>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<float,float>::Volume_Material(RENDER_WORLD<float>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<double,float>::Material(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<double,float>::Volume_Material(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<double,double>::Material(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
template void GENERIC_RENDER_EXAMPLE<double,double>::Volume_Material(RENDER_WORLD<double>& world,const int frame,PARAMETER_LIST& parameters);
