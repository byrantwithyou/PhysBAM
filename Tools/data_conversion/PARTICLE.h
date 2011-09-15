//#####################################################################
// Copyright 2002, 2003, Ronald Fedkiw, Eran Guendelman, Eilene Hao, Neil Molino, 
// Duc Nguyen, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLES  
//##################################################################### 
//
//#####################################################################
// Fedkiw - June 8, 2003
// Molino - August 19, 2002
// Guendelman - May 8, 2002
// Teran - August 19, 2002
// Nguyen - June 30, 2003
// Hao - August 25, 2003
//#####################################################################
#ifndef __PARTICLES__
#define __PARTICLES__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <iostream>
namespace PhysBAM{

template<class T>
class PARTICLES
{
public:
    int number; // total active particles 
    int array_size; // includes dead space
    int array_buffer_size; // preferred number of extra cells in the buffer
    int smallest_inactive_index; // first dead space where a particle could be added
    ARRAYS<VECTOR<bool,1> > active;
    ARRAYS<VECTOR<T,1> > mass,radius,temperature,density,age;
    ARRAYS<VECTOR<int,1> > id;    // for rendering
    int number_of_attributes; //for writing; default is 7, set to 11 if you need the extra set
protected:
    bool store_mass,store_radius,store_temperature,store_density,store_age,store_id;
    bool store_position,store_velocity,store_acceleration,update_position,update_velocity;

public:
    PARTICLES()
                    :number(0),array_size(0),array_buffer_size(0),smallest_inactive_index(1),number_of_attributes(7),
                     store_mass(false),store_radius(false),store_temperature(false),store_density(false),store_age(false),store_id(false),
                     store_position(false),store_velocity(false),store_acceleration(false),update_position(false),update_velocity(false)
    {}

    virtual ~PARTICLES()
    {}

protected:
    void Clean_Up_Memory()
    {number=0;array_size=0;array_buffer_size=0;smallest_inactive_index=1;active.Resize(1,0);Delete_Mass();Delete_Radius();}

public:
    void Set_Array_Buffer_Size(const int array_buffer_size_input)
    {array_buffer_size=array_buffer_size_input;}

    void Store_Mass()
    {if(!store_mass) mass.Resize(1,active.m);store_mass=true;}

    void Store_Radius()
    {if(!store_radius) radius.Resize(1,active.m);store_radius=true;}

    void Store_Temperature()
    {if(!store_temperature) temperature.Resize(1,active.m);store_temperature=true;}

    void Store_Density()
    {if(!store_density) density.Resize(1,active.m);store_density=true;}

    void Store_Age()
    {if(!store_age) age.Resize(1,active.m);store_age=true;}

    void Store_Id()
    {if(!store_id) id.Resize(1,active.m);store_id=true;}

    void Delete_Mass()
    {mass.Resize(1,0);store_mass=false;}

    void Delete_Radius()
    {radius.Resize(1,0);store_radius=false;}

    void Delete_Temperature()
    {temperature.Resize(1,0);store_temperature=false;}

    void Delete_Density()
    {density.Resize(1,0);store_density=false;}

    void Delete_Age()
    {age.Resize(1,0);store_age=false;}

    void Delete_Id()
    {id.Resize(1,0);store_id=false;}
    
    int Add_Particle()
    {if(smallest_inactive_index > array_size) Increase_Array_Size(1+array_buffer_size);
    int new_index=smallest_inactive_index;active(new_index)=true;number++;
    while(smallest_inactive_index <= array_size && active(smallest_inactive_index)) smallest_inactive_index++; 
    return new_index;}

    void Delete_Particle(const int index)
    {assert(index >= 1 && index <= array_size);active(index)=false;number--;
    if(index < smallest_inactive_index) smallest_inactive_index=index;} // update smallest inactive index

    bool Particles_Compacted() // all the particle are in consecutive order
    {if(smallest_inactive_index == number+1) return true;else return false;}

    void Default()
    {std::cout << "THIS PARTICLES FUNCTION IS NOT DEFINED!" << std::endl;}

protected:
    void Read(std::istream &input_stream)
    {Read_State(input_stream);
    if(array_size > 0){
        if(store_mass) mass.Read(input_stream);
        if(store_radius) radius.Read(input_stream);
        if(store_temperature) temperature.Read(input_stream);
        if(store_density) density.Read(input_stream);
        if(store_age) age.Read(input_stream);
        if(store_id) id.Read(input_stream);}}

    void Write(std::ostream& output_stream) const
    {Write_State(output_stream);
    if(array_size > 0){
        if(store_mass) mass.Write(output_stream);
        if(store_radius) radius.Write(output_stream);
        if(store_temperature) temperature.Write(output_stream);
        if(store_density) density.Write(output_stream);
        if(store_age) age.Write(output_stream);
        if(store_id) id.Write(output_stream);}}

public:
    void Read_State(std::istream& input_stream)
    {input_stream.read((char*)&number,sizeof(int));assert(number >= 0);
    input_stream.read((char*)&array_size,sizeof(int));assert(array_size >= 0);
    input_stream.read((char*)&array_buffer_size,sizeof(int));assert(array_buffer_size >= 0);
    input_stream.read((char*)&smallest_inactive_index,sizeof(int));
    //input_stream.read((char*)&number_of_attributes,sizeof(int));
    input_stream.read((char*)&store_position,sizeof(bool));
    input_stream.read((char*)&store_velocity,sizeof(bool));
    input_stream.read((char*)&store_acceleration,sizeof(bool));
    input_stream.read((char*)&store_mass,sizeof(bool));
    input_stream.read((char*)&store_radius,sizeof(bool));
    input_stream.read((char*)&update_position,sizeof(bool));
    input_stream.read((char*)&update_velocity,sizeof(bool));
    if(number_of_attributes >= 11){
        input_stream.read((char*)&store_temperature,sizeof(bool));
        input_stream.read((char*)&store_density,sizeof(bool));
        input_stream.read((char*)&store_age,sizeof(bool));
        input_stream.read((char*)&store_id,sizeof(bool));}
    else{store_temperature=false;store_density=false;store_age=false;store_id=false;}
    if(array_size > 0) active.Read(input_stream);}

    void Write_State(std::ostream& output_stream) const
    {output_stream.write((const char*)&number,sizeof(int));
    output_stream.write((const char*)&array_size,sizeof(int));
    output_stream.write((const char*)&array_buffer_size,sizeof(int));
    output_stream.write((const char*)&smallest_inactive_index,sizeof(int));
    output_stream.write((const char*)&number_of_attributes,sizeof(int));
    output_stream.write((const char*)&store_position,sizeof(bool));
    output_stream.write((const char*)&store_velocity,sizeof(bool));
    output_stream.write((const char*)&store_acceleration,sizeof(bool));
    output_stream.write((const char*)&store_mass,sizeof(bool));
    output_stream.write((const char*)&store_radius,sizeof(bool));
    output_stream.write((const char*)&update_position,sizeof(bool));
    output_stream.write((const char*)&update_velocity,sizeof(bool));
    if(number_of_attributes >= 11){
        output_stream.write((const char*)&store_temperature,sizeof(bool));
        output_stream.write((const char*)&store_density,sizeof(bool));
        output_stream.write((const char*)&store_age,sizeof(bool));
        output_stream.write((const char*)&store_id,sizeof(bool));}
    if(array_size > 0) active.Write(output_stream);}

//#####################################################################
    virtual void Increase_Array_Size(const int number_of_new_indices){Default();}
//#####################################################################
};      
}
#endif

