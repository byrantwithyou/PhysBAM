//#####################################################################
// Copyright 2006, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING_DYNAMIC.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/PARTICLE_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
using namespace PhysBAM;
//#####################################################################
// Function Create_Structure
//#####################################################################
template<class TV> BINDING<TV>* BINDING<TV>::
Create(TYPED_ISTREAM& input,DEFORMABLE_PARTICLES<TV>& particles)
{
    try{
        int name;Read_Binary(input,name);
        BINDING* binding=Create_From_Name(name,particles);
        binding->Read_Helper(input);return binding;}
    catch(std::runtime_error& error){
        throw READ_ERROR(error.what());}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void BINDING<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    ARRAY<int> parents;
    Parents(parents);
    for(int i=0;i<parents.m;i++) dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(parents(i),particle_index));
}
//#####################################################################
// Function Create_From_Name
//#####################################################################
template<class TV> void BINDING<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,Name());
    Write_Helper(output);
}
static bool Initialize_Bindings()
{
#define HELPER(T,d) \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<PARTICLE_BINDING<VECTOR<T,d> > >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<LINEAR_BINDING<VECTOR<T,d>,2> >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<LINEAR_BINDING<VECTOR<T,d>,3> >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<LINEAR_BINDING<VECTOR<T,d>,4> >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<LINEAR_BINDING_DYNAMIC<VECTOR<T,d> > >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<RIGID_BODY_BINDING<VECTOR<T,d> > >();

    HELPER(float,1)
    HELPER(float,2)
    HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    HELPER(double,1)
    HELPER(double,2)
    HELPER(double,3)
#endif
    return true;
}
static bool initialize_bindings=Initialize_Bindings();
//#####################################################################
// Function Create_From_Name
//#####################################################################
template<class TV> BINDING<TV>* BINDING<TV>::
Create_From_Name(const int name,DEFORMABLE_PARTICLES<TV>& particles)
{
    PHYSBAM_ASSERT(initialize_bindings);
    BINDING* binding=BINDING_REGISTRY<TV>::Name_To_Factory(name)->Create(dynamic_cast<GEOMETRY_PARTICLES<TV>&>(particles));
    if(!binding){LOG::cerr<<name<<" has no Create(GEOMETRY_PARTICLES<TV>& particles) function."<<std::endl;PHYSBAM_FATAL_ERROR();}
    return binding;
}
//#####################################################################
template class BINDING<VECTOR<float,1> >;
template class BINDING<VECTOR<float,2> >;
template class BINDING<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BINDING<VECTOR<double,1> >;
template class BINDING<VECTOR<double,2> >;
template class BINDING<VECTOR<double,3> >;
#endif
