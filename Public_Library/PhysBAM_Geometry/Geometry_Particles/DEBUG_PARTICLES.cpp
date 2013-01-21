//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEBUG_PARTICLES<TV>::
DEBUG_PARTICLES()
    :debug_particles(*new GEOMETRY_PARTICLES<TV>),edge_separation(0)
{
    debug_particles.template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    debug_particles.Store_Velocity(true);
    Store_Debug_Particles(this);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEBUG_PARTICLES<TV>::
~DEBUG_PARTICLES()
{
    delete &debug_particles;
}
//#####################################################################
// Function Store_Debug_Particles
//#####################################################################
template<class TV> DEBUG_PARTICLES<TV>* DEBUG_PARTICLES<TV>::
Store_Debug_Particles(DEBUG_PARTICLES<TV>* particle)
{
    static DEBUG_PARTICLES<TV>* stored_particles=0;
    DEBUG_PARTICLES<TV>* tmp=stored_particles;
    if(particle) stored_particles=particle;
    return tmp;
}
//#####################################################################
// Function Write_Debug_Particles
//#####################################################################
template<class TV> void DEBUG_PARTICLES<TV>::
Write_Debug_Particles(STREAM_TYPE stream_type,const std::string& output_directory,int frame) const
{
    FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory.c_str(),frame));
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%i/debug_particles",output_directory.c_str(),frame),debug_particles,debug_objects);
    debug_particles.Delete_All_Elements();
    debug_objects.Remove_All();
}
//#####################################################################
// Function Add_Debug_Particle
//#####################################################################
template<class TV> void PhysBAM::
Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color)
{
    typedef typename TV::SCALAR T;
    DEBUG_PARTICLES<TV>* particles=DEBUG_PARTICLES<TV>::Store_Debug_Particles();
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles->debug_particles.template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    int p=particles->debug_particles.Add_Element();
    particles->debug_particles.X(p)=X;
    (*color_attribute)(p)=color;
}
//#####################################################################
// Function Debug_Particle_Set_Attribute
//#####################################################################
template<class TV,class ATTR> void PhysBAM::
Debug_Particle_Set_Attribute(ATTRIBUTE_ID id,const ATTR& attr)
{
    typedef typename TV::SCALAR T;
    DEBUG_PARTICLES<TV>* particles=DEBUG_PARTICLES<TV>::Store_Debug_Particles();
    ARRAY_VIEW<ATTR>* attribute=particles->debug_particles.template Get_Array<ATTR>(id);
    attribute->Last()=attr;
}
//#####################################################################
// Function Add_Debug_Object
//#####################################################################
template<class TV,int d> void PhysBAM::
Add_Debug_Object(const VECTOR<TV,d>& object,const VECTOR<typename TV::SCALAR,3>& color,const VECTOR<typename TV::SCALAR,3>& bgcolor)
{
    DEBUG_PARTICLES<TV>* dp=DEBUG_PARTICLES<TV>::Store_Debug_Particles();
    DEBUG_OBJECT<TV> obj;
    obj.type=(typename DEBUG_OBJECT<TV>::TYPE)d;
    obj.X=VECTOR<TV,3>(object);
    obj.color=color;
    obj.bgcolor=bgcolor;
    obj.draw_vertices=false;
    obj.separation=dp->edge_separation;
    dp->debug_objects.Append(obj);
}
//#####################################################################
// Function Dump_Surface
//#####################################################################
template<class T_SURFACE,class T> void PhysBAM::
Dump_Surface(const T_SURFACE& surface,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor)
{
    typedef VECTOR<T,T_SURFACE::dimension+1> TV;
    for(int i=0;i<surface.mesh.elements.m;i++)
        Add_Debug_Object(VECTOR<TV,TV::m>(surface.particles.X.Subset(surface.mesh.elements(i))),color,bgcolor);
}
//#####################################################################
// Function Dump_Levelset
//#####################################################################
template<class TV,class TV_INT,class T> void PhysBAM::
Dump_Levelset(const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi,bool node_centered,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor)
{
    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT surface;
    MARCHING_CUBES<TV>::Create_Surface(surface,node_centered?grid:grid.Get_Regular_Grid_At_MAC_Positions(),phi);
    Dump_Surface(surface,color,bgcolor);
}
namespace PhysBAM{
template class DEBUG_PARTICLES<VECTOR<float,1> >;
template class DEBUG_PARTICLES<VECTOR<float,2> >;
template class DEBUG_PARTICLES<VECTOR<float,3> >;
template void Add_Debug_Particle<VECTOR<float,1> >(VECTOR<float,1> const&,VECTOR<float,3> const&);
template void Add_Debug_Particle<VECTOR<float,2> >(VECTOR<float,2> const&,VECTOR<float,3> const&);
template void Add_Debug_Particle<VECTOR<float,3> >(VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,1>,float>(ATTRIBUTE_ID,float const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,2>,float>(ATTRIBUTE_ID,float const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,3>,float>(ATTRIBUTE_ID,float const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,1>,VECTOR<float,1> >(ATTRIBUTE_ID,VECTOR<float,1> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,2>,VECTOR<float,2> >(ATTRIBUTE_ID,VECTOR<float,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,3>,VECTOR<float,3> >(ATTRIBUTE_ID,VECTOR<float,3> const&);
template void Add_Debug_Object<VECTOR<float,3>,2>(const VECTOR<VECTOR<float,3>,2>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,3>,3>(const VECTOR<VECTOR<float,3>,3>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,2>,2>(const VECTOR<VECTOR<float,2>,2>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,2>,3>(const VECTOR<VECTOR<float,2>,3>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template class DEBUG_PARTICLES<VECTOR<double,1> >;
template class DEBUG_PARTICLES<VECTOR<double,2> >;
template class DEBUG_PARTICLES<VECTOR<double,3> >;
template void Add_Debug_Particle<VECTOR<double,1> >(VECTOR<double,1> const&,VECTOR<double,3> const&);
template void Add_Debug_Particle<VECTOR<double,2> >(VECTOR<double,2> const&,VECTOR<double,3> const&);
template void Add_Debug_Particle<VECTOR<double,3> >(VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,1>,double>(ATTRIBUTE_ID,double const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,2>,double>(ATTRIBUTE_ID,double const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,3>,double>(ATTRIBUTE_ID,double const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,1>,VECTOR<double,1> >(ATTRIBUTE_ID,VECTOR<double,1> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,2>,VECTOR<double,2> >(ATTRIBUTE_ID,VECTOR<double,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,3>,VECTOR<double,3> >(ATTRIBUTE_ID,VECTOR<double,3> const&);
template void Add_Debug_Object<VECTOR<double,3>,2>(const VECTOR<VECTOR<double,3>,2>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,3>,3>(const VECTOR<VECTOR<double,3>,3>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,2>,2>(const VECTOR<VECTOR<double,2>,2>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,2>,3>(const VECTOR<VECTOR<double,2>,3>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,1>,2>(const VECTOR<VECTOR<double,1>,2>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,1>,3>(const VECTOR<VECTOR<double,1>,3>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<float,1>,2>(const VECTOR<VECTOR<float,1>,2>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,1>,3>(const VECTOR<VECTOR<float,1>,3>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template void Dump_Levelset<VECTOR<double,1>,VECTOR<int,1>,double>(GRID<VECTOR<double,1> > const&,ARRAY<double,VECTOR<int,1> > const&,bool,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<double,2>,VECTOR<int,2>,double>(GRID<VECTOR<double,2> > const&,ARRAY<double,VECTOR<int,2> > const&,bool,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<double,3>,VECTOR<int,3>,double>(GRID<VECTOR<double,3> > const&,ARRAY<double,VECTOR<int,3> > const&,bool,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<float,1>,VECTOR<int,1>,float>(GRID<VECTOR<float,1> > const&,ARRAY<float,VECTOR<int,1> > const&,bool,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Levelset<VECTOR<float,2>,VECTOR<int,2>,float>(GRID<VECTOR<float,2> > const&,ARRAY<float,VECTOR<int,2> > const&,bool,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Levelset<VECTOR<float,3>,VECTOR<int,3>,float>(GRID<VECTOR<float,3> > const&,ARRAY<float,VECTOR<int,3> > const&,bool,VECTOR<float,3> const&,VECTOR<float,3> const&);
}

