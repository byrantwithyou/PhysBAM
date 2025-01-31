//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEBUG_PARTICLES<TV>::
DEBUG_PARTICLES()
    :debug_particles(*new GEOMETRY_PARTICLES<TV>),edge_separation(0)
{
    debug_particles.template Add_Array<VECTOR<T,3> >("color");
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
// Function Clear_Debug_Particles
//#####################################################################
template<class TV> void DEBUG_PARTICLES<TV>::
Clear_Debug_Particles() const
{
    debug_particles.Delete_All_Elements();
    debug_objects.Remove_All();
    debug_text.Remove_All();
}
//#####################################################################
// Function Write_Debug_Particles
//#####################################################################
template<class TV> void DEBUG_PARTICLES<TV>::
Write_Debug_Particles(STREAM_TYPE stream_type,const VIEWER_DIR& viewer_dir) const
{
    Write_To_File(stream_type,viewer_dir.current_directory+"/debug_particles",
        debug_particles,debug_objects,debug_text);
    Clear_Debug_Particles();
}
//#####################################################################
// Function Add_Debug_Particle
//#####################################################################
template<class TV> void
Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color)
{
    typedef typename TV::SCALAR T;
    DEBUG_PARTICLES<TV>* particles=DEBUG_PARTICLES<TV>::Store_Debug_Particles();
    ARRAY_VIEW<VECTOR<T,3> >* color_attribute=particles->debug_particles.template Get_Array<VECTOR<T,3> >("color");
    int p=particles->debug_particles.Add_Element();
    particles->debug_particles.X(p)=X;
    (*color_attribute)(p)=color;
}
//#####################################################################
// Function Debug_Particle_Set_Attribute
//#####################################################################
template<class TV,class ATTR> void
Debug_Particle_Set_Attribute(const std::string& name,const ATTR& attr)
{
    DEBUG_PARTICLES<TV>* particles=DEBUG_PARTICLES<TV>::Store_Debug_Particles();
    ARRAY_VIEW<ATTR>* attribute=particles->debug_particles.template Get_Array<ATTR>(name);
    attribute->Last()=attr;
}
//#####################################################################
// Function Add_Debug_Object
//#####################################################################
template<class TV,int d> void
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
template<class T_SURFACE,class T> void
Dump_Surface(const T_SURFACE& surface,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor)
{
    typedef VECTOR<T,T_SURFACE::dimension+1> TV;
    for(int i=0;i<surface.mesh.elements.m;i++)
        Add_Debug_Object(VECTOR<TV,TV::m>(surface.particles.X.Subset(surface.mesh.elements(i))),color,bgcolor);
}
//#####################################################################
// Function Dump_Levelset
//#####################################################################
template<class TV,class TV_INT,class T> void
Dump_Levelset(const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor)
{
    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT surface;
    MARCHING_CUBES<TV>::Create_Surface(surface,grid,phi);
    Dump_Surface(surface,color,bgcolor);
}
//#####################################################################
// Function Dump_Levelset
//#####################################################################
template<class TV,class TV_INT,class T> void
Dump_Levelset(const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor,T contour_value)
{
    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT surface;
    MARCHING_CUBES<TV>::Create_Surface(surface,grid,phi,contour_value);
    Dump_Surface(surface,color,bgcolor);
}
//#####################################################################
// Function Dump_Levelset
//#####################################################################
template<class TV,class T> void
Dump_Levelset(const GRID<TV>& grid,const IMPLICIT_OBJECT<TV>& phi,const VECTOR<T,3>& color,const VECTOR<T,3>& bgcolor)
{
    GRID<TV> node_grid(grid.Is_MAC_Grid()?grid.Get_Regular_Grid():grid);
    ARRAY<T,VECTOR<int,TV::m> > phi_array(node_grid.Numbers_Of_Nodes());
    for(NODE_ITERATOR<TV> it(node_grid);it.Valid();it.Next())
        phi_array(it.index)=phi.Extended_Phi(it.Location());
    Dump_Levelset(node_grid,phi_array,color,bgcolor);
}
//#####################################################################
// Function Add_Debug_Text
//#####################################################################
template<class TV> inline void
Add_Debug_Text(const TV& X,const std::string& text,const VECTOR<typename TV::SCALAR,3>& color)
{
    DEBUG_PARTICLES<TV>* dp=DEBUG_PARTICLES<TV>::Store_Debug_Particles();
    DEBUG_TEXT<TV> dt;
    dt.X=X;
    dt.text=text;
    dt.color=color;
    dp->debug_text.Append(dt);
}
//#####################################################################
// Function Get_Debug_Particles
//#####################################################################
template<class TV> DEBUG_PARTICLES<TV>&
Get_Debug_Particles()
{
    static DEBUG_PARTICLES<TV> debug_particles;
    return debug_particles;
}
template class DEBUG_PARTICLES<VECTOR<float,1> >;
template class DEBUG_PARTICLES<VECTOR<float,2> >;
template class DEBUG_PARTICLES<VECTOR<float,3> >;
template class DEBUG_PARTICLES<VECTOR<double,1> >;
template class DEBUG_PARTICLES<VECTOR<double,2> >;
template class DEBUG_PARTICLES<VECTOR<double,3> >;

template void Add_Debug_Object<VECTOR<double,1>,2>(const VECTOR<VECTOR<double,1>,2>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,1>,3>(const VECTOR<VECTOR<double,1>,3>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,2>,2>(const VECTOR<VECTOR<double,2>,2>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,2>,3>(const VECTOR<VECTOR<double,2>,3>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,3>,2>(const VECTOR<VECTOR<double,3>,2>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<double,3>,3>(const VECTOR<VECTOR<double,3>,3>&,const VECTOR<double,3>&,const VECTOR<double,3>&);
template void Add_Debug_Object<VECTOR<float,1>,2>(const VECTOR<VECTOR<float,1>,2>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,1>,3>(const VECTOR<VECTOR<float,1>,3>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,2>,2>(const VECTOR<VECTOR<float,2>,2>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,2>,3>(const VECTOR<VECTOR<float,2>,3>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,3>,2>(const VECTOR<VECTOR<float,3>,2>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template void Add_Debug_Object<VECTOR<float,3>,3>(const VECTOR<VECTOR<float,3>,3>&,const VECTOR<float,3>&,const VECTOR<float,3>&);
template void Add_Debug_Particle<VECTOR<double,1> >(VECTOR<double,1> const&,VECTOR<double,3> const&);
template void Add_Debug_Particle<VECTOR<double,2> >(VECTOR<double,2> const&,VECTOR<double,3> const&);
template void Add_Debug_Particle<VECTOR<double,3> >(VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Add_Debug_Particle<VECTOR<float,1> >(VECTOR<float,1> const&,VECTOR<float,3> const&);
template void Add_Debug_Particle<VECTOR<float,2> >(VECTOR<float,2> const&,VECTOR<float,3> const&);
template void Add_Debug_Particle<VECTOR<float,3> >(VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,1>,VECTOR<double,1> >(const std::string&,VECTOR<double,1> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,1>,double>(const std::string&,double const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,2>,VECTOR<double,2> >(const std::string&,VECTOR<double,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,2>,double>(const std::string&,double const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,3>,VECTOR<double,3> >(const std::string&,VECTOR<double,3> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,3>,double>(const std::string&,double const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,1>,VECTOR<float,1> >(const std::string&,VECTOR<float,1> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,1>,float>(const std::string&,float const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,2>,VECTOR<float,2> >(const std::string&,VECTOR<float,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,2>,float>(const std::string&,float const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,3>,VECTOR<float,3> >(const std::string&,VECTOR<float,3> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,3>,float>(const std::string&,float const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,1>,VECTOR<int,2> >(const std::string&,VECTOR<int,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,2>,VECTOR<int,2> >(const std::string&,VECTOR<int,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<double,3>,VECTOR<int,2> >(const std::string&,VECTOR<int,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,1>,VECTOR<int,2> >(const std::string&,VECTOR<int,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,2>,VECTOR<int,2> >(const std::string&,VECTOR<int,2> const&);
template void Debug_Particle_Set_Attribute<VECTOR<float,3>,VECTOR<int,2> >(const std::string&,VECTOR<int,2> const&);
template void Dump_Levelset<VECTOR<double,1>,VECTOR<int,1>,double>(GRID<VECTOR<double,1> > const&,ARRAY<double,VECTOR<int,1> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<double,2>,VECTOR<int,2>,double>(GRID<VECTOR<double,2> > const&,ARRAY<double,VECTOR<int,2> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<double,3>,VECTOR<int,3>,double>(GRID<VECTOR<double,3> > const&,ARRAY<double,VECTOR<int,3> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<float,1>,VECTOR<int,1>,float>(GRID<VECTOR<float,1> > const&,ARRAY<float,VECTOR<int,1> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Levelset<VECTOR<float,2>,VECTOR<int,2>,float>(GRID<VECTOR<float,2> > const&,ARRAY<float,VECTOR<int,2> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Levelset<VECTOR<float,3>,VECTOR<int,3>,float>(GRID<VECTOR<float,3> > const&,ARRAY<float,VECTOR<int,3> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Levelset<VECTOR<float,1>,VECTOR<int,1>,float>(GRID<VECTOR<float,1> > const&,ARRAY<float,VECTOR<int,1> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&,float);
template void Dump_Levelset<VECTOR<float,2>,VECTOR<int,2>,float>(GRID<VECTOR<float,2> > const&,ARRAY<float,VECTOR<int,2> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&,float);
template void Dump_Levelset<VECTOR<float,3>,VECTOR<int,3>,float>(GRID<VECTOR<float,3> > const&,ARRAY<float,VECTOR<int,3> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&,float);
template void Dump_Levelset<VECTOR<double,1>,VECTOR<int,1>,double>(GRID<VECTOR<double,1> > const&,ARRAY<double,VECTOR<int,1> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&,double);
template void Dump_Levelset<VECTOR<double,2>,VECTOR<int,2>,double>(GRID<VECTOR<double,2> > const&,ARRAY<double,VECTOR<int,2> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&,double);
template void Dump_Levelset<VECTOR<double,3>,VECTOR<int,3>,double>(GRID<VECTOR<double,3> > const&,ARRAY<double,VECTOR<int,3> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&,double);
template void Dump_Levelset<VECTOR<double,1>,double>(GRID<VECTOR<double,1> > const&,IMPLICIT_OBJECT<VECTOR<double,1> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<double,2>,double>(GRID<VECTOR<double,2> > const&,IMPLICIT_OBJECT<VECTOR<double,2> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<double,3>,double>(GRID<VECTOR<double,3> > const&,IMPLICIT_OBJECT<VECTOR<double,3> > const&,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Levelset<VECTOR<float,1>,float>(GRID<VECTOR<float,1> > const&,IMPLICIT_OBJECT<VECTOR<float,1> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Levelset<VECTOR<float,2>,float>(GRID<VECTOR<float,2> > const&,IMPLICIT_OBJECT<VECTOR<float,2> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Levelset<VECTOR<float,3>,float>(GRID<VECTOR<float,3> > const&,IMPLICIT_OBJECT<VECTOR<float,3> > const&,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Surface<TRIANGULATED_SURFACE<float>,float>(TRIANGULATED_SURFACE<float> const&,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Surface<TRIANGULATED_SURFACE<double>,double>(TRIANGULATED_SURFACE<double> const&,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Dump_Surface<SEGMENTED_CURVE_2D<float>,float>(SEGMENTED_CURVE_2D<float> const&,VECTOR<float,3> const&,VECTOR<float,3> const&);
template void Dump_Surface<SEGMENTED_CURVE_2D<double>,double>(SEGMENTED_CURVE_2D<double> const&,VECTOR<double,3> const&,VECTOR<double,3> const&);
template void Add_Debug_Text<VECTOR<double,2> >(VECTOR<double,2> const&,std::string const&,VECTOR<double,3> const&);
template void Add_Debug_Text<VECTOR<float,2> >(VECTOR<float,2> const&,std::string const&,VECTOR<float,3> const&);
template void Add_Debug_Text<VECTOR<double,3> >(VECTOR<double,3> const&,std::string const&,VECTOR<double,3> const&);
template void Add_Debug_Text<VECTOR<float,3> >(VECTOR<float,3> const&,std::string const&,VECTOR<float,3> const&);
template DEBUG_PARTICLES<VECTOR<float,1> >& Get_Debug_Particles<VECTOR<float,1> >();
template DEBUG_PARTICLES<VECTOR<float,2> >& Get_Debug_Particles<VECTOR<float,2> >();
template DEBUG_PARTICLES<VECTOR<float,3> >& Get_Debug_Particles<VECTOR<float,3> >();
template DEBUG_PARTICLES<VECTOR<double,1> >& Get_Debug_Particles<VECTOR<double,1> >();
template DEBUG_PARTICLES<VECTOR<double,2> >& Get_Debug_Particles<VECTOR<double,2> >();
template DEBUG_PARTICLES<VECTOR<double,3> >& Get_Debug_Particles<VECTOR<double,3> >();
}
