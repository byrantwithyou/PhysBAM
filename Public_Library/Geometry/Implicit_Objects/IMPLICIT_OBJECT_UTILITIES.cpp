//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{
//#####################################################################
// Function Initialize_Implicit_Surface_Helper
//#####################################################################
template<class TV,class T_SURFACE> static auto
Initialize_Implicit_Surface_Helper(T_SURFACE& surface,int max_res)
{
    typedef typename TV::SCALAR T;
    // undeformed levelset
    LEVELSET_IMPLICIT_OBJECT<TV>& undeformed_levelset=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    surface.Update_Bounding_Box();
    RANGE<TV>& box=*surface.bounding_box;
    GRID<TV>& grid=undeformed_levelset.levelset.grid;
    ARRAY<T,VECTOR<int,TV::m> >& phi=undeformed_levelset.levelset.phi;
    grid=GRID<TV>::Create_Grid_Given_Cell_Size(box,box.Edge_Lengths().Max()/max_res,false,5);
    phi.Resize(grid.Domain_Indices());
    LEVELSET_MAKER_UNIFORM<TV>::Compute_Level_Set(surface,grid,0,phi);
    undeformed_levelset.Update_Box();
    return &undeformed_levelset;
}
//#####################################################################
// Function Initialize_Implicit_Surface
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >*
Initialize_Implicit_Surface(TRIANGULATED_SURFACE<T>& surface,int max_res)
{
    return Initialize_Implicit_Surface_Helper<VECTOR<T,3> >(surface,max_res);
}
//#####################################################################
// Function Initialize_Implicit_Surface
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,2> >*
Initialize_Implicit_Surface(SEGMENTED_CURVE_2D<T>& surface,int max_res)
{
    return Initialize_Implicit_Surface_Helper<VECTOR<T,2> >(surface,max_res);
}
//#####################################################################
// Function Levelset_From_File
//
// TODO: put this function where it belongs
//#####################################################################
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* 
Levelset_From_Tri_File(const std::string& filename,int max_resolution)
{
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    Read_From_File(STREAM_TYPE(0.f),filename,*surface);
    LOG::printf("Read mesh: %d triangle, %d particles\n",surface->mesh.elements.m,surface->particles.number);
    surface->mesh.Initialize_Adjacent_Elements();
    surface->mesh.Initialize_Neighbor_Nodes();
    surface->mesh.Initialize_Incident_Elements();
    surface->Update_Bounding_Box();
    surface->Initialize_Hierarchy();
    surface->Update_Triangle_List();
    return Initialize_Implicit_Surface(*surface,max_resolution);
}
template LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> >*
Initialize_Implicit_Surface<double>(TRIANGULATED_SURFACE<double>&,int);
template LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> >*
Initialize_Implicit_Surface<float>(TRIANGULATED_SURFACE<float>&,int);
template LEVELSET_IMPLICIT_OBJECT<VECTOR<double,2> >*
Initialize_Implicit_Surface<double>(SEGMENTED_CURVE_2D<double>&,int);
template LEVELSET_IMPLICIT_OBJECT<VECTOR<float,2> >*
Initialize_Implicit_Surface<float>(SEGMENTED_CURVE_2D<float>&,int);
template LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> >*
Levelset_From_Tri_File<double>(std::string const&,int);
template LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> >*
Levelset_From_Tri_File<float>(std::string const&,int);
}
