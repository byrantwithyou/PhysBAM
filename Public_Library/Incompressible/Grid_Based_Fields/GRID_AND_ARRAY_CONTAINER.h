//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_AND_ARRAY_CONTAINER
//#####################################################################
#ifndef __GRID_AND_ARRAY_CONTAINER__    
#define __GRID_AND_ARRAY_CONTAINER__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Vectors/VECTOR_UTILITIES.h>
#include <Grid_PDE/Advection/ADVECTION_FORWARD.h>
#include <Grid_PDE/Advection/ADVECTION_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class TV> struct BOUNDARY_POLICY;
template<class TV,class T> struct BOUNDARY;
template<class TV> struct GRID_ARRAYS_POLICY;

template<class TV,class T2>
class GRID_AND_ARRAY_CONTAINER
{
    typedef typename TV::SCALAR T;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,TV::m> TV_SIDES;typedef VECTOR<int,TV::m> TV_INT;
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T> T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
public:
    GRID<TV>& grid;
    ARRAY<T2,TV_INT> array;
    ADVECTION<TV,T>* advection;
    ADVECTION<TV,T>* advection_maccormack;
    BOUNDARY<TV,T>* boundary;
private:
    T_ADVECTION_SEMI_LAGRANGIAN_SCALAR& advection_default;
protected:
    BOUNDARY<TV,T>& boundary_default; 
    const ARRAY<T,FACE_INDEX<TV::m> >* face_velocities;
    const ARRAY<TV,TV_INT>* cell_velocities;
public:

    GRID_AND_ARRAY_CONTAINER(GRID<TV>& grid_input);
    GRID_AND_ARRAY_CONTAINER(const GRID_AND_ARRAY_CONTAINER&) = delete;
    void operator=(const GRID_AND_ARRAY_CONTAINER&) = delete;
    virtual ~GRID_AND_ARRAY_CONTAINER();

    void Clean_Memory()
    {array.Clean_Memory();}
    
    virtual void Initialize_Array(const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {array.Resize(grid.Cell_Indices(ghost_cells),initialize_new_elements,copy_existing_elements);}
  
    void Initialize_Domain_Boundary_Conditions(const TV_SIDES& domain_walls=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)))
    {boundary->Set_Constant_Extrapolation(Complement(domain_walls));}

    void Set_Custom_Boundary(BOUNDARY<TV,T>& boundary_input)
    {boundary=&boundary_input;}

    void Set_Custom_Advection(ADVECTION<TV,T>& advection_input)
    {advection=&advection_input;}

    void Set_To_Constant_Value(const T2& value)
    {array.Fill(value);}

    void Set_Velocity(const ARRAY<T,FACE_INDEX<TV::m> >* face_velocities_input)
    {face_velocities=face_velocities_input;}
    
    void Set_Velocity(const ARRAY<TV,TV_INT>* cell_velocities_input)
    {cell_velocities=cell_velocities_input;}

    template<class T_ARRAYS_BOOL>
    void Use_Maccormack_Advection(const T_ARRAYS_BOOL& cell_mask)
    {advection_maccormack=new ADVECTION_MACCORMACK_UNIFORM<TV,T,ADVECTION<TV,T> >(*advection,0,&cell_mask,0);
    Set_Custom_Advection(*advection_maccormack);}

//#####################################################################
    virtual void Euler_Step(const T dt,const T time,const int number_of_ghost_cells);
//#####################################################################
};
}
#endif
