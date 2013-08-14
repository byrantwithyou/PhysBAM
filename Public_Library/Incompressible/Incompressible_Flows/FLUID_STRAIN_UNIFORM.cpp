//#####################################################################
// Copyright 2004-2007, Doug Enright, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andy Lutomirski, Paul-James White.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Tools/Matrices/MATRIX_1X1.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/ROTATION.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Incompressible/Incompressible_Flows/EXTERNAL_STRAIN_ADJUSTMENT.h>
#include <Incompressible/Incompressible_Flows/FLUID_STRAIN_UNIFORM.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Incompressible/Level_Sets/LEVELSET_MULTIPLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_STRAIN_UNIFORM<TV>::
FLUID_STRAIN_UNIFORM(const GRID<TV>& grid_input)
    :external_strain_adjustment(0),e_boundary_default(*new BOUNDARY<TV,SYMMETRIC_MATRIX<T,TV::m> >),e_advection_default(*new T_ADVECTION_SEMI_LAGRANGIAN_SYMMETRIC_MATRIX),cfl_called(false)
{
    e_advection=&e_advection_default;
    e_boundary=&e_boundary_default;
    Initialize_Grid(grid_input);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> FLUID_STRAIN_UNIFORM<TV>::
~FLUID_STRAIN_UNIFORM()
{
    delete &e_boundary_default;delete &e_advection_default;
}
//#####################################################################
// Function Update_Strain_Equation_Helper_Cell_Centered
//#####################################################################
template<class TV> void FLUID_STRAIN_UNIFORM<TV>::
Update_Strain_Equation_Helper_Cell_Centered(const T dt,const T time,const T density,const T heaviside_bandwidth,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<TV,TV_INT>& V,
    const ARRAY<T,TV_INT>& phi_ghost,const int number_of_ghost_cells)
{
    T_ARRAYS_SYMMETRIC_MATRIX e_ghost(grid.Domain_Indices(number_of_ghost_cells),false);e_boundary->Fill_Ghost_Cells(grid,e,e_ghost,dt,time,number_of_ghost_cells);

    // update the strain to time n+1
    e_advection->Update_Advection_Equation_Cell(grid,e,e_ghost,face_velocities_ghost,*e_boundary,dt,time);

    TV one_over_two_DX=(T).5*Inverse(grid.dX);T one_over_dimension=1/TV::m;
    for(CELL_ITERATOR<TV> iterator(grid,0);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();MATRIX<T,TV::dimension> VX;
        for(int axis=0;axis<TV::m;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
            VX.Set_Column(axis,one_over_two_DX[axis]*(V(cell+offset)-V(cell-offset)));} // TODO: maybe change this to do a better difference for the derivative in the same axis as the face
        e(cell)=SYMMETRIC_MATRIX<T,TV::m>::Conjugate(MATRIX<T,TV::dimension>::Rotation_Matrix(dt*VX.Antisymmetric_Part_Cross_Product_Vector()),e(cell));
        e(cell)+=dt*VX.Symmetric_Part();
        if(plasticity_alpha){
            SYMMETRIC_MATRIX<T,TV::m> e_prime=e_ghost(cell)-one_over_dimension*e_ghost(cell).Trace();T e_prime_norm=e_prime.Frobenius_Norm();
            if(e_prime_norm>plasticity_gamma) e(cell)-=dt*plasticity_alpha*(e_prime_norm-plasticity_gamma)/e_prime_norm*e_prime;}}
    e_boundary->Apply_Boundary_Condition(grid,e,time+dt);

    // surface boundary condition
    e_boundary->Fill_Ghost_Cells(grid,e,e_ghost,dt,time+dt,number_of_ghost_cells);
    if(external_strain_adjustment) external_strain_adjustment->Adjust_Strain(e_ghost,time);
    T epsilon=heaviside_bandwidth*grid.dX.Max();
    for(CELL_ITERATOR<TV> iterator(grid,2);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        e_ghost(cell)*=1-LEVELSET_UTILITIES<T>::Heaviside(phi_ghost(cell),0*epsilon);} // surface boundary condition

    // update velocity due to strain at time n+1
    TV dt_elastic_modulus_over_density_over_two_DX=dt*elastic_modulus/density*one_over_two_DX;
    V.Fill(TV());
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        for(int axis=0;axis<TV::m;axis++){TV_INT offset=TV_INT::Axis_Vector(axis);
            V(cell)+=dt_elastic_modulus_over_density_over_two_DX[axis]*(e_ghost(cell+offset).Column(axis)-e_ghost(cell-offset).Column(axis));}}
}
//#####################################################################
// Function Update_Strain_Equation
//#####################################################################
template<class TV> void FLUID_STRAIN_UNIFORM<TV>::
Update_Strain_Equation(const T dt,const T time,const T density,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const ARRAY<T,TV_INT>& phi_ghost,const int number_of_ghost_cells)
{
    if(!cfl_called) PHYSBAM_WARNING("Using strain without calling strain CFL");
    ARRAY<TV,TV_INT> V(grid.Domain_Indices(1));
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next())for(int axis=0;axis<TV::m;axis++)
        V(iterator.Cell_Index())[axis]=(T).5*(face_velocities_ghost.Component(axis)(iterator.First_Face_Index(axis))+face_velocities_ghost.Component(axis)(iterator.Second_Face_Index(axis)));
    Update_Strain_Equation_Helper_Cell_Centered(dt,time,density,(T)1.5,face_velocities_ghost,V,phi_ghost,number_of_ghost_cells);
    // apply the velocity update to the face velocities
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis(); // TODO: use a better difference as above
        face_velocities.Component(axis)(iterator.Face_Index())+=(T).5*(V(iterator.First_Cell_Index())[axis]+V(iterator.Second_Cell_Index())[axis]);}
}
//#####################################################################
// Function Update_Strain_Equation
//#####################################################################
template<class TV> void FLUID_STRAIN_UNIFORM<TV>::
Update_Strain_Equation_Multiphase(const T dt,const T time,const T density,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,
    const LEVELSET_MULTIPLE<TV>& levelset,const int region,const int number_of_ghost_cells)
{
    if(!cfl_called) PHYSBAM_WARNING("Using strain without calling strain CFL");
    ARRAY<TV,TV_INT> V(grid.Domain_Indices(1));
    for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next())for(int axis=0;axis<TV::m;axis++)
        V(iterator.Cell_Index())[axis]=(T).5*(face_velocities_ghost.Component(axis)(iterator.First_Face_Index(axis))+face_velocities_ghost.Component(axis)(iterator.Second_Face_Index(axis)));
    Update_Strain_Equation_Helper_Cell_Centered(dt,time,density,0,face_velocities_ghost,V,levelset.phis(region),number_of_ghost_cells);
    // apply the velocity update to the face velocities
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        if(levelset.Inside_Region_Face(axis,iterator.Face_Index())==region)
            face_velocities.Component(axis)(iterator.Face_Index())+=(T).5*(V(iterator.First_Cell_Index())[axis]+V(iterator.Second_Cell_Index())[axis]);}// TODO: use a better difference as above
}
//#####################################################################
// Function Extrapolate_Strain_Across_Interface
//#####################################################################
template<class TV> void FLUID_STRAIN_UNIFORM<TV>::
Extrapolate_Strain_Across_Interface(ARRAY<T,TV_INT>& phi_ghost,const T band_width)
{
    T delta=band_width*grid.dX.Max();
    for(CELL_ITERATOR<TV> iterator(grid,0);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(phi_ghost(cell)>=delta) e(cell)=SYMMETRIC_MATRIX<T,TV::m>();}
    EXTRAPOLATION_UNIFORM<TV,SYMMETRIC_MATRIX<T,TV::m> > extrapolate(grid,phi_ghost,e,3);extrapolate.Set_Band_Width(band_width);extrapolate.Extrapolate();
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR FLUID_STRAIN_UNIFORM<TV>::
CFL(const T density) const
{
    cfl_called=true;
    return elastic_modulus?grid.dX.Min()*sqrt(density/elastic_modulus):(T)FLT_MAX;
}
//#####################################################################
namespace PhysBAM{
template class FLUID_STRAIN_UNIFORM<VECTOR<float,1> >;
template class FLUID_STRAIN_UNIFORM<VECTOR<float,2> >;
template class FLUID_STRAIN_UNIFORM<VECTOR<float,3> >;
template class FLUID_STRAIN_UNIFORM<VECTOR<double,1> >;
template class FLUID_STRAIN_UNIFORM<VECTOR<double,2> >;
template class FLUID_STRAIN_UNIFORM<VECTOR<double,3> >;
}
