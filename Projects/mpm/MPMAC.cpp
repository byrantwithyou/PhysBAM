//#####################################################################
// Copyright 2013
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include "MPMAC_POISSON_SYSTEM.h"
#include "MPMAC_POISSON_VECTOR.h"
#include "MPMAC.h"
namespace PhysBAM{

//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MPMAC<TV>::
Initialize()
{
    mac_grid.Initialize(grid.numbers_of_cells,RANGE<TV>(grid.domain.min_corner,grid.domain.max_corner),true);
    for(int d=0;d<TV::m;d++) face_grid[d]=mac_grid.Get_Face_Grid(d);
    face_velocities.Resize(mac_grid);
    face_velocities_old.Resize(mac_grid);
    face_masses.Resize(mac_grid);
    face_momenta.Resize(mac_grid);
    cell_dirichlet.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    cell_neumann.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    neumann_cell_normal_axis.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    div_u.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    pressure.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    for(int d=0;d<TV::m;d++){
        influence_corner[d].Resize(particles.number);
        weight[d].Resize(particles.number);
        grad_weight[d].Resize(particles.number);
        for(int p=0;p<particles.number;p++){
            weight[d](p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));
            grad_weight[d](p).Resize(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));}}
    min_mass=particles.mass.Min()*(T)1e-5;
    frame=0;
}

//#####################################################################
// Function Reinitialize
//#####################################################################
template<class TV> void MPMAC<TV>::
Reinitialize()
{
    face_velocities.Fill((T)0);
    face_velocities_old.Fill((T)0);
    face_masses.Fill((T)0);
    face_momenta.Fill((T)0);
    cell_dirichlet.Fill(true);
    cell_neumann.Fill(false);
    div_u.Fill((T)0);
    pressure.Fill((T)0);
    max_div=(T)0;
}

//#####################################################################
// Function Weights
//#####################################################################
template<class TV> void MPMAC<TV>::
Weights()
{
#pragma omp parallel for
    for(int p=0;p<particles.number;p++)
        for(int axis=0;axis<TV::m;axis++)
            grid_basis_function.Build_Weights_And_Grad_Weights_Exact(particles.X(p),face_grid[axis],influence_corner[axis](p),weight[axis](p),grad_weight[axis](p));
}

//#####################################################################
// Function Rasterize
//#####################################################################
template<class TV> void MPMAC<TV>::
Rasterize()
{
    for(int p=0;p<particles.number;p++){
        for(int axis=0;axis<TV::m;axis++){
            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
                TV_INT ind=influence_corner[axis](p)+it.index;
                FACE_INDEX<TV::m> face_ind(axis,ind);
                face_masses(face_ind)+=weight[axis](p)(it.index)*particles.mass(p);
                face_momenta(face_ind)+=weight[axis](p)(it.index)*particles.mass(p)*particles.V(p)(axis);}}}
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face_ind=iterator.Full_Index();
        if(face_masses(face_ind)>min_mass){
            face_velocities(face_ind)=face_momenta(face_ind)/face_masses(face_ind);}}
}

//#####################################################################
// Function Advection
//#####################################################################
template<class TV> void MPMAC<TV>::
Advection()
{
    face_velocities_old=face_velocities;
    if(use_gravity){
        for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),face_grid[1].counts));it.Valid();it.Next()){
            FACE_INDEX<TV::m> face_ind(1,it.index);
            if(face_masses(face_ind)>min_mass){
                face_velocities(face_ind)+=dt*(-9.8);
            }
        }
    }
}

//#####################################################################
// Function Identify_Dirichlet
//#####################################################################
template<class TV> void MPMAC<TV>::
Identify_Dirichlet()
{
    for(int p=0;p<particles.number;p++){
        TV_INT cell=mac_grid.Cell(particles.X(p),0);
        cell_dirichlet(cell)=false;}
    
//    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),mac_grid.counts));it.Valid();it.Next()){
//        T total_mass=(T)0;
//        for(int d=0;d<TV::m;d++){
//            total_mass+=face_masses(FACE_INDEX<TV::m>(d,mac_grid.First_Face_Index_In_Cell(d,it.index)));
//            total_mass+=face_masses(FACE_INDEX<TV::m>(d,mac_grid.Second_Face_Index_In_Cell(d,it.index)));}
//        if(total_mass>min_mass*1e2)
//            cell_dirichlet(it.index)=false;}
}

//#####################################################################
// Function Identify_Neumann
//#####################################################################
template<class TV> void MPMAC<TV>::
Identify_Neumann()
{
    cell_neumann.Fill(false);
    // Do it outside
}

//#####################################################################
// Function Build_Velocity_Divergence
//#####################################################################
template<class TV> void MPMAC<TV>::
Build_Velocity_Divergence()
{
    max_div=(T)0;
    div_u.Fill((T)0);
    T one_over_h=(T)1/mac_grid.dX.Min();
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),mac_grid.counts));it.Valid();it.Next()){
        if(!cell_dirichlet(it.index)){
            for(int axis=0;axis<TV::m;axis++)
                div_u(it.index)+=face_velocities(FACE_INDEX<TV::m>(axis,mac_grid.Second_Face_Index_In_Cell(axis,it.index)))
                    -face_velocities(FACE_INDEX<TV::m>(axis,mac_grid.First_Face_Index_In_Cell(axis,it.index)));
            div_u(it.index)*=one_over_h;}
        else div_u(it.index)=(T)0; // dirichlet p cells
        if(!cell_neumann(it.index) && abs(div_u(it.index))>max_div) // record maximum divergence on non neumann cells
            max_div=abs(div_u(it.index));}
}

//#####################################################################
// Function Fix_RHS_Neumann_Cells
//#####################################################################
template<class TV> void MPMAC<TV>::
Fix_RHS_Neumann_Cells(ARRAY<T,TV_INT>& rhs)
{
    T one_over_h=(T)1/mac_grid.dX.Min();
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+mac_grid.counts));it.Valid();it.Next()){
        if(!cell_dirichlet(it.index) && !cell_neumann(it.index)){ // cell is fluid
            for(int d=0;d<TV::m;d++){
                TV_INT left_index=it.index;left_index(d)--;
                TV_INT right_index=it.index;right_index(d)++;
                if(cell_neumann(left_index)){
                        rhs(it.index)-=one_over_h*face_velocities(FACE_INDEX<TV::m>(d,mac_grid.First_Face_Index_In_Cell(d,it.index)));}
                if(cell_neumann(right_index)){
                        rhs(it.index)+=one_over_h*face_velocities(FACE_INDEX<TV::m>(d,mac_grid.Second_Face_Index_In_Cell(d,it.index)));}}}}
}

//#####################################################################
// Function Solve_For_Pressure
// The equation is -div(dt/rho*grad(p))=-div(u)
//#####################################################################
template<class TV> void MPMAC<TV>::
Solve_For_Pressure()
{
    MPMAC_POISSON_SYSTEM<TV> system(debug_cast<MPMAC<TV>&>(*this));
    MPMAC_POISSON_VECTOR<TV> rhs,x;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> vectors;
    rhs.v.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    x.v.Resize(RANGE<TV_INT>(TV_INT(),mac_grid.counts));
    KRYLOV_SOLVER<T>::Ensure_Size(vectors,x,3);
    rhs.v.Copy(-(T)1.0,div_u);
    Fix_RHS_Neumann_Cells(rhs.v);
    x.v=rhs.v;
    system.Test_System(*vectors(0),*vectors(1),*vectors(2));
    CONJUGATE_GRADIENT<T> cg;
    CONJUGATE_RESIDUAL<T> cr;
    KRYLOV_SOLVER<T>* solver=&cg;
    solver->print_residuals=false;
    solver->Solve(system,x,rhs,vectors,(T)1e-7,0,1000);
    pressure=x.v;
    vectors.Delete_Pointers_And_Clean_Memory();
    
    T max_pressure_on_dirichlet=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+mac_grid.counts));it.Valid();it.Next()){
        if(cell_dirichlet(it.index))
            if(sqr(pressure(it.index)>max_pressure_on_dirichlet))
                max_pressure_on_dirichlet=sqr(pressure(it.index));}
    LOG::cout<<"max_pressure_on_dirichlet "<<max_pressure_on_dirichlet<<std::endl;
    
}

//#####################################################################
// Function Do_Projection
// The equation is u=u-dt/rho*grad(p)
//#####################################################################
template<class TV> void MPMAC<TV>::
Do_Projection()
{
    LOG::cout<<"Maximum velocity divergence before projection: "<<max_div<<std::endl;        
    T one_over_h=(T)1/mac_grid.dX.Min();

    if(!uniform_density){
        for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
            FACE_INDEX<TV::m> face_index=iterator.Full_Index();
            int axis=iterator.Axis();
            TV_INT first_cell=iterator.First_Cell_Index();
            TV_INT second_cell=iterator.Second_Cell_Index();
            if(first_cell(axis)>=0&&second_cell(axis)<mac_grid.counts(axis)){ // only deal with non-boundary faces
                if(cell_dirichlet(first_cell)){
                    if(sqr(pressure(first_cell))>1e-4){
                        LOG::cout<<pressure(first_cell)<<std::endl;
                        PHYSBAM_FATAL_ERROR();}}
                if(cell_dirichlet(second_cell)){
                    if(sqr(pressure(second_cell))>1e-4){
                        LOG::cout<<pressure(second_cell)<<std::endl;
                        PHYSBAM_FATAL_ERROR();}}
                if(!cell_neumann(first_cell) && !cell_neumann(second_cell)){
                    T grad_p=(pressure(second_cell)-pressure(first_cell))*one_over_h;
                    if(face_masses(face_index)>min_mass) face_velocities(face_index)-=dt/face_masses(face_index)*grad_p;}}}}
    else{
        for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
            FACE_INDEX<TV::m> face_index=iterator.Full_Index();
            int axis=iterator.Axis();
            TV_INT first_cell=iterator.First_Cell_Index();
            TV_INT second_cell=iterator.Second_Cell_Index();        
            if(first_cell(axis)>=0&&second_cell(axis)<mac_grid.counts(axis)){ // only deal with non-boundary faces
                if(!cell_neumann(first_cell) && !cell_neumann(second_cell)){
                    T grad_p=(pressure(second_cell)-pressure(first_cell))*one_over_h;
                    if(face_masses(face_index)>min_mass) face_velocities(face_index)-=dt*grad_p;}}}}

    // Enforce face velocities for neumann cells
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),mac_grid.counts));it.Valid();it.Next()){
        if(cell_neumann(it.index)){
            int d=neumann_cell_normal_axis(it.index);
            int axis=abs(d)-1;
            // if(face_velocities(FACE_INDEX<TV::m>(axis,mac_grid.First_Face_Index_In_Cell(axis,it.index)))*d<0)
                face_velocities(FACE_INDEX<TV::m>(axis,mac_grid.First_Face_Index_In_Cell(axis,it.index)))=T(0);
            // if(face_velocities(FACE_INDEX<TV::m>(axis,mac_grid.Second_Face_Index_In_Cell(axis,it.index)))*d<0)
                face_velocities(FACE_INDEX<TV::m>(axis,mac_grid.Second_Face_Index_In_Cell(axis,it.index)))=T(0);}}
    // check whether divergence free
    Build_Velocity_Divergence();
    LOG::cout<<"Maximum velocity divergence after projection: "<<max_div<<std::endl;
}

//#####################################################################
// Function Update_Particle_Velocities
//#####################################################################
template<class TV> void MPMAC<TV>::
Update_Particle_Velocities()
{
#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        for(int axis=0;axis<TV::m;axis++){
            T V_PIC=(T)0;
            T V_FLIP=particles.V(p)(axis);
            for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+IN));it.Valid();it.Next()){
                TV_INT ind=influence_corner[axis](p)+it.index;
                FACE_INDEX<TV::m> face_ind(axis,ind);
                T w=weight[axis](p)(it.index);
                V_PIC+=face_velocities(face_ind)*w;
                V_FLIP+=(face_velocities(face_ind)-face_velocities_old(face_ind))*w;}
            particles.V(p)(axis)=((T)1-FLIP_alpha)*V_PIC+FLIP_alpha*V_FLIP;}}
}

//#####################################################################
// Function Particle_Based_Body_Collisions
//#####################################################################
template<class TV> void MPMAC<TV>::
Particle_Based_Body_Collisions()
{
#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        TV x=particles.X(p)+dt*particles.V(p); // candidate position
        for(int d=0;d<TV::m;d++){
            T left_wall=grid.domain.min_corner(d)+4.01*grid.dX.Min(),right_wall=grid.domain.max_corner(d)-4.01*grid.dX.Min();
            if(x(d)<=left_wall){
                TV vt(particles.V(p));vt(d)=(T)0;
                particles.V(p)=vt;}
            if(x(d)>=right_wall){
                TV vt(particles.V(p));vt(d)=(T)0;
                particles.V(p)=vt;}}}
}

//#####################################################################
// Function Update_Particle_Positions
//#####################################################################
template<class TV> void MPMAC<TV>::
Update_Particle_Positions()
{
#pragma omp parallel for
    for(int p=0;p<particles.number;p++){
        if(dt*particles.V(p).Max_Mag()>grid.dX.Min()){
            LOG::cout<<"breaking CFL"<<std::endl;
            PHYSBAM_FATAL_ERROR();}
        particles.X(p)+=dt*particles.V(p);}
}

//#####################################################################
// Function Get_Total_Momentum_On_Faces
//#####################################################################
template<class TV> TV MPMAC<TV>::
Get_Total_Momentum_On_Faces() const
{
    TV momen;
    for(FACE_ITERATOR<TV> iterator(mac_grid);iterator.Valid();iterator.Next()){
        FACE_INDEX<TV::m> face_index=iterator.Full_Index();
        int axis=iterator.Axis();
        momen(axis)+=face_masses(face_index)*face_velocities(face_index);}
    return momen;
}
    
//#####################################################################
// Function Get_Total_Momentum_On_Particles
//#####################################################################
template<class TV> TV MPMAC<TV>::
Get_Total_Momentum_On_Particles() const
{
    TV momen;
    for(int p=0;p<particles.number;p++)
        momen+=particles.mass(p)*particles.V(p);
    return momen;
}
    
//#####################################################################
template class MPMAC<VECTOR<float,2> >;
template class MPMAC<VECTOR<float,3> >;
template class MPMAC<VECTOR<double,2> >;
template class MPMAC<VECTOR<double,3> >;
}

