//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <Incompressible/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> POISSON_COLLIDABLE_UNIFORM<TV>::
POISSON_COLLIDABLE_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input)
    :BASE(grid_input,u_input,initialize_grid,multiphase_input,enforce_compatibility_input),levelset_multiple(0),
    levelset_multiple_default(grid,phis_default,false),dt(0),dt_is_set(false)
{
    Initialize_Grid(grid_input);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> POISSON_COLLIDABLE_UNIFORM<TV>::
POISSON_COLLIDABLE_UNIFORM(const GRID<TV>& grid_input,ARRAY<T,TV_INT>& u_input,LEVELSET<TV>& cell_centered_levelset,const bool initialize_grid,const bool multiphase_input,
    const bool enforce_compatibility_input)
    :BASE(grid_input,u_input,initialize_grid,multiphase_input,enforce_compatibility_input),levelset_multiple(0),
    levelset_multiple_default(grid,phis_default,false),dt(0),dt_is_set(false)
{
    levelset=&cell_centered_levelset;
    Initialize_Grid(grid_input);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> POISSON_COLLIDABLE_UNIFORM<TV>::
~POISSON_COLLIDABLE_UNIFORM()
{
}
//#####################################################################
// Function Compute_beta_And_Add_Jumps_To_b
//#####################################################################
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Compute_beta_And_Add_Jumps_To_b(const T dt,const T time)
{
    int ghost_cells=3;
    if(!multiphase){
        ARRAY<T,TV_INT> phi_ghost;
        if((!use_variable_beta && !beta_given_on_faces) || u_jumps || beta_un_jumps){
            assert(levelset);phi_ghost.Resize(grid.Domain_Indices(ghost_cells),false);levelset->boundary->Fill_Ghost_Cells(grid,levelset->phi,phi_ghost,dt,time,ghost_cells);}
        if(!beta_given_on_faces){if(use_variable_beta) Find_Variable_beta();else Find_Constant_beta(phi_ghost);}
        if(u_jumps) Add_Jump_To_b(phi_ghost);
        if(beta_un_jumps) Add_Derivative_Jump_To_b(phi_ghost);}
    else{
        ARRAY<ARRAY<T,TV_INT>> phis_ghost;
        if((!use_variable_beta && !beta_given_on_faces) || u_jumps || beta_un_jumps){assert(levelset_multiple);
            phis_ghost.Resize(levelset_multiple->phis.m);
            for(int i=0;i<phis_ghost.m;i++){phis_ghost(i).Resize(grid.Domain_Indices(ghost_cells),false);
                levelset_multiple->levelsets(i)->boundary->Fill_Ghost_Cells(grid,levelset_multiple->phis(i),phis_ghost(i),dt,time,ghost_cells);}}
        if(!beta_given_on_faces){if(use_variable_beta) Find_Variable_beta();else Find_Constant_beta_Multiphase(phis_ghost);}
        if(u_jumps) Add_Jump_To_b_Multiphase(phis_ghost);
        if(beta_un_jumps) PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Find_Constant_beta
//#####################################################################
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Find_Constant_beta(const ARRAY<T,TV_INT>& phi_ghost)
{
    Find_Constant_beta(beta_face,phi_ghost);
}
//#####################################################################
// Function Find_Constant_beta
//#####################################################################
// only set up for jump conditons - doesn't work for Dirichlet boundary conditions yet 
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Find_Constant_beta(ARRAY<T,FACE_INDEX<TV::m> >& beta_face,const ARRAY<T,TV_INT>& phi_ghost)
{
    assert(!multiphase);
    
    // interior
    T half_width=(T).5*number_of_interface_cells*grid.dX.Max();
    if(GFM || smear_beta){
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
            beta_face.Component(iterator.Axis())(iterator.Face_Index())=
                LEVELSET_UTILITIES<T>::Heaviside((T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index())),beta_minus,beta_plus,half_width);}
    else{ // smear 1/beta for the delta function method
        T rho_minus=1/beta_minus,rho_plus=1/beta_plus;
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
            beta_face.Component(iterator.Axis())(iterator.Face_Index())=
                1/LEVELSET_UTILITIES<T>::Heaviside((T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index())),rho_minus,rho_plus,half_width);}

    if(GFM){ // adjust beta near interface
        T beta_minus_times_beta_plus=beta_minus*beta_plus;
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(iterator.First_Cell_Index()),phi_ghost(iterator.Second_Cell_Index())))
            beta_face.Component(iterator.Axis())(iterator.Face_Index())=
                beta_minus_times_beta_plus/LEVELSET_UTILITIES<T>::Convex_Average(phi_ghost(iterator.First_Cell_Index()),phi_ghost(iterator.Second_Cell_Index()),beta_minus,beta_plus);}
}
//#####################################################################
// Function Find_Constant_beta_Multiphase
//#####################################################################
// only set up for jump conditons - doesn't work for Dirichlet boundary conditions yet 
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Find_Constant_beta_Multiphase(ARRAY<ARRAY<T,TV_INT>>& phis_ghost)
{
    LEVELSET_MULTIPLE<TV> levelset(grid,phis_ghost);levelset.Recreate_Levelsets();

    T half_width=(T).5*number_of_interface_cells*grid.dX.Max();
    if(GFM || smear_beta){
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            // use the beta from the non dirichlet region for dirichlet boundary conditions
            if(psi_D(iterator.First_Cell_Index())) beta_face.Component(iterator.Axis())(iterator.Face_Index())=levelset.Region_Value(iterator.Second_Cell_Index(),beta_multiphase); 
            else if(psi_D(iterator.Second_Cell_Index())) beta_face.Component(iterator.Axis())(iterator.Face_Index())=levelset.Region_Value(iterator.First_Cell_Index(),beta_multiphase);
            else beta_face.Component(iterator.Axis())(iterator.Face_Index())=levelset.Heaviside(iterator.First_Cell_Index(),iterator.Second_Cell_Index(),beta_multiphase,half_width);}}
    else{ // smear 1/beta for the delta function method
        ARRAY<T> rho_multiphase(beta_multiphase.m);for(int i=0;i<rho_multiphase.m;i++)rho_multiphase(i)=1/beta_multiphase(i);
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            // use the beta from the non dirichlet region for dirichlet boundary conditions
            if(psi_D(iterator.First_Cell_Index())) beta_face.Component(iterator.Axis())(iterator.Face_Index())=levelset.Region_Value(iterator.Second_Cell_Index(),beta_multiphase); 
            else if(psi_D(iterator.Second_Cell_Index())) beta_face.Component(iterator.Axis())(iterator.Face_Index())=levelset.Region_Value(iterator.First_Cell_Index(),beta_multiphase);
            else beta_face.Component(iterator.Axis())(iterator.Face_Index())=1/levelset.Heaviside(iterator.First_Cell_Index(),iterator.Second_Cell_Index(),rho_multiphase,half_width);}}
    if(GFM){ // adjust beta near interface
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) 
            if(!(psi_D(iterator.First_Cell_Index())||psi_D(iterator.Second_Cell_Index()))&&levelset.Interface(iterator.First_Cell_Index(),iterator.Second_Cell_Index())){
                T denominator=levelset.Convex_Average(iterator.First_Cell_Index(),iterator.Second_Cell_Index(),beta_multiphase);
                if(denominator == 0) beta_face.Component(iterator.Axis())(iterator.Face_Index())=0; // fix for overloading this class for implicit viscosity (with viscosity=0)
                else beta_face.Component(iterator.Axis())(iterator.Face_Index())=levelset.Region_Value(iterator.First_Cell_Index(),beta_multiphase)*
                    levelset.Region_Value(iterator.Second_Cell_Index(),beta_multiphase)/denominator;}}
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    BASE::Find_A_Part_Two(domain,A_array,b_array,cell_index_to_matrix_index);
    if(second_order_cut_cell_method) Apply_Second_Order_Cut_Cell_Method(domain,A_array,b_array,cell_index_to_matrix_index);
}
//#####################################################################
// Function Add_Jump_To_b
//#####################################################################
// b is negative
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Add_Jump_To_b(const ARRAY<T,TV_INT>& phi_ghost)
{
    assert(!multiphase);
    TV one_over_dx2=Inverse(grid.dX*grid.dX);

    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index(),face_index=iterator.Face_Index();int axis=iterator.Axis();
        if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(first_cell_index),phi_ghost(second_cell_index)) && !psi_N.Component(axis)(face_index) && !(psi_D(first_cell_index) && psi_D(second_cell_index))){
            T jump=beta_face.Component(axis)(face_index)*one_over_dx2[axis]*
                LEVELSET_UTILITIES<T>::Average(phi_ghost(first_cell_index),u_jump(first_cell_index),phi_ghost(second_cell_index),u_jump(second_cell_index));
            f(first_cell_index)-=LEVELSET_UTILITIES<T>::Sign(phi_ghost(first_cell_index))*jump;f(second_cell_index)-=LEVELSET_UTILITIES<T>::Sign(phi_ghost(second_cell_index))*jump;}}
}
//#####################################################################
// Function Add_Jump_To_b_Multiphase
//#####################################################################
// b is negative
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Add_Jump_To_b_Multiphase(ARRAY<ARRAY<T,TV_INT>>& phis_ghost)
{
    TV one_over_dx2=Inverse(grid.dX*grid.dX);
    LEVELSET_MULTIPLE<TV> levelset(grid,phis_ghost);levelset.Recreate_Levelsets();
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index(),face_index=iterator.Face_Index();int axis=iterator.Axis();
        if(levelset.Interface(first_cell_index,second_cell_index) && !psi_N.Component(axis)(face_index) && !(psi_D(first_cell_index) && psi_D(second_cell_index))){
            int region_1,region_2;T phi_1,phi_2;levelset.Minimum_Regions(first_cell_index,second_cell_index,region_1,region_2,phi_1,phi_2);
            T jump=beta_face.Component(axis)(face_index)*one_over_dx2[axis]*u_jump_face.Component(axis)(face_index);
            f(first_cell_index)-=LEVELSET_MULTIPLE<TV>::Sign(region_1,region_2)*jump;f(second_cell_index)-=LEVELSET_MULTIPLE<TV>::Sign(region_2,region_1)*jump;}}
}
//#####################################################################
// Function Add_Derivative_Jump_To_b
//#####################################################################
// b is negative
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Add_Derivative_Jump_To_b(const ARRAY<T,TV_INT>& phi_ghost)
{
    assert(!multiphase);
    TV one_over_dx=Inverse(grid.dX);
    bool normals_defined=(levelset->normals!=0);levelset->Compute_Normals();

    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index(),face_index=iterator.Face_Index();int axis=iterator.Axis();
        if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(first_cell_index),phi_ghost(second_cell_index)) && !psi_N.Component(axis)(face_index) && !(psi_D(first_cell_index) && psi_D(second_cell_index))){
            T jump=LEVELSET_UTILITIES<T>::Sign(phi_ghost(first_cell_index))*beta_face.Component(axis)(face_index)*one_over_dx[axis]*
                  LEVELSET_UTILITIES<T>::Average(phi_ghost(first_cell_index),beta_un_jump(first_cell_index)*
                  (*levelset->normals)(first_cell_index)[axis],phi_ghost(second_cell_index),beta_un_jump(second_cell_index)*(*levelset->normals)(second_cell_index)[axis]);
            T theta=LEVELSET_UTILITIES<T>::Theta(phi_ghost(first_cell_index),phi_ghost(second_cell_index));
            f(first_cell_index)-=(1-theta)*jump/LEVELSET_UTILITIES<T>::Heaviside(phi_ghost(second_cell_index),beta_minus,beta_plus);
            f(second_cell_index)-=theta*jump/LEVELSET_UTILITIES<T>::Heaviside(phi_ghost(first_cell_index),beta_minus,beta_plus);}}

    if(!normals_defined){delete levelset->normals;levelset->normals=0;}
}
//#####################################################################
// Function Apply_Second_Order_Cut_Cell_Method
//#####################################################################
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Apply_Second_Order_Cut_Cell_Method(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<ARRAY<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index)
{
    assert(levelset);
    TV minus_one_over_dx_squared=(T)-1*Inverse(grid.dX*grid.dX);

    if(use_variable_beta && !beta_given_on_faces) for(int i=0;i<TV::dimension;i++){
        RANGE<TV_INT> face_domain(domain);
        for(int axis=0;axis<TV::dimension;axis++){
            if(face_domain.min_corner(axis)==grid.Domain_Indices(1).min_corner(axis)) face_domain.min_corner(axis)+=1;
            if(face_domain.max_corner(axis)==grid.Domain_Indices(1).max_corner(axis)) face_domain.max_corner(axis)-=1;}
        if(face_domain.min_corner(i)==grid.Domain_Indices().min_corner(i)) face_domain.min_corner(i)+=1;
        if(face_domain.max_corner(i)!=grid.Domain_Indices().max_corner(i)) face_domain.max_corner(i)+=1;
        for(FACE_ITERATOR<TV> iterator(grid,face_domain,i);iterator.Valid();iterator.Next()){
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index(),face_index=iterator.Face_Index();int axis=iterator.Axis();
            if(!psi_N.Component(axis)(face_index) && LEVELSET_UTILITIES<T>::Interface(levelset->phi(first_cell_index),levelset->phi(second_cell_index))){
                T theta=LEVELSET_UTILITIES<T>::Theta(levelset->phi(first_cell_index),levelset->phi(second_cell_index));
                if(psi_D(first_cell_index) && !psi_D(second_cell_index) && domain.Lazy_Inside_Half_Open(second_cell_index)){ // interface is to the negative side of second cell
                    int color=filled_region_colors(second_cell_index);int matrix_index=cell_index_to_matrix_index(second_cell_index);
                    T A_right_i=beta_face.Component(axis)(face_index)*minus_one_over_dx_squared[axis];A_array(color).Add_Element(matrix_index,matrix_index,-A_right_i);
                    b_array(color)(matrix_index)-=A_right_i*u(first_cell_index); 
                    T beta_ghost=(beta_interface_face.Component(axis)(face_index)-theta*variable_beta(second_cell_index))/max(1-theta,second_order_cut_cell_threshold);
                    T beta_new=(T).5*(beta_ghost+variable_beta(second_cell_index));
                    A_right_i=beta_new/max((1-theta),second_order_cut_cell_threshold)*minus_one_over_dx_squared[axis];
                    A_array(color).Add_Element(matrix_index,matrix_index,A_right_i);b_array(color)(matrix_index)+=A_right_i*u_interface.Component(axis)(face_index);}
                else if(psi_D(second_cell_index) && !psi_D(first_cell_index) && domain.Lazy_Inside_Half_Open(first_cell_index)){ // interface is to the positive side of first cell
                    int color=filled_region_colors(first_cell_index);int matrix_index=cell_index_to_matrix_index(first_cell_index);
                    T A_right_i=beta_face.Component(axis)(face_index)*minus_one_over_dx_squared[axis];A_array(color).Add_Element(matrix_index,matrix_index,-A_right_i);
                    b_array(color)(matrix_index)-=A_right_i*u(second_cell_index);
                    T beta_ghost=(beta_interface_face.Component(axis)(face_index)+(theta-1)*variable_beta(first_cell_index))/max(theta,second_order_cut_cell_threshold);
                    T beta_new=(T).5*(beta_ghost+variable_beta(first_cell_index));
                    A_right_i=beta_new/max(theta,second_order_cut_cell_threshold)*minus_one_over_dx_squared[axis];
                    A_array(color).Add_Element(matrix_index,matrix_index,A_right_i);b_array(color)(matrix_index)+=A_right_i*u_interface.Component(axis)(face_index);}}}}
    else for(int i=0;i<TV::dimension;i++){
        RANGE<TV_INT> face_domain(domain);
        for(int axis=0;axis<TV::dimension;axis++){
            if(face_domain.min_corner(axis)==grid.Domain_Indices(1).min_corner(axis)) face_domain.min_corner(axis)+=1;
            if(face_domain.max_corner(axis)==grid.Domain_Indices(1).max_corner(axis)) face_domain.max_corner(axis)-=1;}
        if(face_domain.min_corner(i)==grid.Domain_Indices().min_corner(i)) face_domain.min_corner(i)+=1;
        if(face_domain.max_corner(i)!=grid.Domain_Indices().max_corner(i)) face_domain.max_corner(i)+=1;
        for(FACE_ITERATOR<TV> iterator(grid,face_domain,i);iterator.Valid();iterator.Next()){
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index(),face_index=iterator.Face_Index();int axis=iterator.Axis();
            if(!psi_N.Component(axis)(face_index) && LEVELSET_UTILITIES<T>::Interface(levelset->phi(first_cell_index),levelset->phi(second_cell_index))){
                T theta=LEVELSET_UTILITIES<T>::Theta(levelset->phi(first_cell_index),levelset->phi(second_cell_index));
                if(psi_D(first_cell_index) && !psi_D(second_cell_index) && domain.Lazy_Inside_Half_Open(second_cell_index)){ // interface is to the negative side of second cell
                    int color=filled_region_colors(second_cell_index);int matrix_index=cell_index_to_matrix_index(second_cell_index);
                    T A_right_i=beta_face.Component(axis)(face_index)*minus_one_over_dx_squared[axis];A_array(color).Add_Element(matrix_index,matrix_index,-A_right_i);
                    b_array(color)(matrix_index)-=A_right_i*u(first_cell_index); 
                    A_right_i/=max((1-theta),second_order_cut_cell_threshold);
                    A_array(color).Add_Element(matrix_index,matrix_index,A_right_i);b_array(color)(matrix_index)+=A_right_i*u_interface.Component(axis)(face_index);}
                else if(psi_D(second_cell_index) && !psi_D(first_cell_index) && domain.Lazy_Inside_Half_Open(first_cell_index)){ // interface is to the positive side of first cell
                    int color=filled_region_colors(first_cell_index);int matrix_index=cell_index_to_matrix_index(first_cell_index);
                    T A_right_i=beta_face.Component(axis)(face_index)*minus_one_over_dx_squared[axis];A_array(color).Add_Element(matrix_index,matrix_index,-A_right_i);
                    b_array(color)(matrix_index)-=A_right_i*u(second_cell_index);
                    A_right_i/=max(theta,second_order_cut_cell_threshold);
                    A_array(color).Add_Element(matrix_index,matrix_index,A_right_i);b_array(color)(matrix_index)+=A_right_i*u_interface.Component(axis)(face_index);}}}}
}
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Use_Internal_Level_Set(const int number_of_regions)
{
    assert(multiphase);
    assert(number_of_regions>=2);
    phis_default.Resize(number_of_regions);
    for(int i=0;i<phis_default.m;i++)phis_default(i).Resize(grid.Domain_Indices(3));
    levelset_multiple=&levelset_multiple_default;
    levelset_multiple->Recreate_Levelsets();
}
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Update_Internal_Level_Set(LEVELSET_MULTIPLE<TV>& levelset_multiple_input)
{
    for(int k=0;k<levelset_multiple_input.phis.m;k++) levelset_multiple->phis(k)=levelset_multiple_input.phis(k);
}
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input)
{
    second_order_cut_cell_method=use_second_order_cut_cell_method_input;
    if(second_order_cut_cell_method){u_interface.Resize(grid);beta_interface_face.Resize(grid);}
    else{u_interface.Clean_Memory();beta_interface_face.Clean_Memory();}
}
template<class TV> void POISSON_COLLIDABLE_UNIFORM<TV>::
Initialize_Grid(const GRID<TV>& grid_input)
{
    BASE::Initialize_Grid(grid_input);
    if(u_jumps)u_jump.Resize(grid.Domain_Indices(1));
    if(beta_un_jumps)beta_un_jump.Resize(grid.Domain_Indices(1));
}
//#####################################################################
template class POISSON_COLLIDABLE_UNIFORM<VECTOR<float,1> >;
template class POISSON_COLLIDABLE_UNIFORM<VECTOR<float,2> >;
template class POISSON_COLLIDABLE_UNIFORM<VECTOR<float,3> >;
template class POISSON_COLLIDABLE_UNIFORM<VECTOR<double,1> >;
template class POISSON_COLLIDABLE_UNIFORM<VECTOR<double,2> >;
template class POISSON_COLLIDABLE_UNIFORM<VECTOR<double,3> >;
}
