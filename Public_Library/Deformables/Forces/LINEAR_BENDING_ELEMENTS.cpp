//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Arrays/CONSTANT_ARRAY.h>
#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Data_Structures/SPARSE_UNION_FIND.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cyclic_shift.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Deformables/Forces/LINEAR_BENDING_ELEMENTS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LINEAR_BENDING_ELEMENTS<TV>::
LINEAR_BENDING_ELEMENTS(DEFORMABLE_PARTICLES<TV>& particles,T_MESH& mesh)
    :DEFORMABLES_FORCES<TV>(particles),mesh(mesh),stiffness(0),damping(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LINEAR_BENDING_ELEMENTS<TV>::
~LINEAR_BENDING_ELEMENTS()
{}
//#####################################################################
// Function Compute_Stiffness_Matrix
//#####################################################################
namespace{
template<class T,class TV> void Compute_Stiffness_Matrix_Helper(SEGMENT_MESH& mesh,ARRAY_VIEW<const TV> X,ARRAY<T>& stiffness_matrix_diagonal,SPARSE_MATRIX_FLAT_MXN<T>& stiffness_matrix_upper)
{
    if(!mesh.neighbor_nodes) mesh.Initialize_Neighbor_Nodes();

    // compute row lengths
    ARRAY<int> row_lengths(mesh.number_nodes);
    for(int s=0;s<mesh.elements.m;s++) row_lengths(mesh.elements(s).Min())++;
    for(int p=0;p<mesh.number_nodes;p++){const ARRAY<int>& neighbors=(*mesh.neighbor_nodes)(p);
        for(int i=0;i<neighbors.m;i++) for(int j=i+1;j<neighbors.m;j++) row_lengths(min(neighbors(i),neighbors(j)))++;}

    // ensure no empty rows
    ARRAY<int> empty_rows;for(int p=0;p<mesh.number_nodes;p++) if(!row_lengths(p)) empty_rows.Append(p);
    row_lengths.Subset(empty_rows)+=1;
    stiffness_matrix_upper.Set_Row_Lengths(row_lengths);
    for(int i=0;i<empty_rows.m;i++){int p=empty_rows(i);stiffness_matrix_upper(p,p)=0;}

    // compute stiffness matrix
    stiffness_matrix_diagonal=CONSTANT_ARRAY<T>(mesh.number_nodes,0);
    for(int p=0;p<mesh.number_nodes;p++){const ARRAY<int>& neighbors=(*mesh.neighbor_nodes)(p);
        for(int i=0;i<neighbors.m;i++) for(int j=i+1;j<neighbors.m;j++){
            const VECTOR<int,3> nodes(neighbors(i),p,neighbors(j));
            TV X0=X(nodes[0]),X1=X(nodes[1]),X2=X(nodes[2]);
            TV e01=X1-X0,e12=X2-X1;
            T length01=e01.Magnitude(),length12=e12.Magnitude();
            // If edge lengths remain constant, we have |sin psi/2| = |(X0-X1)/length01 + (X2-X1)/length12|.
            // We'll define energy as stiffness/(2*length_scale)*(sin psi/2)^2
            T scale=length01+length12;
            VECTOR<T,3> c(1/length01,-1/length01-1/length12,1/length12);
            for(int i=0;i<nodes.m;i++){
                stiffness_matrix_diagonal(nodes[i])+=scale*sqr(c[i]);
                for(int j=i+1;j<nodes.m;j++){int a=nodes[i],b=nodes[j];exchange_sort(a,b);
                    stiffness_matrix_upper(a,b)+=scale*c[i]*c[j];}}}}
}
template<class T,class TV> void Compute_Stiffness_Matrix_Helper(TRIANGLE_MESH& mesh,ARRAY_VIEW<const TV> X,ARRAY<T>& stiffness_matrix_diagonal,SPARSE_MATRIX_FLAT_MXN<T>& stiffness_matrix_upper)
{
    // compute bending quadruples
    if(!mesh.adjacent_elements) mesh.Initialize_Adjacent_Elements();
    ARRAY<VECTOR<int,4> > bending_quadruples;
    for(int t=0;t<mesh.elements.m;t++){
        int t1,t2,t3;mesh.elements(t).Get(t1,t2,t3);
        for(int a=0;a<(*mesh.adjacent_elements)(t).m;a++){
            int s=(*mesh.adjacent_elements)(t)(a);
            if(s>t){
                int s1,s2,s3;mesh.elements(s).Get(s1,s2,s3);
                if(t1==s1 || t1==s2 || t1==s3){cyclic_shift(t1,t2,t3);if(t1==s1 || t1==s2 || t1==s3) cyclic_shift(t1,t2,t3);}
                bending_quadruples.Append(VECTOR<int,4>(t1,t2,t3,mesh.Other_Node(t2,t3,s)));}}}

    // compute row lengths
    SEGMENT_MESH& segment_mesh=mesh.Get_Segment_Mesh();
    ARRAY<int> row_lengths(mesh.number_nodes);
    for(int s=0;s<segment_mesh.elements.m;s++) row_lengths(segment_mesh.elements(s).Min())++;
    for(int q=0;q<bending_quadruples.m;q++) row_lengths(min(bending_quadruples(q)[0],bending_quadruples(q)[3]))++;

    // ensure no empty rows
    ARRAY<int> empty_rows;for(int p=0;p<mesh.number_nodes;p++) if(!row_lengths(p)) empty_rows.Append(p);
    row_lengths.Subset(empty_rows)+=1;
    stiffness_matrix_upper.Set_Row_Lengths(row_lengths);
    for(int i=0;i<empty_rows.m;i++){int p=empty_rows(i);stiffness_matrix_upper(p,p)=0;}

    // compute stiffness matrix
    // for details see Wardetzky et al., "Discrete Quadratic Curvature Energies", Computer Aided Geometric Design, 2007
    stiffness_matrix_diagonal=CONSTANT_ARRAY<T>(mesh.number_nodes,0);
    for(int q=0;q<bending_quadruples.m;q++){const VECTOR<int,4>& nodes=bending_quadruples(q);
        TV X0=X(nodes[0]),X1=X(nodes[1]),X2=X(nodes[2]),X3=X(nodes[3]);
        TV e0=X2-X1,e1=X3-X1,e2=X0-X1,e3=X3-X2,e4=X0-X2; // edge numbering matches Wardetzky et al. p. 16
        T cross01=TV::Cross_Product(e0,e1).Magnitude(),cross02=TV::Cross_Product(e0,e2).Magnitude(),
            cross03=TV::Cross_Product(e0,e3).Magnitude(),cross04=TV::Cross_Product(e0,e4).Magnitude();
        T dot01=TV::Dot_Product(e0,e1),dot02=TV::Dot_Product(e0,e2),dot03=-TV::Dot_Product(e0,e3),dot04=-TV::Dot_Product(e0,e4);
        T cot01=dot01/cross01,cot02=dot02/cross02,cot03=dot03/cross03,cot04=dot04/cross04;
        T max_cot=max(cot01,cot02,cot03,cot04);
        T scale=3/(cross01+cross02);
        if(max_cot>5) scale*=25/sqr(max_cot);
        VECTOR<T,4> c;c[1]=cot03+cot04;c[2]=cot01+cot02;c[3]=-cot01-cot03;c[0]=-cot02-cot04; // node numbering matches bending quadruple
        for(int i=0;i<nodes.m;i++){
            stiffness_matrix_diagonal(nodes[i])+=scale*sqr(c[i]);
            for(int j=i+1;j<nodes.m;j++){int a=nodes[i],b=nodes[j];exchange_sort(a,b);
                stiffness_matrix_upper(a,b)+=scale*c[i]*c[j];}}}
}
}
template<class TV> void LINEAR_BENDING_ELEMENTS<TV>::
Compute_Stiffness_Matrix(ARRAY_VIEW<const TV> X)
{
    Compute_Stiffness_Matrix_Helper(mesh,X,stiffness_matrix_diagonal,stiffness_matrix_upper);

    LOG::cout<<"max diagonal element = "<<stiffness_matrix_diagonal.Max()<<std::endl;

    // verify that all edges occur in a triple/quadruple
    SEGMENT_MESH& segment_mesh=mesh.Get_Segment_Mesh();
    for(int s=0;s<segment_mesh.elements.m;s++){int i,j;segment_mesh.elements(s).Get(i,j);
        if(i<j && !stiffness_matrix_upper.Element_Present(i,j)) PHYSBAM_FATAL_ERROR("not all edges occur in bending quadruples");}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void LINEAR_BENDING_ELEMENTS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void LINEAR_BENDING_ELEMENTS<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> void LINEAR_BENDING_ELEMENTS<TV>::
Add_Force(const T scale,ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> F) const
{
    if(!stiffness_matrix_diagonal.Size()) PHYSBAM_FATAL_ERROR();

    T twice_scale=2*scale;
    // diagonal
    for(int p=0;p<X.Size();p++) F(p)-=twice_scale*stiffness_matrix_diagonal(p)*X(p);
    // offdiagonal
    const ARRAY<int>& offsets=stiffness_matrix_upper.offsets;
    const ARRAY<SPARSE_MATRIX_ENTRY<T> >& A=stiffness_matrix_upper.A;
    int index=offsets(0);
    for(int i=1;i<offsets.Size();i++){
        int end=offsets(i+1);
        for(;index<end;index++){
            int j=A(index).j;T twice_entry=twice_scale*A(index).a;
            F(i)-=twice_entry*X(j);F(j)-=twice_entry*X(i);}}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LINEAR_BENDING_ELEMENTS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    Add_Force(stiffness,particles.X,F);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void LINEAR_BENDING_ELEMENTS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    PHYSBAM_ASSERT(V.Size()==particles.Size());
    Add_Force(damping,V,F);
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_BENDING_ELEMENTS<TV>::
Compute_Energy() const
{
    if(!stiffness_matrix_diagonal.Size()) PHYSBAM_FATAL_ERROR();

    T energy=0;
    ARRAY_VIEW<const TV> X(particles.X);
    // diagonal
    for(int p=0;p<X.Size();p++) energy+=stiffness*stiffness_matrix_diagonal(p)*X(p).Magnitude_Squared();
    // offdiagonal
    const ARRAY<int>& offsets=stiffness_matrix_upper.offsets;
    const ARRAY<SPARSE_MATRIX_ENTRY<T> >& A=stiffness_matrix_upper.A;
    int index=offsets(0);
    for(int i=1;i<offsets.Size();i++){
        int end=offsets(i+1);
        for(;index<end;index++){
            int j=A(index).j;T entry=stiffness*A(index).a;
            energy+=2*entry*TV::Dot_Product(X(i),X(j));}}
    return energy;
}
//#####################################################################
namespace PhysBAM{
template class LINEAR_BENDING_ELEMENTS<VECTOR<float,2> >;
template class LINEAR_BENDING_ELEMENTS<VECTOR<float,3> >;
template class LINEAR_BENDING_ELEMENTS<VECTOR<double,2> >;
template class LINEAR_BENDING_ELEMENTS<VECTOR<double,3> >;
}
