//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENSUBDIV_SURFACE
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/IDENTITY_ARRAY.h>
#include <Core/Matrices/BANDED_MATRIX.h>
#include <Core/Matrices/MATRIX.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <sstream>
#include <string>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int gauss_order> OPENSUBDIV_SURFACE<TV,gauss_order>::
OPENSUBDIV_SURFACE()
    :need_destroy_particles(true),particles(*new GEOMETRY_PARTICLES<TV>)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int gauss_order> OPENSUBDIV_SURFACE<TV,gauss_order>::
OPENSUBDIV_SURFACE(GEOMETRY_PARTICLES<TV>& particles)
    :need_destroy_particles(false),particles(particles)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int gauss_order> OPENSUBDIV_SURFACE<TV,gauss_order>::
~OPENSUBDIV_SURFACE()
{
    if(need_destroy_particles) delete &particles;
}
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,int gauss_order> OPENSUBDIV_SURFACE<TV,gauss_order>* OPENSUBDIV_SURFACE<TV,gauss_order>::
Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const
{
    OPENSUBDIV_SURFACE* surf=Create(new_particles);
    surf->need_destroy_particles=false;
    int offset=new_particles.Size();
    new_particles.Append(particles);
    if(particle_indices) particle_indices->Append_Elements(IDENTITY_ARRAY<>(particles.Size())+offset);
    surf->control_points=control_points+offset;
    surf->m=m;
    surf->mesh=mesh;
    surf->face_data=face_data;
    surf->thickness=thickness;
    return surf;
}
template<class TV,int gauss_order> void OPENSUBDIV_SURFACE<TV,gauss_order>::FACE_DATA::
Read(TYPED_ISTREAM& input)
{
    int size;
    Read_Binary(input,size);
    nodes.Clean_Memory();
    nodes.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,nodes.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,nodes.Get_Array_Pointer(),size);

    int size_confirm;
    Read_Binary(input,size_confirm);
    PHYSBAM_ASSERT(size==size_confirm);
    w.Clean_Memory();
    w.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,w.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,w.Get_Array_Pointer(),size);

    Read_Binary(input,size_confirm);
    PHYSBAM_ASSERT(size==size_confirm);
    A.Clean_Memory();
    A.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,A.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,A.Get_Array_Pointer(),size);
}
template<class TV,int gauss_order> void OPENSUBDIV_SURFACE<TV,gauss_order>::FACE_DATA::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,nodes,w,A);
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV,int gauss_order> void OPENSUBDIV_SURFACE<TV,gauss_order>::
Read(TYPED_ISTREAM& input)
{
    int num_verts;
    Read_Binary(input,m,num_verts,gauss_order);

    mesh.Clean_Memory();
    mesh.Resize(m);
    int size;
    Read_Binary(input,size);
    PHYSBAM_ASSERT(size==m);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,mesh.Get_Array_Pointer(),m);
    else Read_Binary_Array<float>(input.stream,mesh.Get_Array_Pointer(),m);

    for(int i=0;i<m;i++)
        for(int j=0;j<4;j++)
            PHYSBAM_ASSERT(mesh(i)(j)<num_verts);

    particles.Clean_Memory();
    Read_Binary(input,size);
    particles.Resize(size);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,particles.X.Get_Array_Pointer(),size);
    else Read_Binary_Array<float>(input.stream,particles.X.Get_Array_Pointer(),size);
    
    control_points.Clean_Memory();
    control_points.Resize(num_verts);
    Read_Binary(input,size);
    PHYSBAM_ASSERT(size==num_verts);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,control_points.Get_Array_Pointer(),num_verts);
    else Read_Binary_Array<float>(input.stream,control_points.Get_Array_Pointer(),num_verts);

    Read_Binary(input,size);
    PHYSBAM_ASSERT(size==m);
    face_data.Exact_Resize(m);
    if(input.type.use_doubles) Read_Binary_Array<double>(input.stream,face_data.Get_Array_Pointer(),m);
    else Read_Binary_Array<float>(input.stream,face_data.Get_Array_Pointer(),m);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV,int gauss_order> void OPENSUBDIV_SURFACE<TV,gauss_order>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,m,control_points.m,gauss_order);
    Write_Binary(output,mesh);
    Write_Binary(output,particles.X);
    Write_Binary(output,control_points);
    Write_Binary(output,face_data);
}
namespace{
double quadrature_weight[8][7]=
{
    {},
    {1},
    {0.5,0.5},
    {0.444444444444444444,0.277777777777777778,0.277777777777777778},
    {0.173927422568726929,0.173927422568726929,0.326072577431273071,0.326072577431273071},
    {0.284444444444444444,0.239314335249683234,0.239314335249683234,0.118463442528094544,0.118463442528094544},
    {0.0856622461895851732,0.0856622461895851732,0.233956967286345525,0.233956967286345525,0.180380786524069304,0.180380786524069304},
    {0.208979591836734694,0.0647424830844348507,0.0647424830844348507,0.190915025252559459,0.190915025252559459,0.139852695744638342,0.139852695744638342}
};
}
template<class TV,int gauss_order> void OPENSUBDIV_SURFACE<TV,gauss_order>::
Compute_G0()
{
    for(int face=0;face<m;face++){
        const ARRAY<int> nodes(control_points.Subset(face_data(face).nodes));
        const ARRAY<MATRIX<VECTOR<T,5>,gauss_order> >& A=face_data(face).A;
        
        for(int i=0;i<gauss_order;i++){
            for(int j=0;j<gauss_order;j++){
                VECTOR<TM,3>& G0_inv=face_data(face).G0_inv(i,j);
                VECTOR<T,3>& G0_det=face_data(face).G0_det(i,j);
                
                VECTOR<TV,5> a;
                for(int five=0;five<5;five++)
                    for(int al=0;al<nodes.m;al++)
                        a(five)+=A(al)(i,j)(five)*particles.X(nodes(al));
                
                TV a3=a(0).Cross(a(1));
                TV da31=a(0).Cross(a(3))+a(2).Cross(a(1));
                TV da32=a(0).Cross(a(4))+a(3).Cross(a(1));
                
                T m_sqr=a3.Magnitude_Squared();
                T dm_sqr1=2*da31.Dot(a3);
                T dm_sqr2=2*da32.Dot(a3);
                
                T J0=sqrt(m_sqr);
                T beta=J0/m_sqr;
                
                TV lambda_n=beta*a3;
                TV dlambda_n1=beta*(da31-(dm_sqr1/m_sqr)*a3);
                TV dlambda_n2=beta*(da31-(dm_sqr2/m_sqr)*a3);
                
                const VECTOR<T,3> simpson_weights=TV(1,4,1);
                for(int k=0;k<3;k++){
                    G0_inv(k)=TM(a(0)+(k-1)*thickness/2*dlambda_n1,a(1)+(k-1)*thickness/2*dlambda_n2,lambda_n);
                    G0_det(k)=abs(G0_inv(k).Determinant());
                    G0_inv(k).Invert();}}}}
}
//#####################################################################
// Function Set_Mass
//#####################################################################
template<class TV,int gauss_order> void OPENSUBDIV_SURFACE<TV,gauss_order>::
Set_Mass(T density,bool use_constant_mass) const
{
    ARRAY_VIEW<T>* mass_attr=particles.template Get_Array<T>(ATTRIBUTE_ID_MASS);
    PHYSBAM_ASSERT(mass_attr);
    if(use_constant_mass&&control_points.m)
        PHYSBAM_NOT_IMPLEMENTED(); // would need Total_Size() for this class.
    else{
        mass_attr->Subset(control_points).Fill((T)0);
    for(int face=0;face<m;face++){
        const ARRAY<int> nodes(control_points.Subset(face_data(face).nodes));
        for(int i=0;i<gauss_order;i++){
            for(int j=0;j<gauss_order;j++){
                const ARRAY<MATRIX<T,gauss_order> >& w=face_data(face).w;
                const VECTOR<T,3>& G0_det=face_data(face).G0_det(i,j);
                T mass=0;
                const VECTOR<T,3> simpson_weights=TV(1,4,1);
                for(int k=0;k<3;k++)
                    mass+=simpson_weights(k)*G0_det(k);
                
                mass*=density*thickness*quadrature_weight[gauss_order][i]*quadrature_weight[gauss_order][j]/6;
                PHYSBAM_ASSERT(mass>0); // will fail here if G0 hasn't been computed yet.
                
                for(int a=0;a<nodes.m;a++)
                    (*mass_attr)(nodes(a))+=mass*w(a)(i,j);}}}}
}
//#####################################################################
// Function Create_Triangulated_Surface
//#####################################################################
template<class TV,int gauss_order> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,2>::OBJECT* PhysBAM::
Create_Triangulated_Surface(const OPENSUBDIV_SURFACE<TV,gauss_order>& surf,bool same_particles)
{
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,2>::OBJECT T_TRIANGULATED_SURFACE;
    T_TRIANGULATED_SURFACE* tri_patch=0;
    if(same_particles) tri_patch=T_TRIANGULATED_SURFACE::Create(surf.particles);
    else tri_patch=T_TRIANGULATED_SURFACE::Create();
    for(int f=0;f<surf.m;f++){
        int a00=surf.control_points(surf.mesh(f)(0)),a01=surf.control_points(surf.mesh(f)(1)),
            a10=surf.control_points(surf.mesh(f)(3)),a11=surf.control_points(surf.mesh(f)(2));
        tri_patch->mesh.elements.Append(VECTOR<int,3>(a00,a01,a11));
        tri_patch->mesh.elements.Append(VECTOR<int,3>(a00,a11,a10));}

    if(!same_particles){
        ARRAY<int> map(surf.particles.X.m,true,-1);
        ARRAY_VIEW<int> av=tri_patch->mesh.elements.Flattened();
        int next=0;
        LOG::printf("av.m = %P\n",av.m);
        for(int i=0;i<av.m;i++)
            if(map(av(i))<0){
                map(av(i))=next++;
                tri_patch->particles.Add_Element();
                tri_patch->particles.X.Last()=surf.particles.X(av(i));}
        av=map.Subset(av);}
    tri_patch->Update_Number_Nodes();
    return tri_patch;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int gauss_order> OPENSUBDIV_SURFACE<TV,gauss_order>* OPENSUBDIV_SURFACE<TV,gauss_order>::
Create()
{
    return new OPENSUBDIV_SURFACE;
}
//#####################################################################
// Function Create
//#####################################################################
template<class TV,int gauss_order> OPENSUBDIV_SURFACE<TV,gauss_order>* OPENSUBDIV_SURFACE<TV,gauss_order>::
Create(GEOMETRY_PARTICLES<TV>& particles)
{
    return new OPENSUBDIV_SURFACE(particles);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV,int gauss_order> void OPENSUBDIV_SURFACE<TV,gauss_order>::
Initialize(const std::string& filename,T thickness_in)
{
    std::istream* input_raw=FILE_UTILITIES::Safe_Open_Input(filename);
    TYPED_ISTREAM input(*input_raw,STREAM_TYPE((RW())));
    Read(input);
    delete input_raw;
    thickness=thickness_in;
    Compute_G0();
}
//#####################################################################
// Function Evaluate
//#####################################################################
template<class TV,int gauss_order> TV OPENSUBDIV_SURFACE<TV,gauss_order>::
Evaluate(T s,T t,VECTOR<TV,2>* tangents) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV,int gauss_order> std::string OPENSUBDIV_SURFACE<TV,gauss_order>::
Name() const
{
    return Static_Name();
}
//#####################################################################
// Function Static_Name
//#####################################################################
template<class TV,int gauss_order> std::string OPENSUBDIV_SURFACE<TV,gauss_order>::
Static_Name()
{
    return LOG::sprintf("OPENSUBDIV_SURFACE<VECTOR<T,%d>>",TV::dimension);
}
//#####################################################################
// Function Extension
//#####################################################################
template<class TV,int gauss_order> std::string OPENSUBDIV_SURFACE<TV,gauss_order>::
Extension() const
{
    return Static_Extension();
}
//#####################################################################
// Function Static_Extension
//#####################################################################
template<class TV,int gauss_order> std::string OPENSUBDIV_SURFACE<TV,gauss_order>::
Static_Extension()
{
    return "";
}
namespace PhysBAM{
template class OPENSUBDIV_SURFACE<VECTOR<float,3>>;
template class OPENSUBDIV_SURFACE<VECTOR<double,3>>;

template TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<float,3>,2>::OBJECT* Create_Triangulated_Surface<VECTOR<float,3>,3>(OPENSUBDIV_SURFACE<VECTOR<float,3>,3> const&,bool);
template TOPOLOGY_BASED_SIMPLEX_POLICY<VECTOR<double,3>,2>::OBJECT* Create_Triangulated_Surface<VECTOR<double,3>,3>(OPENSUBDIV_SURFACE<VECTOR<double,3>,3> const&,bool);
}
