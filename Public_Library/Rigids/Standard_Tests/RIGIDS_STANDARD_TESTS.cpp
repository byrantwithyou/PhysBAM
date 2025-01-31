//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_STANDARD_TESTS
//#####################################################################
#include <Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <Geometry/Basic_Geometry/BOWL.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/LINE_2D.h>
#include <Geometry/Basic_Geometry/RING.h>
#include <Geometry/Basic_Geometry/SMOOTH_GEAR.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Tessellation/BOWL_TESSELLATION.h>
#include <Geometry/Tessellation/CYLINDER_TESSELLATION.h>
#include <Geometry/Tessellation/GEAR_TESSELLATION.h>
#include <Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <Geometry/Tessellation/PLANE_TESSELLATION.h>
#include <Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <Geometry/Tessellation/RING_TESSELLATION.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Tessellation/TORUS_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Joints/JOINT_MESH.h>
#include <Rigids/Joints/POINT_JOINT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_STANDARD_TESTS<TV>::
RIGIDS_STANDARD_TESTS(STREAM_TYPE stream_type,const std::string& data_directory,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :stream_type(stream_type),data_directory(data_directory),rigid_body_collection(rigid_body_collection_input)
{
}
//#####################################################################
// Function Add_Analytic_Sphere
//#####################################################################
namespace{
template<class T,class TV>
static IMPLICIT_OBJECT<TV>* Add_Analytic_Sphere(const TV&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(SPHERE<TV>(TV(),scaling_factor));
}
//#####################################################################
// Function Add_Analytic_Cylinder
//#####################################################################
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,2> >* Add_Analytic_Cylinder(const VECTOR<T,2>&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<RANGE<VECTOR<T,2> > >(RANGE<VECTOR<T,2> >(VECTOR<T,2>(-scaling_factor/2,-scaling_factor/2),VECTOR<T,2>(scaling_factor/2,scaling_factor/2)));
}
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,3> >* Add_Analytic_Cylinder(const VECTOR<T,3>&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(CYLINDER<T>(VECTOR<T,3>(0,-scaling_factor,0),VECTOR<T,3>(0,scaling_factor,0),(T).1*scaling_factor));
}
//#####################################################################
// Function Add_Analytic_Ring
//#####################################################################
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,2> >* Add_Analytic_Ring(const VECTOR<T,2>&,const T scaling_factor)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,3> >* Add_Analytic_Ring(const VECTOR<T,3>&,const T scaling_factor)
{
    T half_height=(T).5*scaling_factor;
    return new ANALYTIC_IMPLICIT_OBJECT<RING<T> >(RING<T>(VECTOR<T,3>(0,-half_height,0),VECTOR<T,3>(0,half_height,0),(T)3*scaling_factor,(T)2*scaling_factor));
}
//#####################################################################
// Function Add_Analytic_Smooth_Gear
//#####################################################################
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,2> >* Add_Analytic_Smooth_Gear(T r,T s,int n)
{
    return new ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<VECTOR<T,2> > >(SMOOTH_GEAR<VECTOR<T,2> >(r,s,n));
}
template<class T>
static IMPLICIT_OBJECT<VECTOR<T,3> >* Add_Analytic_Smooth_Gear(T r,T s,T w,int n)
{
    return new ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<VECTOR<T,3> > >(SMOOTH_GEAR<VECTOR<T,3> >(r,s,n,w));
}
//#####################################################################
// Function Add_Analytic_Ground
//#####################################################################
template<class TV,class T>
static IMPLICIT_OBJECT<TV>* Add_Analytic_Ground(const TV&,const T scaling_factor)
{
    return new ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >(BOUNDED_HORIZONTAL_PLANE<TV>(100*scaling_factor));
}
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Rigid_Body(const std::string& name,const T scaling_factor,const T friction,const bool read_implicit,const bool always_read_object)
{
    std::string rigid_directory=data_directory+"/Rigid_Bodies"+(TV::m==3?"":"_2D"),basename=rigid_directory+"/"+name;
    bool read_simplicial=true;
    if(!always_read_object){
        IMPLICIT_OBJECT<TV>* implicit=0;
        if(name=="sphere") implicit=::Add_Analytic_Sphere(TV(),scaling_factor);
        else if(name=="skinnycyllink") implicit=::Add_Analytic_Cylinder(TV(),scaling_factor);
        else if(name=="Rings_Test/ring_revolve") implicit=Add_Analytic_Ring(TV(),scaling_factor);
        else if(name=="ground") implicit=Add_Analytic_Ground(TV(),scaling_factor);
        if(implicit){
            if(TV::m==3){
                STRUCTURE<TV>* structure=choice<TV::m-2>(implicit,TESSELLATION::Generate_Triangles(*implicit));
                if(!rigid_body_collection.Register_Analytic_Replacement_Structure(basename+".tri",scaling_factor,structure))
                    delete structure;}
            else
                read_simplicial=false;
            if(!rigid_body_collection.Register_Analytic_Replacement_Structure(basename+(TV::m==3?".phi":".phi2d"),scaling_factor,implicit))
                delete implicit;}}
    int id=rigid_body_collection.Add_Rigid_Body(basename,scaling_factor,read_simplicial,read_implicit);
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);
    rigid_body.coefficient_of_friction=friction;
    rigid_body.name=name;
    return rigid_body;
}
//#####################################################################
// Function Add_Ground
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Ground(const T friction,const T height,const T coefficient_of_restitution,const T scale)
{
    PHYSBAM_ASSERT(scale>0);
    RIGID_BODY<TV>& ground=Add_Rigid_Body("ground",scale,friction);
    ground.Frame().t.y=height;
    ground.is_static=true;
    ground.rigid_body_collection.collision_body_list->Get_Collision_Geometry(ground.particle_index)->add_to_spatial_partition=false;
    ground.Set_Coefficient_Of_Restitution(coefficient_of_restitution);
    ground.name="ground";
    return ground;
}
//#####################################################################
// Function Add_Analytic_Box
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Box(const VECTOR<T,1>& scaling_factor)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    RANGE<TV> box((T)-.5*scaling_factor,(T).5*scaling_factor);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(box));
    GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
    POINT_SIMPLICES_1D<T>& simplicial_object=*POINT_SIMPLICES_1D<T>::Create(particles);
    particles.Add_Elements(2);
    POINT_SIMPLEX_MESH& segment_mesh=simplicial_object.mesh;
    segment_mesh.number_nodes=2;
    segment_mesh.elements.Preallocate(1);
    particles.X(0)=VECTOR<T,1>(box.min_corner.x);
    particles.X(1)=VECTOR<T,1>(box.max_corner.x);
    simplicial_object.mesh.elements.Append(VECTOR<int,1>(0));
    simplicial_object.mesh.directions.Append(false);
    simplicial_object.mesh.elements.Append(VECTOR<int,1>(1));
    simplicial_object.mesh.directions.Append(true);
    rigid_body.Add_Structure(simplicial_object);
    simplicial_object.Update_Point_Simplex_List();
    assert(simplicial_object.point_simplex_list);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;   
}
//#####################################################################
// Function Add_Analytic_Box
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Box(const VECTOR<T,2>& scaling_factor,int segments_per_side,T density)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    RANGE<TV> box((T)-.5*scaling_factor,(T).5*scaling_factor);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(box));
    GEOMETRY_PARTICLES<TV>& particles=*new GEOMETRY_PARTICLES<TV>;
    SEGMENTED_CURVE_2D<T>& simplicial_object=*SEGMENTED_CURVE_2D<T>::Create(particles);
    particles.Add_Elements(4*segments_per_side);
    SEGMENT_MESH& segment_mesh=simplicial_object.mesh;
    segment_mesh.number_nodes=4*segments_per_side;
    segment_mesh.elements.Preallocate(4*segments_per_side);
    int last_node=4*segments_per_side-1;
    VECTOR<T,2> position=box.min_corner,e=box.Edge_Lengths(),dir[4]={{e.x,0},{0,e.y},{-e.x,0},{0,-e.y}};
    for(int side=0;side<4;side++)
        for(int vertex=0;vertex<segments_per_side;vertex++){
            int current_node=side*segments_per_side+vertex;
            particles.X(current_node)=position;
            position+=dir[side];
            simplicial_object.mesh.elements.Append(VECTOR<int,2>(last_node,current_node));
            last_node=current_node;}
    rigid_body.Add_Structure(simplicial_object);
    simplicial_object.Update_Segment_List();
    assert(simplicial_object.segment_list);
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    if(density){
        rigid_body.Mass()=density*box.Size();
        rigid_body.Inertia_Tensor()(0)=rigid_body.Mass()/12*box.Edge_Lengths().Magnitude_Squared();}
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Box
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Box(const VECTOR<T,3>& scaling_factor,const VECTOR<int,3>& resolution,typename TV::SCALAR density)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    RANGE<TV> box((T)-.5*scaling_factor,(T).5*scaling_factor);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(box));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(box,resolution));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    if(density){
        rigid_body.Mass()=density*box.Size();
        TV r2=box.Edge_Lengths()*box.Edge_Lengths();
        rigid_body.Inertia_Tensor()=rigid_body.Mass()/12*(r2.Sum()-DIAGONAL_MATRIX<T,3>(r2));}
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Torus
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Torus(const T inner_radius,const T outer_radius,int inner_resolution,int outer_resolution,T density)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    TORUS<T> torus(TV(),TV(0,0,1),inner_radius,outer_radius);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >(torus));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(torus,inner_resolution,outer_resolution));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    if(density){
        if(TV::m==2){
            rigid_body.Mass()=density*pi*(sqr(outer_radius)-sqr(inner_radius));
            rigid_body.Inertia_Tensor()(0)=rigid_body.Mass()*.5*(sqr(outer_radius)+sqr(inner_radius));}
        else if(TV::m==3){
            T a=0.5*(outer_radius-inner_radius),c=0.5*(outer_radius+inner_radius);
            rigid_body.Mass()=density*2*sqr(pi*a)*c;
            rigid_body.Inertia_Tensor()(0)=rigid_body.Inertia_Tensor()(1)=rigid_body.Mass()*(5*sqr(a)/8+.5*sqr(c));
            rigid_body.Inertia_Tensor()(2)=rigid_body.Mass()*(3*sqr(a)/4+sqr(c));}}
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Cylinder
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Cylinder(const T height,const T radius,int resolution_radius,int resolution_height,T density)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    CYLINDER<T> cylinder(TV(0,0,-height/2),TV(0,0,height/2),radius);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >(cylinder));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(cylinder,resolution_height,resolution_radius));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    if(density){
        if(TV::m==2){
            rigid_body.Mass()=density*2*radius*height;
            rigid_body.Inertia_Tensor()(0)=rigid_body.Mass()/12*(sqr(2*radius)+sqr(height));}
        else if(TV::m==3){
            rigid_body.Mass()=density*pi*sqr(radius)*height;
            rigid_body.Inertia_Tensor()(0)=rigid_body.Inertia_Tensor()(1)=rigid_body.Mass()/12*(3*sqr(radius)+sqr(height));
            rigid_body.Inertia_Tensor()(2)=rigid_body.Mass()/2*sqr(radius);}}
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Shell
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Shell(const T height,const T outer_radius,const T inner_radius,int resolution,T density)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    RING<T> ring(TV(0,0,-height/2),TV(0,0,height/2),outer_radius,inner_radius);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<RING<T> >(ring));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(ring,resolution));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    if(TV::m==2){
        rigid_body.Mass()=density*pi*(sqr(outer_radius)-sqr(inner_radius));
        rigid_body.Inertia_Tensor()(0)=rigid_body.Mass()*.5*(sqr(outer_radius)+sqr(inner_radius));}
    else if(TV::m==3){
        rigid_body.Mass()=density*pi*(sqr(outer_radius)-sqr(inner_radius))*height;
        rigid_body.Inertia_Tensor()(0)=rigid_body.Inertia_Tensor()(1)=
            rigid_body.Mass()/12*(3*(sqr(inner_radius)+sqr(outer_radius))+sqr(height));
        rigid_body.Inertia_Tensor()(2)=rigid_body.Mass()*0.5*(sqr(inner_radius)+sqr(outer_radius));}
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Bowl
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Bowl(const TV& location,const TV& axis,const T hole_radius,const T depth,const T thickness,int res_radial, int res_vertical)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);

    BOWL<T> bowl(location,axis,hole_radius,depth,thickness);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<BOWL<T> >(bowl));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(bowl,res_radial,res_vertical));

    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Smooth_Gear
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Smooth_Gear(const TV& dimensions,int cogs,int levels)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    SMOOTH_GEAR<TV> gear(dimensions,cogs);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SMOOTH_GEAR<TV> >(gear));
    rigid_body.Add_Structure(*TESSELLATION::Generate_Triangles(gear,levels));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Smooth_Gear
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Plane(const TV& normal,const TV& X,T surface_size,int elements_per_side)
{
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE T_OBJECT;
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    T_OBJECT plane(normal,X);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<T_OBJECT>(plane));
    rigid_body.Add_Structure(*TESSELLATION::Tessellate_Boundary(plane,surface_size,elements_per_side));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    return rigid_body;
}
//#####################################################################
// Function Add_Analytic_Sphere
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGIDS_STANDARD_TESTS<TV>::
Add_Analytic_Sphere(const T radius,const T density,int levels)
{
    RIGID_BODY<TV>& rigid_body=*new RIGID_BODY<TV>(rigid_body_collection,true);
    SPHERE<TV> sphere(TV(),radius);
    rigid_body.Add_Structure(*new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(sphere));
    rigid_body.Add_Structure(*TESSELLATION::Tessellate_Boundary(sphere,levels));
    rigid_body_collection.Add_Rigid_Body_And_Geometry(&rigid_body);
    T r2=radius*radius,m2d=r2*(T)pi*density;
    if(TV::m==2){
        rigid_body.Mass()=m2d;
        rigid_body.Inertia_Tensor()(0)=rigid_body.Mass()*r2*(T).5;}
    else if(TV::m==3){
        rigid_body.Mass()=m2d*radius*((T)4/3);
        rigid_body.Inertia_Tensor()*=0;
        rigid_body.Inertia_Tensor()+=rigid_body.Mass()*r2*(T).4;}
    return rigid_body;
}
//#####################################################################
// Function Make_Lathe_Chain
//#####################################################################
template<class TV> void RIGIDS_STANDARD_TESTS<TV>::
Make_Lathe_Chain(const FRAME<TV>& frame,const T scale,const T friction,const T cor)
{
    PHYSBAM_ASSERT(scale>0);
    ARTICULATED_RIGID_BODY<TV>& arb=rigid_body_collection.articulated_rigid_body;
    RIGID_BODY<TV>* links[7];

    for(int i=0;i<6;i++){
        RIGID_BODY<TV>& rigid_body=Add_Rigid_Body("ARB/lathe_object",scale,friction);
        rigid_body.Set_Coefficient_Of_Restitution(cor);
        links[i]=&rigid_body;
        switch(i){
            case 0:rigid_body.Frame()=frame*FRAME<TV>(TV(0,4*sin((T)pi/3),0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-(T)pi,TV(1,0,0)));break;
            case 1:rigid_body.Frame()=frame*FRAME<TV>(TV(2+2*cos((T)pi/3),4*sin((T)pi/3)*(T).5,0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-2*(T)pi/3,TV(1,0,0)));break;
            case 2:rigid_body.Frame()=frame*FRAME<TV>(TV(2+2*cos((T)pi/3),-4*sin((T)pi/3)*(T).5,0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-(T)pi/3,TV(1,0,0)));break;
            case 3:rigid_body.Frame()=frame*FRAME<TV>(TV(0,-4*sin((T)pi/3),0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0)));break;
            case 4:rigid_body.Frame()=frame*FRAME<TV>(TV(-2-2*cos((T)pi/3),-4*sin((T)pi/3)*(T).5,0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/3,TV(1,0,0)));break;
            case 5:rigid_body.Frame()=frame*FRAME<TV>(TV(-2-2*cos((T)pi/3),4*sin((T)pi/3)*(T).5,0)*scale,ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(2*(T)pi/3,TV(1,0,0)));break;}}

    links[6]=links[0];
    for(int i=0;i<6;i++){
        JOINT<TV>* joint=new POINT_JOINT<TV>();
        arb.joint_mesh.Add_Articulation(links[i]->particle_index,links[i]->particle_index,joint);
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,2*scale)));
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,0,-2*scale)));}
}
template<class TV> void RIGIDS_STANDARD_TESTS<TV>::
Set_Joint_Frames(JOINT_ID id,const TV& location)
{
    RIGID_BODY<TV>& parent=*rigid_body_collection.articulated_rigid_body.Parent(id);
    rigid_body_collection.articulated_rigid_body.joint_mesh(id)->Set_Joint_To_Parent_Frame(FRAME<TV>(-parent.Object_Space_Point(location)));
    rigid_body_collection.articulated_rigid_body.Set_Consistent_Child_Frame(id);
}
template<class TV> JOINT_ID RIGIDS_STANDARD_TESTS<TV>::
Connect_With_Point_Joint(RIGID_BODY<TV>& parent,RIGID_BODY<TV>& child,const TV& location)
{
    POINT_JOINT<TV>* joint=new POINT_JOINT<TV>;
    rigid_body_collection.articulated_rigid_body.joint_mesh.Add_Articulation(parent.particle_index,child.particle_index,joint);
    Set_Joint_Frames(joint->id_number,location);
    return joint->id_number;
}
//#####################################################################
namespace PhysBAM{
template JOINT_ID RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Connect_With_Point_Joint(RIGID_BODY<VECTOR<double,3> >&,RIGID_BODY<VECTOR<double,3> >&,VECTOR<double,3> const&);
template RIGID_BODY<VECTOR<double,1> >& RIGIDS_STANDARD_TESTS<VECTOR<double,1> >::Add_Analytic_Box(VECTOR<double,1> const&);
template RIGID_BODY<VECTOR<double,2> >& RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Analytic_Box(VECTOR<double,2> const&,int,double);
template RIGID_BODY<VECTOR<double,2> >& RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Ground(double,double,double,double);
template RIGID_BODY<VECTOR<double,2> >& RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Rigid_Body(std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,double,double,bool,bool);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Bowl(VECTOR<double,3> const&,VECTOR<double,3> const&,double,double,double,int,int);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Box(VECTOR<double,3> const&,VECTOR<int,3> const&,double);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Cylinder(double,double,int,int,double);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Shell(double,double,double,int,double);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Smooth_Gear(VECTOR<double,3> const&,int,int);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Sphere(double,double,int);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Torus(double,double,int,int,double);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Ground(double,double,double,double);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Rigid_Body(std::string const&,double,double,bool,bool);
template void RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Make_Lathe_Chain(FRAME<VECTOR<double,3> > const&,double,double,double);
template RIGID_BODY<VECTOR<double,2> >& RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Analytic_Smooth_Gear(VECTOR<double,2> const&,int,int);
template RIGID_BODY<VECTOR<double,2> >& RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Analytic_Sphere(double,double,int);
template RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::RIGIDS_STANDARD_TESTS(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,RIGID_BODY_COLLECTION<VECTOR<double,2> >&);
template JOINT_ID RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Connect_With_Point_Joint(RIGID_BODY<VECTOR<float,3> >&,RIGID_BODY<VECTOR<float,3> >&,VECTOR<float,3> const&);
template RIGID_BODY<VECTOR<float,1> >& RIGIDS_STANDARD_TESTS<VECTOR<float,1> >::Add_Analytic_Box(VECTOR<float,1> const&);
template RIGID_BODY<VECTOR<float,2> >& RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Analytic_Box(VECTOR<float,2> const&,int,float);
template RIGID_BODY<VECTOR<float,2> >& RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Ground(float,float,float,float);
template RIGID_BODY<VECTOR<float,2> >& RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Rigid_Body(std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,float,float,bool,bool);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Bowl(VECTOR<float,3> const&,VECTOR<float,3> const&,float,float,float,int,int);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Box(VECTOR<float,3> const&,VECTOR<int,3> const&,float);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Cylinder(float,float,int,int,float);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Shell(float,float,float,int,float);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Smooth_Gear(VECTOR<float,3> const&,int,int);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Sphere(float,float,int);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Torus(float,float,int,int,float);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Ground(float,float,float,float);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Rigid_Body(std::string const&,float,float,bool,bool);
template void RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Make_Lathe_Chain(FRAME<VECTOR<float,3> > const&,float,float,float);
template RIGID_BODY<VECTOR<float,2> >& RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Analytic_Smooth_Gear(VECTOR<float,2> const&,int,int);
template RIGID_BODY<VECTOR<float,2> >& RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Analytic_Sphere(float,float,int);
template RIGID_BODY<VECTOR<float,3> >& RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::Add_Analytic_Plane(VECTOR<float,3> const&,VECTOR<float,3> const&,float,int);
template RIGID_BODY<VECTOR<double,3> >& RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::Add_Analytic_Plane(VECTOR<double,3> const&,VECTOR<double,3> const&,double,int);
template RIGID_BODY<VECTOR<float,2> >& RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::Add_Analytic_Plane(VECTOR<float,2> const&,VECTOR<float,2> const&,float,int);
template RIGID_BODY<VECTOR<double,2> >& RIGIDS_STANDARD_TESTS<VECTOR<double,2> >::Add_Analytic_Plane(VECTOR<double,2> const&,VECTOR<double,2> const&,double,int);
template RIGIDS_STANDARD_TESTS<VECTOR<double,1> >::RIGIDS_STANDARD_TESTS(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,RIGID_BODY_COLLECTION<VECTOR<double,1> >&);
template RIGIDS_STANDARD_TESTS<VECTOR<double,3> >::RIGIDS_STANDARD_TESTS(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,RIGID_BODY_COLLECTION<VECTOR<double,3> >&);
template RIGIDS_STANDARD_TESTS<VECTOR<float,1> >::RIGIDS_STANDARD_TESTS(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,RIGID_BODY_COLLECTION<VECTOR<float,1> >&);
template RIGIDS_STANDARD_TESTS<VECTOR<float,2> >::RIGIDS_STANDARD_TESTS(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,RIGID_BODY_COLLECTION<VECTOR<float,2> >&);
template RIGIDS_STANDARD_TESTS<VECTOR<float,3> >::RIGIDS_STANDARD_TESTS(STREAM_TYPE,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,RIGID_BODY_COLLECTION<VECTOR<float,3> >&);
}



