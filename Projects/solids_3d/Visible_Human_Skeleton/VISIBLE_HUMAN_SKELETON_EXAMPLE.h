//#####################################################################
// Copyright 2004, Cynthia Lau, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VISIBLE_HUMAN_SKELETON_EXAMPLE
//##################################################################### 
#ifndef __VISIBLE_HUMAN_SKELETON_EXAMPLE__
#define __VISIBLE_HUMAN_SKELETON_EXAMPLE__

#include <PhysBAM_Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <Rigid_Bodies/MASS_PROPERTIES_3D.h>

namespace PhysBAM{

template<class T,class RW>
class VISIBLE_HUMAN_SKELETON_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    T time1,time2,time3,time4,time5;
    T kinematic_body_vel;

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

    VISIBLE_HUMAN_SKELETON_EXAMPLE():SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE)
    {
        last_frame=600;
        output_directory="Visible_Human_Skeleton/output";
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
    time1=(T)10/24;
    time2=(T)200/24;
    time3=(T)350/24;
    time4=(T)400/24;

    RIGID_BODY<TV> *rigid_body=0;
    T scale=1;
    T mass=1;
    
    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    rigid_body->position.y -= 2;

    //sternum
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/sternum_ribs", scale);
    rigid_body->Set_Mass(mass);
    //clavicle
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/clavicle_right", scale);
    rigid_body->Set_Mass(mass);
    //scapula
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/scapula_right", scale);
    rigid_body->Set_Mass(mass);
    //humerus
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/humerus_right", scale);
    rigid_body->Set_Mass(mass);
    //ulna
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/ulna_right", scale);
    rigid_body->Set_Mass(mass);
    //radius
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/radius_right", scale);
    rigid_body->Set_Mass(mass);

    //reflected clavicle
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/clavicle_right_reflected", scale);
    rigid_body->Set_Mass(mass);
    //reflected scapula
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/scapula_right_reflected", scale);
    rigid_body->Set_Mass(mass);
    //reflected humerus
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/humerus_right_reflected", scale);
    rigid_body->Set_Mass(mass);
    //reflected ulna
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/ulna_right_reflected", scale);
    rigid_body->Set_Mass(mass);
    //reflected radius
    rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/radius_right_reflected", scale);
    rigid_body->Set_Mass(mass);
/*
    //sternum
    rigid_body=Create_Rigid_Body("/VH_Bones/sternum_ribs", scale, T(0.005));
    //rigid_body=Initialize_Rigid_Body("Visible_Human_Bones/sternum_and_ribs_right", scale);
    rigid_body->Set_Mass(mass);
    //clavicle
    rigid_body=Create_Rigid_Body("/VH_Bones/clavicle_right", scale, T(0.007));
    rigid_body->Set_Mass(mass);
    //scapula
    rigid_body=Create_Rigid_Body("/VH_Bones/scapula_right", scale, T(0.006));
    rigid_body->Set_Mass(mass);
    //humerus
    rigid_body=Create_Rigid_Body("/VH_Bones/humerus_right", scale, T(0.006));
    rigid_body->Set_Mass(mass);
    //ulna
    rigid_body=Create_Rigid_Body("/VH_Bones/ulna_right", scale, T(0.006));
    rigid_body->Set_Mass(mass);
    //radius
    rigid_body=Create_Rigid_Body("/VH_Bones/radius_right", scale, T(0.006));
    rigid_body->Set_Mass(mass);

    //reflected clavicle
    rigid_body=Create_Rigid_Body("/VH_Bones/clavicle_right_reflected", scale, T(0.007));
    rigid_body->Set_Mass(mass);
    //reflected scapula
    rigid_body=Create_Rigid_Body("/VH_Bones/scapula_right_reflected", scale, T(0.006));
    rigid_body->Set_Mass(mass);
    //reflected humerus
    rigid_body=Create_Rigid_Body("/VH_Bones/humerus_right_reflected", scale, T(0.006));
    rigid_body->Set_Mass(mass);
    //reflected ulna
    rigid_body=Create_Rigid_Body("/VH_Bones/ulna_right_reflected", scale, T(0.006));
    rigid_body->Set_Mass(mass);
    //reflected radius
    rigid_body=Create_Rigid_Body("/VH_Bones/radius_right_reflected", scale, T(0.006));
    rigid_body->Set_Mass(mass);
*/
    kinematic_body_vel=T(5);
}

//#####################################################################
// Function Initialize_Rigid_Body
//#####################################################################
RIGID_BODY<TV>* Initialize_Rigid_Body(const std::string& basename,T scaling_factor=1)
{
    std::string fullbasename=data_directory+"/Rigid_Bodies/"+basename;
    int index=solid_body_collection.deformable_object.rigid_body_particles.template Add_Rigid_Body<T>(fullbasename, scaling_factor);

    RIGID_BODY<TV> *rigid_body=&solid_body_collection.deformable_object.rigid_body_particles.Rigid_Body(index);

    // If this is a ground plane (based on filename which is kind of a hack) then we fix the bounding box here to only be extended 
    // downwards (rather than ymin-=1, ymax+=1 as is done in RIGID_BODIES_DRIVER).  This reduces bounding volume collisions between 
    // objects near the plane and the plane.
    if(basename=="ground"){
        BOX_3D<T>& box=*rigid_body->triangulated_surface->bounding_box;
        if(box.ymin == box.ymax){ 
            std::cout<<"Fixing ground bounding box"<<std::endl;
            std::cout<<"before: "<<box<<std::endl;
            box.ymin-=2;
            std::cout<<"after: "<<box<<std::endl;}}

    rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    return rigid_body;
}
/*
//#####################################################################
// Function Transform_Rigid_Body
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Transform_Rigid_Body(RIGID_BODY<TV>* rigid_body, MATRIX<T,4>& left_transform)
{
    //transform triangulated surface
    TRIANGULATED_SURFACE<T> *surface=rigid_body->triangulated_surface;
    for(int i=0;i<surface->particles.array_size;i++)
    {
        surface->particles.X(i)=left_transform * surface->particles.X(i);
    }
    int temp;
    for(int j=0;j<surface->triangle_mesh.triangles.m;j++)
    {
        temp=surface->triangle_mesh.triangles(1,j);
        surface->triangle_mesh.triangles(1,j)=surface->triangle_mesh.triangles(3,j);
        surface->triangle_mesh.triangles(3,j)=temp;
    }
    //transform level set
    //transform rigid body
}

*/
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    T baseboxsize=2;
    if(time<time1){frame.t=VECTOR_3D<T>(0,baseboxsize,0);frame.r=QUATERNION<T>();}
    else if(time<time2){frame.t=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time-time1),0);frame.r=QUATERNION<T>();}
    else if(time<time3){frame.t=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);frame.r=QUATERNION<T>();}
    else if(time<time4){frame.t=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);frame.r=QUATERNION<T>(-(T)0.4*(time-time3),VECTOR_3D<T>(0,0,1));}
    else{frame.t=VECTOR_3D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1),0);frame.r=QUATERNION<T>(-(T)0.4*(time4-time3),VECTOR_3D<T>(0,0,1));}
}



//#####################################################################
// Function Create_Rigid_Body
//#####################################################################
RIGID_BODY<TV> *Create_Rigid_Body(std::string basename, T scaling_factor, T bone_levelset_dx_percentage)
{
    std::fstream input;std::string filename,fullbasename,filename2;
    RIGID_BODY<TV> *rigid_body=new RIGID_BODY<TV>();
    VECTOR_3D<T> local_coordinate,temp;
    BOX_3D<T> box;T largest_dimension,dx;int m,n,mn,index;
    SYMMETRIC_MATRIX<T,3> inertia_tensor;
    FRAME<T> frame;
    DIAGONAL_MATRIX<T,3> eigenvalues;
    MATRIX<T,3> eigenvectors;

    fullbasename=data_directory+basename;
    // triangulated surface
    index=rg_list->triangulated_surface_list.Add_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>* tri_surface=rg_list->triangulated_surface_list.triangulated_surface(index);
    filename=fullbasename+".tri";input.open(filename.c_str(),std::ios::in|std::ios::binary);
    rg_list->triangulated_surface_list.triangulated_surface(index)->template Read<RW>(input);input.close();
    rg_list->triangulated_surface_hash.Set(filename,index);
    MASS_PROPERTIES_3D<T> mass_props(*(rg_list->triangulated_surface_list.triangulated_surface(index)));
    mass_props.Get_Center_Of_Mass_And_Inertia_Tensor(frame.t,inertia_tensor);
    inertia_tensor.Solve_Eigenproblem(eigenvalues, eigenvectors);
    frame.r=QUATERNION<T>(eigenvectors);
    // compute bone local coordinates
    tri_surface->Clean_Memory();
    for(int p=0;p<tri_surface->particles.array_collection->Size();p++){
        local_coordinate=frame.Local_Coordinate(tri_surface->particles.X(p));
        tri_surface->particles.X(p)=local_coordinate;}
    MASS_PROPERTIES_3D<T> mass_props2(*tri_surface);
    mass_props2.Get_Center_Of_Mass_And_Inertia_Tensor(temp,inertia_tensor);
    // rigid body
    filename=data_directory+"/Rigid_Bodies/ground.rgd";input.open(filename.c_str(),std::ios::in|std::ios::binary);//dummy file
    rigid_body->template Read<RW>(input);input.close();
    rigid_body->position=frame.t;
    rigid_body->orientation=frame.r;
    rigid_body->inertia_tensor=DIAGONAL_MATRIX<T,3>(inertia_tensor.x11,inertia_tensor.x22,inertia_tensor.x33);
    rigid_body->coefficient_of_friction=(T).3;
    rigid_body->Initialize_Triangulated_Surface(*tri_surface);
    // implicit surface
    tri_surface->Update_Bounding_Box();
    box=*(tri_surface->bounding_box);
    largest_dimension=max(box.xmax-box.xmin,box.ymax-box.ymin,box.zmax-box.zmin);dx=bone_levelset_dx_percentage*largest_dimension;
    m=int((box.xmax-box.xmin)/dx);n=int((box.ymax-box.ymin)/dx);mn=int((box.zmax-box.zmin)/dx);
    T extra_padding=(T)10*dx;
    LEVELSET_IMPLICIT_SURFACE<T>* implicit_surface=LEVELSET_IMPLICIT_SURFACE<T>::Create();
    implicit_surface->levelset.phi.Resize(1,m,1,n,1,mn);
    std::cout<<"levelset "<<basename<<" "<<m<<" "<<n<<" "<<mn<<std::endl;
    implicit_surface->levelset.grid.Initialize(m,n,mn,box.xmin-extra_padding,box.xmax+extra_padding,box.ymin-extra_padding,box.ymax+extra_padding,box.zmin-extra_padding,box.zmax+extra_padding);
    SIGNED_DISTANCE::Calculate(*tri_surface,implicit_surface->levelset.grid,implicit_surface->levelset.phi,true);
    implicit_surface->Update_Box();
    filename=fullbasename+".phi";
    index=rg_list->implicit_surface_list.Add_Implicit_Surface(implicit_surface);
    rg_list->implicit_surface_hash.Set(filename,index);
    rigid_body->Initialize_Implicit_Surface(*implicit_surface);
    rg_list->Add_Rigid_Body(rigid_body);
    rigid_body->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    rigid_body->Update_Angular_Velocity();
    std::fstream output;
    filename2=output_directory+basename.substr(basename.rfind("/", basename.length()),basename.length())+".tri";
    output.open(filename2.c_str(),std::ios::out|std::ios::binary);
    tri_surface->template Write<RW>(output);output.close();
    filename2=output_directory+basename.substr(basename.rfind("/", basename.length()),basename.length())+".phi";
    output.open(filename2.c_str(),std::ios::out|std::ios::binary);
    implicit_surface->template Write<RW>(output);output.close();
    filename2=output_directory+basename.substr(basename.rfind("/", basename.length()),basename.length())+".rgd";
    output.open(filename2.c_str(),std::ios::out|std::ios::binary);
    rigid_body->template Write<RW>(output);output.close();
    return rigid_body;
}

//#####################################################################
};
}
#endif
