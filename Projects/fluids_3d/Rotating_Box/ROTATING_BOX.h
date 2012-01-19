//#####################################################################
// Class ROTATING_BOX 
//##################################################################### 
//
//#####################################################################
// Fedkiw - February 14, 2002
// Koltakov - February 26, 2003
//#####################################################################
#ifndef __ROTATING_BOX__
#define __ROTATING_BOX__    

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_BOX.h>
#include "../SMOKE_3D_EXAMPLE.h"
#include "STRING_INT_MAP.h"
#include <Geometry/IMPLICIT_SURFACE.h>

using namespace PhysBAM;

class ROTATING_BOX : public SMOKE_3D_EXAMPLE
{
public:    
    ROTATING_BOX()
    {
        time = 0; final_time = 1;
        m=100,n=100,mn=100;
        xmin=-1.5;xmax=2.5;ymin=-1.5;ymax=2.5;zmin=-1.5;zmax=2.5;
        Initialize();
    }

    ~ROTATING_BOX()
    {}

//#####################################################################
// Function Initialize_Main_Program
//#####################################################################
void Initialize()
{            

    RIGID_BODY<TV> *rigid_body = 0;

    // SPHERE
    rigid_body=Initialize_Rigid_Body("sphere");
    rigid_body->position=VECTOR_3D(0,1,.5);
    rigid_body->velocity=VECTOR_3D(10,0,0);
//    rigid_body->angular_velocity=VECTOR_3D(1,0,0);
//    rigid_body->Update_Angular_Velocity();
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("sphere (coeff 1)");
    rigid_bodies.Append(rigid_body);
}                    

//#####################################################################
// Function Update_Geometry
//#####################################################################
void Update_Geometry_Old()
{    
    int i,j,ij;

    // source
    if(time == 0){ 
        ELLIPSOID ellipsoid(VECTOR_3D(.5,.5,.5),VECTOR_3D(.4,.2,.2));
        for(i=0;i<m;i++) for(j=0;j<n;j++) for(ij=0;ij<mn;ij++)
            if(ellipsoid.Inside(VECTOR_3D(grid->x(i),grid->y(j),grid->z(ij)))){(*density)(i,j,ij)=1;(*T)(i,j,ij)=1;}}

//    CYLINDER cylinder();
//    cylinder->Set_Radius(.1);cylinder->Set_Height(.1);
//    double x_location=.1+2.4*time/9;
//    cylinder->Transform(Translation_Matrix(VECTOR_3D(x_location,.1,(zmax+zmin)/2)));
//    for(i=0;i<m;i++) for(j=0;j<n;j++) for(ij=0;ij<mn;ij++)
//        if(cylinder->Inside(VECTOR_3D(grid->x(i),grid->y(j),grid->z(ij)))){
//            (*density)(i,j,ij)=1;(*T)(i,j,ij)=1;(*psi_N)(i,j,ij)=1;/*(*u_fixed)(i,j,ij)=.2;*/(*v_fixed)(i,j,ij)=.8;}
        
    // box
    RENDERING_BOX box(.2,.8,.45,.55,.45,.55);
    double speed=2*pi,rotation=speed*time;
    box.transform.MATRIX_4X4::Translation_Matrix(VECTOR_3D(.5,0,.5))*box.transform.MATRIX_4X4::Rotation_Matrix_Y_Axis(rotation)*box.transform.MATRIX_4X4::Translation_Matrix(VECTOR_3D(-.5,0,-.5));
    for(i=0;i<m;i++) for(j=0;j<n;j++) for(ij=0;ij<mn;ij++)
        if(box.Inside(VECTOR_3D(grid->x(i),grid->y(j),grid->z(ij)))){
            (*psi_N)(i,j,ij)=1;
            (*u_fixed)(i,j,ij)=speed*(grid->z(ij)-.5);
            (*v_fixed)(i,j,ij)=0;
            (*w_fixed)(i,j,ij)=-speed*(grid->x(i)-.5);}    
}
//#####################################################################

void Update_Geometry()//ARRAY<RIGID_BODY<TV>*>& rigid_bodies)
{
    for(int index=0;index<rigid_bodies.m;index++) {
        for(int i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++)
        if(rigid_bodies(index)->Implicit_Surface_Inside(VECTOR_3D(grid->x(i),grid->y(j),grid->z(ij)))){
            (*density)(i,j,ij)=1;(*T)(i,j,ij)=1;
            (*psi_N)(i,j,ij)=1;
            VECTOR_3D point_velocity = rigid_bodies(index)->Pointwise_Object_Velocity(VECTOR_3D(grid->x(i),grid->y(j),grid->z(ij)));
            (*u_fixed)(i,j,ij)= point_velocity.x;
            (*v_fixed)(i,j,ij)= point_velocity.y;
            (*w_fixed)(i,j,ij)= point_velocity.z;}
    }
}

private:
//#####################################################################
// Function Initialize_Rigid_Body
//#####################################################################
// This function hashes triangulated surface and implicit surface filenames so the same file doesn't get loaded more than once.
// It now takes into account rescaling when hashing the surfaces.
RIGID_BODY<TV>* Initialize_Rigid_Body(const char* basename,double scaling_factor=1)
{
    char filename[256],hashname[256], rigid_bodies_repository[256];int index;
    bool verbose = true;
    STRING_INT_MAP triangulated_surface_hash,implicit_surface_hash;
    strcpy(rigid_bodies_repository,"../../Public_Library/Data/Rigid_Bodies");
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>();

    // triangulated surface
    sprintf(filename,"%s/%s.tri",rigid_bodies_repository,basename);
    TRIANGULATED_SURFACE* triangulated_surface=0;
    if(scaling_factor != 1) sprintf(hashname,"%s@%.6f",filename,scaling_factor); // mangle hash name if we're rescaling it
    else strcpy(hashname,filename);
    if(triangulated_surface_hash.Find(hashname,index)){ // already read in triangulated surface
        if(verbose) std::cout << "Using hashed " << hashname << std::endl;
        triangulated_surface_list.Add_Triangulated_Surface(index);
        triangulated_surface=triangulated_surface_list.triangulated_surface(index);}
    else{ // read in triangulated surface for the first time
        std::ifstream input(filename,std::ios::binary);if(!input){std::cerr << "Error loading file " << filename << std::endl;exit(-1);}
        index=triangulated_surface_list.Add_Triangulated_Surface();
        triangulated_surface=triangulated_surface_list.triangulated_surface(index);
        triangulated_surface->Read(input);
        if(scaling_factor != 1) triangulated_surface->Rescale(scaling_factor);
        triangulated_surface_hash.Add_Pair(hashname,index);}
    rigid_body->Initialize_Triangulated_Surface(*triangulated_surface);

    // implicit surface
    sprintf(filename,"%s/%s.phi",rigid_bodies_repository,basename);
    IMPLICIT_SURFACE* implicit_surface=0;
    if(scaling_factor != 1) sprintf(hashname,"%s@%.6f",filename,scaling_factor); // mangle hash name if we're rescaling it
    else strcpy(hashname,filename);
    if(implicit_surface_hash.Find(hashname,index)){    // already read in implicit surface
        if(verbose) std::cout << "Using hashed " << hashname << std::endl;
        implicit_surface_list.Add_Levelset_Implicit_Surface(index);
        implicit_surface=implicit_surface_list.implicit_surface(index);}
    else{ // read in implicit surface for the first time
        std::ifstream input(filename, std::ios::binary);if(!input){std::cerr << "Error loading file " << filename << std::endl;exit(-1);}
        index=implicit_surface_list.Add_Levelset_Implicit_Surface();
        implicit_surface=implicit_surface_list.implicit_surface(index);
        implicit_surface->Read(input);
        if(scaling_factor != 1) implicit_surface->Rescale(scaling_factor);
        implicit_surface_hash.Add_Pair(hashname,index);}
    rigid_body->Initialize_Implicit_Surface(*implicit_surface);

    // rigid body
    sprintf(filename, "%s/%s.rgd", rigid_bodies_repository,basename);
    std::ifstream input(filename,std::ios::binary);if(!input){std::cerr << "Error loading file " << filename << std::endl;exit(-1);}
    rigid_body->Read(input);
    if(scaling_factor != 1) rigid_body->Rescale(scaling_factor);

    // If this is a ground plane (based on filename which is kind of a hack) then we fix the bounding box here to only be extended 
    // downwards (rather than ymin-=1, ymax+=1 as is done in RIGID_BODIES_DRIVER).  This reduces bounding volume collisions between 
    // objects near the plane and the plane.
    if(!strcmp(basename,"ground")){
        BOX_3D& box=*rigid_body->triangulated_surface->bounding_box;
        if(box.ymin == box.ymax){ 
            std::cout << "Fixing ground bounding box" << std::endl;
            std::cout << "before: " << box << std::endl;
            box.ymin-=2;
            std::cout << "after: " << box << std::endl;}}

    return rigid_body;
}
};      
#endif
