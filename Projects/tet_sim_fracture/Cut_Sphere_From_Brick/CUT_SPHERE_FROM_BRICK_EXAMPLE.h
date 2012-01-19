//#####################################################################
// Copyright 2003, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUT_SPHERE_FROM_BRICK_EXAMPLE
//#####################################################################
#ifndef __CUT_SPHERE_FROM_BRICK_EXAMPLE__
#define __CUT_SPHERE_FROM_BRICK_EXAMPLE__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include "../TET_SIM_EMBEDDED_EXAMPLE.h"
using namespace PhysBAM;

template <class T>
class CUT_SPHERE_FROM_BRICK_EXAMPLE:public TET_SIM_EMBEDDED_EXAMPLE<T>
{
public:
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    int m_input,n_input,mn_input;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input; 
    GRID<TV> mattress_grid;
    bool cube_mesh;
    ARRAY<int> constrained_nodes;


    CUT_SPHERE_FROM_BRICK_EXAMPLE()
                                :initial_height((T)1.2),initial_orientation(0,VECTOR_3D<T>(1,0,0)),
                                initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
                                m_input(18),n_input(18),mn_input(18),           
                                xmin_input((T)-.5),xmax_input((T)-xmin_input),ymin_input((T)-.5),ymax_input((T)-ymin_input),zmin_input((T)-.5),zmax_input((T)-zmin_input),
                                mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
                                cube_mesh(true)
    {
       restart_step_number=0;
        final_time=5;                                         
        frame_rate=24;
        cfl_number=(T).5;
        cg_iterations=250;
        youngs_modulus=.001*(7.9e7);poissons_ratio=(T).4;Rayleigh_coefficient=5*(T).02;
        strcpy(output_directory,"Cut_Sphere_From_Brick/output");
        
        perform_fracture=true;
        T scale=1;
        fracture_threshold_1=(T)1e4*scale;fracture_threshold_2=(T)5e4*scale;fracture_threshold_3=(T)8e4*scale;
        compressive_threshold_1=(T)-1e4*scale;compressive_threshold_2=(T)-5e4*scale;compressive_threshold_3=(T)-8e4*scale;

        interpolation_fraction_threshhold=(T).0002;
        max_number_of_cuts=3;
        number_of_fracture_initiation_points=1000;

        bias_stress=true;



    }

    virtual ~CUT_SPHERE_FROM_BRICK_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Tetrahedralized_Volume
//#####################################################################
void Initialize_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    tetrahedralized_volume.particles.Update_Position();
    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);

    for(int i=0;i<mattress_grid.m;i++)for(int k=0;k<mattress_grid.mn;k++){
        constrained_nodes.Append(tetrahedralized_volume.particles.array_collection->Size()+1-(1+i+k*mattress_grid.m*mattress_grid.n));}


    tetrahedralized_volume.particles.Delete_Velocity_And_Acceleration();
    tetrahedralized_volume.particles.Delete_Mass(); // in case they're accidently stored in the .tet file
    tetrahedralized_volume.particles.Update_Position_And_Velocity();
    tetrahedralized_volume.particles.Store_Mass();

    tetrahedralized_volume.tetrahedron_mesh.Initialize_Incident_Tetrahedrons();

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}
    std::cout << "total vertices = " << tetrahedralized_volume.particles.array_collection->Size() << std::endl;
    std::cout << "total tets = " << tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m << std::endl;
    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
virtual void Set_External_Velocities(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=0;i<constrained_nodes.m;i++) V(constrained_nodes(i))=VECTOR_3D<T>(0,0,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
virtual void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> ,VECTOR<int,1> >& V,const T time)
{
    for(int i=0;i<constrained_nodes.m;i++) V(constrained_nodes(i))=VECTOR_3D<T>(0,0,0);
} 
//#####################################################################
// Function Initialize_Bias_Stress_Constants
//#####################################################################
virtual void Initialize_Bias_Stress_Constants(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& etv,FRACTURE_TETRAHEDRALIZED_VOLUME<T>& rftv)
{
    TET_SIM_EMBEDDED_EXAMPLE<T>::Initialize_Bias_Stress_Constants(etv,rftv);

    etv.tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T>center=etv.tetrahedralized_volume.bounding_box->Center();
    T radius=(T).25*min(mattress_grid.xmax-mattress_grid.xmin,
                     mattress_grid.ymax-mattress_grid.ymin,
                     mattress_grid.zmax-mattress_grid.zmin);
    center.y-=(T).5*(mattress_grid.ymax - mattress_grid.ymin);
    SPHERE<T> sphere(center,radius);

    for(int t=0;t<etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;t++){
        int i,j,k,l;etv.tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.Get(t,i,j,k,l);
        T phi1=sphere.Signed_Distance(etv.tetrahedralized_volume.particles.X(i));
        T phi2=sphere.Signed_Distance(etv.tetrahedralized_volume.particles.X(j));
        T phi3=sphere.Signed_Distance(etv.tetrahedralized_volume.particles.X(k));
        T phi4=sphere.Signed_Distance(etv.tetrahedralized_volume.particles.X(l));
        int positive_count=(phi1 > 0)+(phi2 > 0)+(phi3 > 0)+(phi4 > 0);
        if(positive_count == 1 || positive_count == 2 || positive_count ==3){
            rftv.Set_Phi(t,-phi1,-phi2,-phi3,-phi4);
            rftv.fracture_bias_stress_scaling(t)=4e4;            
            //rftv.fracture_bias_direction(t)=sphere.Normal(etv.tetrahedralized_volume.Centroid(t));
            //rftv.fracture_bias_magnitude(t)=4e4;
        }
        if(positive_count == 4)rftv.fracture_bias_stress_scaling(t)=0;
        if(positive_count == 0)rftv.fracture_bias_stress_scaling(t)=0;
    }

}
//##################################################################
};
#endif
