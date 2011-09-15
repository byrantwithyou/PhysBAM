#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include "../../Personal_Libraries/Eran_Library/POLYGONAL_SURFACE.h"
#include "../../Public_Library/Rigid_Bodies/MASS_PROPERTIES_3D.h"
#include "RIGID_BODY_BUILDER.h"

using namespace PhysBAM;
using namespace std;

template class RIGID_BODY_BUILDER<double>;
template class RIGID_BODY_BUILDER<float>;

//#####################################################################
// Constructor
//#####################################################################
template<class T>
RIGID_BODY_BUILDER<T>::RIGID_BODY_BUILDER()
    : given_mass(true), mass(1), density(0), surface_roughness(1e-8)
{
}

//#####################################################################
// Function Build_Rigid_Body
//#####################################################################
template<class T> RIGID_BODY_3D<T> *RIGID_BODY_BUILDER<T>::
Build_Rigid_Body(TRIANGULATED_SURFACE<T> *triangulated_surface, bool thin_shell) const
{
    cout << "Building rigid body using following parameters:" << endl;
    if (given_mass)
        cout << "  Specified mass = " << mass << endl;
    else
        cout << "  Specified density = " << density << endl;
    cout << "  Surface roughness = " << surface_roughness << endl;
    cout << endl;

    cout << "Computing center of mass and inertia tensor" << endl;
    VECTOR_3D<T> center_of_mass;
    SYMMETRIC_MATRIX_3X3<T> inertia_tensor;
    T approx_volume;
    if(thin_shell)
        Calculate_Center_Of_Mass_And_Inertia_Tensor_Thin_Shell(*triangulated_surface,center_of_mass,
                                                               inertia_tensor,approx_volume);
    else
        Calculate_Center_Of_Mass_And_Inertia_Tensor(*triangulated_surface,center_of_mass,
                                                    inertia_tensor,approx_volume);

    DIAGONAL_MATRIX_3X3<T> moments_of_inertia;
    QUATERNION<T> orientation;
    Transform_To_Object_Space(*triangulated_surface,
                              center_of_mass,
                              inertia_tensor,
                              orientation,
                              moments_of_inertia);

    cout << "  Got approx volume " << approx_volume << std::endl;
    cout << "  Got center of mass " << center_of_mass << endl;
    cout << "  Got moments of inertia " << moments_of_inertia << endl;
    cout << endl;

    // Compute mass from density if necessary
    T approx_mass = given_mass ? mass : approx_volume * density;

    RIGID_BODY_3D<T> *rigid_body = new RIGID_BODY_3D<T>();

    rigid_body->mass = approx_mass;
    rigid_body->inertia_tensor = moments_of_inertia * approx_mass; // scale by actual mass since we assumed mass=1 below
    rigid_body->surface_roughness = surface_roughness;
    rigid_body->position = center_of_mass;
    rigid_body->orientation = orientation;

    rigid_body->velocity = VECTOR_3D<T>(0,0,0);
    rigid_body->angular_momentum = VECTOR_3D<T>(0,0,0);

    rigid_body->Initialize_Triangulated_Surface(*triangulated_surface);

    return rigid_body;
}

//#####################################################################
// Function Create_Triangulated_Surface
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T> *RIGID_BODY_BUILDER<T>::
Create_Triangulated_Surface(const POLYGONAL_SURFACE *polygonal_surface)
{
    cout << "Creating triangulated surface" << endl;

    int i;

    // Initialize particles
    SOLIDS_PARTICLES<T,VECTOR_3D<T> > *particles = new SOLIDS_PARTICLES<T,VECTOR_3D<T> >();
    particles->Increase_Array_Size(polygonal_surface->vertices.m);
    for (i = 1; i <= polygonal_surface->vertices.m; i++)
    {
        int idx = particles->Add_Particle();
        particles->X(idx) = polygonal_surface->vertices(i);
    }
    
    // Initialize triangle mesh
    // Triangulate non-triangular polygons (assumes the polygons
    // are planar and doesn't try to triangulate smartly)
    TRIANGLE_MESH *triangle_mesh = new TRIANGLE_MESH();
    triangle_mesh->number_nodes = polygonal_surface->vertices.m;

    int number_triangles = 0;
    for (i = 1; i <= polygonal_surface->polygons.m; i++)
        number_triangles += polygonal_surface->polygons(i).m - 2;

    triangle_mesh->triangles.Resize(3, 1, number_triangles);
    int triidx = 0;
    for (i = 1; i <= polygonal_surface->polygons.m; i++)
    {
        for (int j = 1; j <= polygonal_surface->polygons(i).m - 2; j++)
        {
            triidx++;

            triangle_mesh->triangles(1, triidx) = polygonal_surface->polygons(i)(1);
            triangle_mesh->triangles(2, triidx) = polygonal_surface->polygons(i)(j+1);
            triangle_mesh->triangles(3, triidx) = polygonal_surface->polygons(i)(j+2);
        }
    }

    TRIANGULATED_SURFACE<T> *triangulated_surface = 
        new TRIANGULATED_SURFACE<T>(*triangle_mesh, *particles);

    cout << "  Got " << triangulated_surface->triangle_mesh.triangles.m << " triangles" << endl;

    return triangulated_surface;
}

//#####################################################################
// Function Create_Triangulated_Surface
//#####################################################################
// Assumes a uniform density and total mass 1
// center_of_mass is set to the center of mass
// inertia_tensor is set to the inertia tensor about the center of mass
template<class T> void RIGID_BODY_BUILDER<T>::
Calculate_Center_Of_Mass_And_Inertia_Tensor(TRIANGULATED_SURFACE<T> &triangulated_surface,
                                            VECTOR_3D<T> &center_of_mass,
                                            SYMMETRIC_MATRIX_3X3<T> &inertia_tensor,
                                            T &approx_volume) const
{
    assert(triangulated_surface.bounding_box);
    MASS_PROPERTIES_3D<T> mass_properties(triangulated_surface);
    mass_properties.Get_Center_Of_Mass_And_Inertia_Tensor(center_of_mass, inertia_tensor);
    approx_volume = mass_properties.Get_Volume();
}

//#####################################################################
// Function Create_Triangulated_Surface
//#####################################################################
// Assumes a uniform density and total mass 1
// center_of_mass is set to the center of mass
// inertia_tensor is set to the inertia tensor about the center of mass
template<class T> void RIGID_BODY_BUILDER<T>::
Calculate_Center_Of_Mass_And_Inertia_Tensor_Thin_Shell(TRIANGULATED_SURFACE<T> &triangulated_surface,
                                                       VECTOR_3D<T> &center_of_mass,
                                                       SYMMETRIC_MATRIX_3X3<T> &inertia_tensor,
                                                       T &approx_volume) const
{
    // compute for unit mass
    T total_mass=0;
    center_of_mass=VECTOR_3D<T>(0,0,0);
    inertia_tensor=SYMMETRIC_MATRIX_3X3<T>();
    TRIANGLE_MESH &triangle_mesh=triangulated_surface.triangle_mesh;
    for(int t=1;t<=triangle_mesh.triangles.m;t++){
        TRIANGLE_3D<T> triangle=triangulated_surface.Get_Triangle(t);
        VECTOR_3D<T> sample_point=triangle.Center();
        T area = triangle.Area();
        total_mass += area;
        center_of_mass += area * sample_point;
        inertia_tensor += area * (sample_point.Magnitude_Squared() * SYMMETRIC_MATRIX_3X3<T>::Identity_Matrix() -
                                                            SYMMETRIC_MATRIX_3X3<T>::Outer_Product(sample_point)); 
    }

    center_of_mass /= total_mass;
    inertia_tensor /= total_mass;

    // translate inertia tensor to center of mass
    inertia_tensor-=(VECTOR_3D<T>::Dot_Product(center_of_mass,center_of_mass)*SYMMETRIC_MATRIX_3X3<T>::Identity_Matrix()-SYMMETRIC_MATRIX_3X3<T>::Outer_Product(center_of_mass));

    approx_volume = 0;
}

//#####################################################################
// Function Transform_To_Object_Frame
//#####################################################################
// Given the triangulated surface with given center of mass and inertia 
// tensor we create an object frame in which the origin is the surface's 
// center of mass, and the axes are aligned with the principal axes of 
// the body.
// orientation is set to the orientation of the object frame in world space
// moments_of_inertia is also set
template<class T> void RIGID_BODY_BUILDER<T>::
Transform_To_Object_Space(TRIANGULATED_SURFACE<T> &triangulated_surface,
                          const VECTOR_3D<T> &center_of_mass,
                          const SYMMETRIC_MATRIX_3X3<T> &inertia_tensor,
                          QUATERNION<T> &orientation,
                          DIAGONAL_MATRIX_3X3<T> &moments_of_inertia) const
{
    // Figure out the orientation which gives us a diagonal inertia tensor,
    // and get the inertia_tensor
    MATRIX_3X3<T> rotation_matrix;
    DIAGONAL_MATRIX_3X3<T> diag_eigenvalues;
    inertia_tensor.Solve_Eigenproblem(diag_eigenvalues,rotation_matrix);
    moments_of_inertia=DIAGONAL_MATRIX_3X3<T>(diag_eigenvalues.x11,diag_eigenvalues.x22,diag_eigenvalues.x33);
    orientation = QUATERNION<T>(rotation_matrix);

    // Transform the particles in triangulated_surface to be relative to the new frame
    for (int i = 1; i <= triangulated_surface.particles.array_size; i++)
    {
        VECTOR_3D<T> original_pos = triangulated_surface.particles.X(i);
        VECTOR_3D<T> new_pos = orientation.Inverse_Rotate(original_pos - center_of_mass);
        triangulated_surface.particles.X(i) = new_pos;
    }

    triangulated_surface.Update_Bounding_Box();
}
