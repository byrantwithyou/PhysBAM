#ifndef __RIGID_BODY_BUILDER__
#define __RIGID_BODY_BUILDER__

namespace PhysBAM
{
class POLYGONAL_SURFACE;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class RIGID_BODY_3D;
template <class T> class VECTOR_3D;
template<class T> class SYMMETRIC_MATRIX_3X3;
template<class T> class DIAGONAL_MATRIX_3X3;
template<class T> class QUATERNION;

template<class T>
class RIGID_BODY_BUILDER
{
public:
    RIGID_BODY_BUILDER();
    
    void    Set_Mass(T input_mass)
            { given_mass = true; mass = input_mass; }
    void    Set_Density(T input_density)
            { given_mass = false; density = input_density; }
    void    Set_Surface_Roughness(T input_roughness) 
            { surface_roughness = input_roughness; }
    void    Set_Write_Inside_Points(bool input_write)
            { write_inside_points = input_write; }

    static TRIANGULATED_SURFACE<T> *Create_Triangulated_Surface
                                    (const POLYGONAL_SURFACE *polygonal_surface);

    RIGID_BODY_3D<T> *Build_Rigid_Body(TRIANGULATED_SURFACE<T> *triangulated_surface, bool thin_shell=false) const;

private:
    void Calculate_Center_Of_Mass_And_Inertia_Tensor(TRIANGULATED_SURFACE<T> &triangulated_surface,
                                                    VECTOR_3D<T> &center_of_mass,
                                                    SYMMETRIC_MATRIX_3X3<T> &inertia_tensor,
                                                    T &approx_volume) const;
    void Calculate_Center_Of_Mass_And_Inertia_Tensor_Thin_Shell(TRIANGULATED_SURFACE<T> &triangulated_surface,
                                                                VECTOR_3D<T> &center_of_mass,
                                                                SYMMETRIC_MATRIX_3X3<T> &inertia_tensor,
                                                                T &approx_volume) const;
    void Transform_To_Object_Space(TRIANGULATED_SURFACE<T> &triangulated_surface,
                                   const VECTOR_3D<T> &center_of_mass,
                                   const SYMMETRIC_MATRIX_3X3<T> &inertia_tensor,
                                   QUATERNION<T> &orientation,
                                   DIAGONAL_MATRIX_3X3<T> &moments_of_inertia) const;

    bool    given_mass;
    T        mass, density;
    T        surface_roughness;
    bool    write_inside_points;
};

}

#endif
