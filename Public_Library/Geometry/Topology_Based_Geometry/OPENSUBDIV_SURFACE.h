//#####################################################################
// Copyright 2014
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENSUBDIV_SURFACE
//##################################################################### 
#ifndef __OPENSUBDIV_SURFACE__
#define __OPENSUBDIV_SURFACE__

#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{
template<class TV> class GEOMETRY_PARTICLES;

template<class TV,int gauss_order=3>
class OPENSUBDIV_SURFACE:public STRUCTURE<TV>
{
    typedef float RW;
    typedef typename TV::SCALAR T;
    typedef MATRIX<T,3> TM;
    bool need_destroy_particles;
public:
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<VECTOR<int,4> > mesh;
    ARRAY<int> control_points;
    int m;
    T thickness;
    struct FACE_DATA
    {
        typedef int HAS_TYPED_READ_WRITE;
        ARRAY<int> nodes;
        ARRAY<MATRIX<T,gauss_order> > w;
        ARRAY<MATRIX<VECTOR<T,5>,gauss_order> > A;
        MATRIX<VECTOR<TM,3>,gauss_order> G0_inv;
        MATRIX<VECTOR<T,3>,gauss_order> G0_det;
        void Read(TYPED_ISTREAM& input);
        void Write(TYPED_OSTREAM& output) const;
    };
    ARRAY<FACE_DATA> face_data; // one slot per face.
    
    OPENSUBDIV_SURFACE();
    OPENSUBDIV_SURFACE(GEOMETRY_PARTICLES<TV>& particles);
    virtual ~OPENSUBDIV_SURFACE();

//#####################################################################
    static OPENSUBDIV_SURFACE* Create(); 
    static OPENSUBDIV_SURFACE* Create(GEOMETRY_PARTICLES<TV>& particles);
    void Initialize(const std::string& filename,T thickness_in);
    TV Evaluate(T s,T t,VECTOR<TV,2>* tangents=0) const;
    OPENSUBDIV_SURFACE<TV,gauss_order>* Append_Particles_And_Create_Copy(GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) const override;
    void Read(TYPED_ISTREAM& input) override;
    void Write(TYPED_OSTREAM& output) const override;
    void Compute_G0();
    void Set_Mass(T density,bool use_constant_mass=false) const;
    std::string Name() const override;
    static std::string Static_Name();
    std::string Extension() const override;
    static std::string Static_Extension();
//#####################################################################
};
template<class TV,int gauss_order> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,2>::OBJECT*
    Create_Triangulated_Surface(const OPENSUBDIV_SURFACE<TV,gauss_order>& surf,bool same_particles);
}
#endif
