//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_VOLUME_COLLISIONS
//#####################################################################
#ifndef __LEVELSET_VOLUME_COLLISIONS__
#define __LEVELSET_VOLUME_COLLISIONS__

#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Deformables/Forces/COLLISION_FORCE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;

struct LEVELSET_VOLUME_COLLISIONS_POLYTOPE
{
    enum WORKAROUND{max_pts=14};
    VECTOR<int,max_pts> poly;
    int size;

    LEVELSET_VOLUME_COLLISIONS_POLYTOPE():size(0){}

    std::string s() const
    {
        std::string x="(";
        for(int i=0;i<size;i++){
            if(i) x+=' ';
            for(int j=0;j<8;j++)
                if(poly(i)&(1<<j))
                    x+="abcdrstu"[j];}
        return x+")";
    }
};

template<class TV>
class LEVELSET_VOLUME_COLLISIONS:public COLLISION_FORCE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<int,TV::m+1> SIMPLEX_NODES;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE SIMPLEX_FACE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT OBJECT;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT OBJECT_BOUNDARY;
    typedef VECTOR<TV,2*TV::m+2> X_VECTOR;
    typedef VECTOR<T,TV::m+1> PHI_VECTOR;
    typedef LEVELSET_VOLUME_COLLISIONS_POLYTOPE POLYTOPE;
    enum WORKAROUND{max_pts=POLYTOPE::max_pts};

    struct HYPER_PLANE
    {
        int plane;
        TV normal;
        T s;
        T Signed_Distance(const TV& location) const
        {
            return normal.Dot(location)-s;
        }
    };

public:
    using DEFORMABLES_FORCES<TV>::particles;using COLLISION_FORCE<TV>::coefficient_of_friction;
    ARRAY<OBJECT*> collision_bodies;
    ARRAY<T> undeformed_phi;
    ARRAY<ARRAY<int> > simplex_for_face;

    T stiffness;

    ARRAY<VECTOR<int,2*TV::m+2> > overlapping_particles;
    ARRAY<VECTOR<int,TV::m+1> > interior_overlapping_particles;
    T pe;
    ARRAY<VECTOR<TV,2*TV::m+2> > grad_pe;
    ARRAY<VECTOR<TV,TV::m+1> > interior_grad_pe;
    ARRAY<MATRIX<MATRIX<T,TV::m>,2*TV::m+2> > H_pe;
    ARRAY<MATRIX<MATRIX<T,TV::m>,TV::m+1> > interior_H_pe;


    LEVELSET_VOLUME_COLLISIONS(DEFORMABLE_PARTICLES<TV>& particles,T stiffness=(T)1e4);
    virtual ~LEVELSET_VOLUME_COLLISIONS();

    void Add_Mesh(OBJECT& object,const IMPLICIT_OBJECT<TV>& implicit_surface);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Update_Position_Based_State_Pair(const OBJECT& o0,OBJECT& o1);
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency) override;
    T Potential_Energy(const T time) const override;
    void Apply_Friction(ARRAY_VIEW<TV> V,const T time,const T dt) const override;
//#####################################################################
    void Simplex_Intersection(const VECTOR<TV,TV::m+1>& s,const ARRAY<HYPER_PLANE>& f,POLYTOPE& polytope) const;
    void Integrate_Simplex(const VECTOR<int,TV::m+1>& simplex,const X_VECTOR& X,const PHI_VECTOR& nodewise_undeformed_phi,VECTOR<TV,2*TV::m+2>& df,MATRIX<MATRIX<T,TV::m>,2*TV::m+2>& ddf);
};
}
#endif
