//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTIC_ETHER_DRAG
//#####################################################################
#ifndef __ELASTIC_ETHER_DRAG__
#define __ELASTIC_ETHER_DRAG__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Deformables/Forces/POINTWISE_DEFORMABLE_FORCE.h>
namespace PhysBAM{

template<class TV>
class ELASTIC_ETHER_DRAG:public POINTWISE_DEFORMABLE_FORCE<TV>
{
private:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef POINTWISE_DEFORMABLE_FORCE<TV> BASE;
    using BASE::particles;using BASE::mpi_solids;using BASE::force_particles;
    T coefficient,dt_dv_over_dx;
    T& dt;

    ELASTIC_ETHER_DRAG(DEFORMABLE_PARTICLES<TV>& particles_input,ARRAY<int>* influenced_particles_input,T coefficient,T dt_dv_over_dx,T& dt)
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,influenced_particles_input),coefficient(coefficient),dt_dv_over_dx(dt_dv_over_dx),dt(dt)
    {}

    ELASTIC_ETHER_DRAG(DEFORMABLE_PARTICLES<TV>& particles_input,const bool influence_all_particles_input,T coefficient,T dt_dv_over_dx,T& dt)
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,influence_all_particles_input),coefficient(coefficient),dt_dv_over_dx(dt_dv_over_dx),dt(dt)
    {}

    virtual ~ELASTIC_ETHER_DRAG()
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    T Potential_Energy(const T time) const override;
//#####################################################################
};
}
#endif
