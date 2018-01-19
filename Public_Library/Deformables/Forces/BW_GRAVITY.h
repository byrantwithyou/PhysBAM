//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_GRAVITY
//#####################################################################
#ifndef __BW_GRAVITY__
#define __BW_GRAVITY__

#include <Core/Data_Structures/FORCE_ELEMENTS.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class BW_GRAVITY:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA DEFORMABLE_FREQUENCY_DATA;

    TV gravity;
    FORCE_ELEMENTS force_particles;
public:

    BW_GRAVITY(DEFORMABLE_PARTICLES<TV>& particles_input,const TV& gravity_input)
        :DEFORMABLES_FORCES<TV>(particles_input),gravity(gravity_input)
    {}

    virtual ~BW_GRAVITY()
    {}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override
    {}

    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) override
    {}

    int Velocity_Dependent_Forces_Size() const override
    {return 0;}

    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const override
    {}

    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const override
    {}

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override
    {}

    void Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency) override
    {}

    T CFL_Strain_Rate() const override
    {return FLT_MAX;}

protected:
    template<class T_ARRAY>
    T_ARRAY Get_Particle_List(const T_ARRAY& array)
    {return array;}

public:
//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    T Potential_Energy(const T time) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
//#####################################################################
};
}
#endif
