//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonthan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_FORCE_COLLECTION
//#####################################################################
#ifndef __DEFORMABLE_FORCE_COLLECTION__
#define __DEFORMABLE_FORCE_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class DEFORMABLE_FORCE_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_DEFORMABLE;
public:
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    ARRAY<DEFORMABLES_FORCES<TV>*> deformables_forces;

    ARRAY<T_FREQUENCY_DEFORMABLE> frequency; // hertz for deformable CFL
    T cfl_number;
    T cfl_elastic,cfl_damping;
    bool implicit_damping;

    bool print_energy;

    DEFORMABLE_FORCE_COLLECTION(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection);
    virtual ~DEFORMABLE_FORCE_COLLECTION();

    template<class T_FORCE> T_FORCE
    Find_Force(const int index=0)
    {return Find_Type<T_FORCE>(deformables_forces,index);}

    template<class T_FORCE> const T_FORCE
    Find_Force(const int index=0) const
    {return Find_Type<T_FORCE>(deformables_forces,index);}

//#####################################################################
    void Set_CFL_Number(const T cfl_number_input=.5);
    int Add_Force(DEFORMABLES_FORCES<TV>* force);
    void Update_CFL();
    T CFL(const bool verbose=false);
    T CFL_Elastic_And_Damping();
    T CFL_Elastic();
    T CFL_Damping();
    T CFL_Strain_Rate();

    void Update_Position_Based_State(const T time,const bool is_position_update);
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F_full,const T time) const;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T time) const;
    void Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T scale,const T time) const;

    void Test_Energy(const T time);
    void Test_Force_Derivatives(const T time);
//#####################################################################
};
}
#endif
