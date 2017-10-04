//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/Robust_Arithmetic.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_PARTICLES<TV>::
DEFORMABLE_PARTICLES()
    :mass(0,0),one_over_mass(0,0),effective_mass(0,0),one_over_effective_mass(0,0),store_mass(false)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_PARTICLES<TV>::
~DEFORMABLE_PARTICLES()
{}
//#####################################################################
// Function Center_Of_Mass
//#####################################################################
template<class TV> TV DEFORMABLE_PARTICLES<TV>::
Center_Of_Mass() const
{
    if(this->store_mass){
        T total=mass.Sum();
        return total?X.Weighted_Sum(mass)/total:TV();}
    return X.Average(); // default to treating mass as one
}
//#####################################################################
// Function Compute_Auxiliary_Attributes
//#####################################################################
template<class TV> void DEFORMABLE_PARTICLES<TV>::
Compute_Auxiliary_Attributes(const SOFT_BINDINGS<TV>& soft_bindings)
{
    Compute_Auxiliary_Attributes(soft_bindings,IDENTITY_ARRAY<>(soft_bindings.particles.Size()),false);
}
//#####################################################################
// Function Compute_Auxiliary_Attributes
//#####################################################################
template<class TV> template<class T_INDICES> void DEFORMABLE_PARTICLES<TV>::
Compute_Auxiliary_Attributes(const SOFT_BINDINGS<TV>& soft_bindings,const T_INDICES& indices,const bool copy_existing_elements)
{
    if(this->template Get_Array<T>("one_over_mass"))
        Remove_Array_Using_Index(Get_Attribute_Index("one_over_mass"));
    this->template Add_Array<T>("one_over_mass",&one_over_mass);
    if(this->template Get_Array<T>("effective_mass"))
        Remove_Array_Using_Index(Get_Attribute_Index("effective_mass"));
    this->template Add_Array<T>("effective_mass",&effective_mass);
    if(this->template Get_Array<T>("ove_over_effective_mass"))
        Remove_Array_Using_Index(Get_Attribute_Index("ove_over_effective_mass"));
    this->template Add_Array<T>("ove_over_effective_mass",&one_over_effective_mass);
    for(int i=0;i<indices.Size();i++){int p=indices(i);
        one_over_mass(p)=Robust_Inverse(mass(p));}
    for(int i=0;i<indices.Size();i++){int p=indices(i);
        one_over_effective_mass(p)=soft_bindings.One_Over_Effective_Mass(p);
        effective_mass(p)=Robust_Inverse(one_over_effective_mass(p));}
}
//#####################################################################
template class DEFORMABLE_PARTICLES<VECTOR<float,1> >;
template class DEFORMABLE_PARTICLES<VECTOR<float,2> >;
template class DEFORMABLE_PARTICLES<VECTOR<float,3> >;
template class DEFORMABLE_PARTICLES<VECTOR<double,1> >;
template class DEFORMABLE_PARTICLES<VECTOR<double,2> >;
template class DEFORMABLE_PARTICLES<VECTOR<double,3> >;
}
