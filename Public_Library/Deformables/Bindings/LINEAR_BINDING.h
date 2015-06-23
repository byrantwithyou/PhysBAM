//#####################################################################
// Copyright 2006-2008, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_BINDING
//##################################################################### 
#ifndef __LINEAR_BINDING__
#define __LINEAR_BINDING__

#include <Tools/Vectors/VECTOR.h>
#include <Deformables/Bindings/BINDING.h>
namespace PhysBAM{

template<class TV,int d>
class LINEAR_BINDING:public BINDING<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef BINDING<TV> BASE;
    using BASE::particles;using BASE::particle_index;
    using BASE::One_Over_Effective_Mass; // silence -Woverloaded-virtual

    VECTOR<int,d> parents;
    VECTOR<T,d> weights; // weights should sum to 1

    LINEAR_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input);
    LINEAR_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input,const int particle_index_input,const VECTOR<int,d>& parents_input,const VECTOR<T,d>& weights_input);
    LINEAR_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input,const int particle_index_input,const VECTOR<int,d>& parents_input,const VECTOR<T,d-1>& weights_input);

    static int Static_Name()
    {return 2+(d<<16);}

    template<class T_ARRAY_T2>
    typename T_ARRAY_T2::ELEMENT Embedded_Value(const T_ARRAY_T2& u) const
    {assert(u.Size()==particles.Size());typedef typename T_ARRAY_T2::ELEMENT T2;
    T2 result=T2();for(int i=0;i<d;i++) result+=weights[i]*u(parents[i]);
    return result;}

    static LINEAR_BINDING* Create(GEOMETRY_PARTICLES<TV>& particles);
    virtual int Name() const override;
    TV Embedded_Position() const override;
    TV Embedded_Position(ARRAY_VIEW<const TV> X) const override;
    TV Embedded_Velocity() const  override;
    TV Embedded_Velocity(ARRAY_VIEW<const TV> V) const  override;
    TV Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const  override;
    TV Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench) const override;
    T One_Over_Effective_Mass() const override;
    void Apply_Impulse(const TV& impulse) override;
    void Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V) const override;
    void Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > rigid_V) const override;
    void Apply_Push(const TV& impulse) override;
    void Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle) override;
    void Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle) override;
    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const override;
    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const override;
    void Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const override;
    void Parents(ARRAY<int>& p) const override;
    void Weights(ARRAY<T>& w) const override;
    SYMMETRIC_MATRIX<T,TV::m> Impulse_Factor() const override;
private:
    void Read_Helper(TYPED_ISTREAM& input) override;
    void Write_Helper(TYPED_OSTREAM& output) const override;
//#####################################################################
}; 
}
#endif
