//#####################################################################
// Copyright 2006-2008, Ron Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_BINDING
//##################################################################### 
#ifndef __PARTICLE_BINDING__
#define __PARTICLE_BINDING__

#include <Core/Arrays/ARRAY.h>
#include <Deformables/Bindings/BINDING.h>
namespace PhysBAM{

template<class TV>
class PARTICLE_BINDING:public BINDING<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef BINDING<TV> BASE;
    using BASE::particles;using BASE::particle_index;
    using BASE::One_Over_Effective_Mass; // silence -Woverloaded-virtual

    int parent;

    PARTICLE_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input)
        :BINDING<TV>(particles_input),parent(0)
    {}

    PARTICLE_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input,const int particle_index_input,const int parent_input)
        :BINDING<TV>(particles_input,particle_index_input),parent(parent_input)
    {}

    static PARTICLE_BINDING* Create(GEOMETRY_PARTICLES<TV>& particles)
    {return new PARTICLE_BINDING(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(particles));}

    virtual int Name() const override {return Static_Name();}
    static int Static_Name()
    {return 1;}

    TV Embedded_Position() const override
    {return particles.X(parent);}

    TV Embedded_Position(ARRAY_VIEW<const TV> X) const override
    {return X(parent);}

    TV Embedded_Velocity() const override
    {return particles.V(parent);}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V) const override
    {return V(parent);}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const override
    {return V(parent);}

    TV Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench) const override
    {return particles.one_over_mass(parent)*F(parent);}

    T One_Over_Effective_Mass() const override
    {return particles.one_over_mass(parent);}

    SYMMETRIC_MATRIX<T,TV::m> Impulse_Factor() const override
    {return SYMMETRIC_MATRIX<T,TV::m>()+PARTICLE_BINDING::One_Over_Effective_Mass();}

    void Apply_Impulse(const TV& impulse) override
    {PARTICLE_BINDING::Apply_Impulse(impulse,particles.V);}

    void Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V) const override
    {V(parent)+=particles.one_over_mass(parent)*impulse;}

    void Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > rigid_V) const override
    {PARTICLE_BINDING::Apply_Impulse(impulse,particles.V);}

    void Apply_Push(const TV& impulse) override
    {particles.X(parent)+=particles.one_over_mass(parent)*impulse;}

    void Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle) override
    {if(!skip_particle || !(*skip_particle)(parent)) particles.X(parent)+=dX;}

    void Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle) override
    {if(!skip_particle || !(*skip_particle)(parent)) particles.V(parent)+=dV;}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const override
    {F_full(parent)+=force;}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const override
    {F_full(parent)+=force;}

    void Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const override
    {mass_full(parent)+=particles.mass(particle_index);}

    void Parents(ARRAY<int>& parents) const override
    {parents.Append(parent);}

    void Weights(ARRAY<T>& weights) const override
    {weights.Append(1);}

private:
    void Read_Helper(TYPED_ISTREAM& input) override
    {BINDING<TV>::Read_Helper(input);Read_Binary(input,parent);}

    void Write_Helper(TYPED_OSTREAM& output) const override
    {BINDING<TV>::Write_Helper(output);Write_Binary(output,parent);}

//#####################################################################
}; 
}
#endif
