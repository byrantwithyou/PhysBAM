//#####################################################################
// Copyright 2006-2008, Ron Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_BINDING
//##################################################################### 
#ifndef __PARTICLE_BINDING__
#define __PARTICLE_BINDING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING.h>
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

    virtual int Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static int Static_Name()
    {return 1;}

    TV Embedded_Position() const PHYSBAM_OVERRIDE
    {return particles.X(parent);}

    TV Embedded_Position(ARRAY_VIEW<const TV> X) const PHYSBAM_OVERRIDE
    {return X(parent);}

    TV Embedded_Velocity() const PHYSBAM_OVERRIDE
    {return particles.V(parent);}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V) const PHYSBAM_OVERRIDE
    {return V(parent);}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const PHYSBAM_OVERRIDE
    {return V(parent);}

    TV Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench) const PHYSBAM_OVERRIDE
    {return particles.one_over_mass(parent)*F(parent);}

    T One_Over_Effective_Mass() const PHYSBAM_OVERRIDE
    {return particles.one_over_mass(parent);}

    SYMMETRIC_MATRIX<T,TV::m> Impulse_Factor() const PHYSBAM_OVERRIDE
    {return SYMMETRIC_MATRIX<T,TV::m>()+PARTICLE_BINDING::One_Over_Effective_Mass();}

    void Apply_Impulse(const TV& impulse) PHYSBAM_OVERRIDE
    {PARTICLE_BINDING::Apply_Impulse(impulse,particles.V);}

    void Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V) const PHYSBAM_OVERRIDE
    {V(parent)+=particles.one_over_mass(parent)*impulse;}

    void Apply_Push(const TV& impulse) PHYSBAM_OVERRIDE
    {particles.X(parent)+=particles.one_over_mass(parent)*impulse;}

    void Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE
    {if(!skip_particle || !(*skip_particle)(parent)) particles.X(parent)+=dX;}

    void Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE
    {if(!skip_particle || !(*skip_particle)(parent)) particles.V(parent)+=dV;}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const PHYSBAM_OVERRIDE
    {F_full(parent)+=force;}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const PHYSBAM_OVERRIDE
    {F_full(parent)+=force;}

    void Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const PHYSBAM_OVERRIDE
    {mass_full(parent)+=particles.mass(particle_index);}

    void Parents(ARRAY<int>& parents) const PHYSBAM_OVERRIDE
    {parents.Append(parent);}

    void Weights(ARRAY<T>& weights) const PHYSBAM_OVERRIDE
    {weights.Append(1);}

private:
    void Read_Helper(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE
    {BINDING<TV>::Read_Helper(input);Read_Binary(input,parent);}

    void Write_Helper(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE
    {BINDING<TV>::Write_Helper(output);Write_Binary(output,parent);}

//#####################################################################
}; 
}
#endif
