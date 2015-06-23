 //#####################################################################
// Copyright 2006, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDING
//#####################################################################
#ifndef __EMBEDDING__
#define __EMBEDDING__

#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class TV> class BINDING_LIST;

template<class TV>
class EMBEDDING:public STRUCTURE<TV>
{
    STATIC_ASSERT((TV::m>1));

    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
public:
    GEOMETRY_PARTICLES<TV>& particles;
    TRIANGLE_MESH material_surface_mesh;
    T_TRIANGULATED_OBJECT material_surface;
private:
    bool need_destroy_particles;
public:
    typedef int HAS_TYPED_READ_WRITE;

    EMBEDDING(GEOMETRY_PARTICLES<TV>& particles_input);

    virtual ~EMBEDDING()
    {if(need_destroy_particles) delete &particles;}

    static EMBEDDING<TV>* Create()
    {EMBEDDING<TV>* embedding=new EMBEDDING<TV>(*new GEOMETRY_PARTICLES<TV>);embedding->need_destroy_particles=true;return embedding;}

    static EMBEDDING<TV>* Create(GEOMETRY_PARTICLES<TV>& new_particles)
    {EMBEDDING<TV>* embedding=new EMBEDDING<TV>(new_particles);return embedding;}

    virtual std::string Name() const override {return Static_Name();}
    static std::string Static_Name()
    {return LOG::sprintf("EMBEDDING<VECTOR<T,%d> >",TV::dimension);}

    void Update_Number_Nodes() override
    {material_surface.Update_Number_Nodes();}

    void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const override
    {material_surface.Mark_Nodes_Referenced(marks,mark);}

    void Read(TYPED_ISTREAM& input) override
    {material_surface.Clean_Memory();Read_Binary(input,material_surface_mesh);}

    void Write(TYPED_OSTREAM& output) const override
    {Write_Binary(output,material_surface_mesh);}

//#####################################################################
};

template<class T_INPUT>
class EMBEDDING<VECTOR<T_INPUT,1> >:public STRUCTURE<VECTOR<T_INPUT,1> >
{
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<VECTOR<T_INPUT,1> >::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
public:
    TRIANGLE_MESH material_surface_mesh;
    T_TRIANGULATED_OBJECT material_surface;
//#####################################################################
};
}
#endif
