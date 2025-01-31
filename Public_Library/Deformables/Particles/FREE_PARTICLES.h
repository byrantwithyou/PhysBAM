//#####################################################################
// Copyright 2006, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//######################################################################
// Class FREE_PARTICLES
//######################################################################
#ifndef __FREE_PARTICLES__
#define __FREE_PARTICLES__

#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
namespace PhysBAM{

template<class TV>
class FREE_PARTICLES:public STRUCTURE<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_TYPED_READ_WRITE;
    ARRAY<int> nodes;

    FREE_PARTICLES();
    virtual ~FREE_PARTICLES();

    virtual std::string Name() const override {return Static_Name();}
    static std::string Static_Name()
    {return LOG::sprintf("FREE_PARTICLES<T,VECTOR<T,%d> >",TV::m);}

    static FREE_PARTICLES* Create(GEOMETRY_PARTICLES<TV>& particles)
    {return Create();}

    void Read(TYPED_ISTREAM input) override
    {Read_Binary(input,nodes);}

    void Write(TYPED_OSTREAM output) const override
    {Write_Binary(output,nodes);}

//######################################################################
    static FREE_PARTICLES* Create();
    void Mark_Nodes_Referenced(ARRAY<int>& marks,const int mark) const override;
//######################################################################
};
}
#endif
