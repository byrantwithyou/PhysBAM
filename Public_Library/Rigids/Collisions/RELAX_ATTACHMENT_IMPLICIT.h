//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RELAX_ATTACHMENT_IMPLICIT
//#####################################################################
#ifndef __RELAX_ATTACHMENT_IMPLICIT__
#define __RELAX_ATTACHMENT_IMPLICIT__

#include <Core/Matrices/MATRIX.h>
namespace PhysBAM{
template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class RELAX_ATTACHMENT_IMPLICIT
{
    typedef typename TV::SCALAR T;
public:
    TV K;
    MATRIX<T,TV::m> dKdZ,dKdX,dKdW,dKdN;
    TV dKdphi;
    bool dynamic;

    void Relax(const TV& Z,const TV& X,const TV& W,T mu);
    void Relax_Search(const TV& Z,const TV& X,const TV& W,const IMPLICIT_OBJECT<TV>* io,T mu);
};
}
#endif
