//#####################################################################
// Copyright 2016, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_GRID_OBJECT
//##################################################################### 
#ifndef __OPENGL_GRID_OBJECT__
#define __OPENGL_GRID_OBJECT__

#include <OpenGL/OpenGL/OPENGL_OBJECT.h>

namespace PhysBAM
{

template<class TV>
class OPENGL_GRID_OBJECT
{
public:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    virtual void Print_Cell_Selection_Info(std::ostream& stream,const TV_INT& cell) const = 0;
    virtual void Print_Node_Selection_Info(std::ostream& stream,const TV_INT& node) const = 0;
};

}

#endif
