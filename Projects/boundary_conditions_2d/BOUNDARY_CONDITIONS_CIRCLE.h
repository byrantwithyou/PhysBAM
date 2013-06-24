#ifndef __BOUNDARY_CONDITIONS_CIRCLE__
#define __BOUNDARY_CONDITIONS_CIRCLE__
#include <Geometry/Basic_Geometry/SPHERE.h>
#include "BOUNDARY_CONDITIONS.h"
using namespace PhysBAM;

template<class TV>
struct BOUNDARY_CONDITIONS_CIRCLE:public BOUNDARY_CONDITIONS<TV>
{
    typedef typename TV::SCALAR T;
    enum {source=10};
    using BOUNDARY_CONDITIONS<TV>::bounding_box;using BOUNDARY_CONDITIONS<TV>::base_domain;using BOUNDARY_CONDITIONS<TV>::grid;

    SPHERE<TV> sphere;
    int type;
    mutable ARRAY<T,FACE_INDEX<TV::m> >* saved_u;

    BOUNDARY_CONDITIONS_CIRCLE(GRID<TV>& grid_input,int bc_type);
    virtual ~BOUNDARY_CONDITIONS_CIRCLE();

    virtual T Theta(const TV& X) const;
    virtual bool Inside(const TV& X) const;
    virtual int Boundary_Condition(const TV& X,const TV& Y,T& theta,TV& value,T time) const;
    virtual void Initialize_Velocity_Field(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const;
    virtual TV Analytic_Velocity(const TV& X,T time) const;
    void Add_Initial_Error(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const;
};

#endif


