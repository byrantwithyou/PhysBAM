#ifndef __BOUNDARY_CONDITIONS_BOX__
#define __BOUNDARY_CONDITIONS_BOX__
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include "BOUNDARY_CONDITIONS.h"
using namespace PhysBAM;

template<class TV>
struct BOUNDARY_CONDITIONS_BOX:public BOUNDARY_CONDITIONS<TV>
{
    typedef typename TV::SCALAR T;
    enum {source=10};
    using BOUNDARY_CONDITIONS<TV>::bounding_box;using BOUNDARY_CONDITIONS<TV>::base_domain;using BOUNDARY_CONDITIONS<TV>::grid;

    LINE_2D<T> planes[4];
    int type[4];
    TV corners[5];
    T boundary_gap;
    POLYGON<T> poly;
    mutable ARRAY<T,FACE_INDEX<TV::m> >* saved_u;

    BOUNDARY_CONDITIONS_BOX(GRID<TV>& grid_input,T offset,T angle,std::string& bc_types);
    virtual ~BOUNDARY_CONDITIONS_BOX();

    LINE_2D<T> Bounding_Edge_From_Endpoints(const TV& p1,const TV& p2) const;
    void Set_Planes_From_Corners();
    void Set_Corners_From_Box(const RANGE<TV>& box);
    void Set_Planes_From_Angle_And_Box(const RANGE<TV>& box,T angle);
    void Set_Planes_From_Offset_And_Box(const RANGE<TV>& box,T offset);
    virtual T Theta(const TV& X) const;
    virtual bool Inside(const TV& X) const;
    virtual int Boundary_Condition(const TV& X,const TV& Y,T& theta,TV& value,T time) const;
    virtual void Initialize_Velocity_Field(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const;
    virtual TV Analytic_Velocity(const TV& X,T time) const;
};

#endif


