#ifndef __BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE__
#define __BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE__
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include "BOUNDARY_CONDITIONS.h"
#include "PARAMETERS_COMMON.h"
using namespace PhysBAM;

template<class TV>
struct BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE:public BOUNDARY_CONDITIONS<TV>
{
    typedef typename TV::SCALAR T;
    enum {source=10};
    using BOUNDARY_CONDITIONS<TV>::bounding_box;using BOUNDARY_CONDITIONS<TV>::base_domain;using BOUNDARY_CONDITIONS<TV>::grid;
    using BOUNDARY_CONDITIONS<TV>::use_analytic_solution;

    BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE(GRID<TV>& grid_input);
    virtual ~BOUNDARY_CONDITIONS_GIBOU_ANALYTIC_ONE();

    virtual TV Analytic_Velocity(const TV& X,T time) const;
    virtual T Theta(const TV& X) const;
    virtual int Boundary_Condition(const TV& X,const TV& Y,T& theta,TV& value,T time) const;
    virtual void Initialize_Velocity_Field(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const;
    virtual void Update_Parameters(PARAMETERS_COMMON<T>& param);
    void Add_Initial_Error(ARRAY<T,FACE_INDEX<TV::m> >& u,T time) const;
};

#endif


