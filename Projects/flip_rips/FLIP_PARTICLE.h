//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FLIP_PARTICLE__
#define __FLIP_PARTICLE__

#include <Tools/Vectors/VECTOR.h>

namespace PhysBAM{

template<class TV,int w>
class FLIP_PARTICLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef FACE_INDEX<TV::m> T_FACE;
    typedef GRID<TV> T_GRID;
    
public:

    TV X,V;

    FLIP_PARTICLE(){}
    ~FLIP_PARTICLE(){}

    VECTOR<TV_INT,TV::m+1> base_index;
    VECTOR<VECTOR<VECTOR<T,TV::m>,w>,TV::m+1> weights;
    VECTOR<VECTOR<VECTOR<T,TV::m>,w>,TV::m+1> dweight;

    enum WORKAROUND{ghost=3};

    void Update_Base_And_Weights(const T_GRID& grid)
    {
        for(int j=0;j<TV::m+1;j++){
            const TV offset=(j==TV::m)?TV():TV::Axis_Vector(j);
            base_index(j)=grid.Cell(X+(T).5*(offset+1-w)*grid.dX,ghost);
            const T_FACE face=T_FACE(j,base_index(j));
            const TV X_eval=X-((j==TV::m)?grid.Center(base_index(j)):grid.Face(face));
            for(int k=0;k<TV::m;k++){
                const T x=X_eval(k);
                const T one_over_dx=(T)1/grid.DX()(k);
                Compute_Weights(x*one_over_dx,one_over_dx,k,j);}}
    }

    void Compute_Weights(T x,const T one_over_dx,const int k,const int j)
    {
        switch(w){
            case 2:{
                weights(j)(0)(k)=(T)1-x;
                dweight(j)(0)(k)=-one_over_dx;
                x-=(T)1;
                weights(j)(1)(k)=(T)1+x;
                dweight(j)(1)(k)=one_over_dx;
            } break;
            case 3:{
                weights(j)(0)(k)=(T).5*x*x-(T)1.5*x+(T)1.125;
                dweight(j)(0)(k)=(x-(T)1.5)*one_over_dx;
                x-=(T)1;
                weights(j)(1)(k)=-x*x+(T).75;
                dweight(j)(1)(k)=(T)(-2)*x*one_over_dx;
                x-=(T)1;
                weights(j)(2)(k)=(T).5*x*x+(T)1.5*x+(T)1.125;
                dweight(j)(2)(k)=(x+(T)1.5)*one_over_dx;
            } break;
            case 4:{
                T x2,x3;
                x2=x*x;x3=x2*x;
                weights(j)(0)(k)=-((T)1/6)*x3+x2-(T)2*x+(T)4/3;
                dweight(j)(0)(k)=(-(T).5*x2+(T)2*x-(T)2)*one_over_dx;
                x-=(T)1;
                x2=x*x;x3=x2*x;
                weights(j)(1)(k)=(T).5*x3-x2+(T)2/3;
                dweight(j)(1)(k)=((T)1.5*x2-(T)2*x)*one_over_dx;
                x-=(T)1;
                x2=x*x;x3=x2*x;
                weights(j)(2)(k)=-(T).5*x3-x2+(T)2/3;
                dweight(j)(2)(k)=(-(T)1.5*x2-(T)2*x)*one_over_dx;
                x-=(T)1;
                x2=x*x;x3=x2*x;
                weights(j)(3)(k)=((T)1/6)*x3+x2+(T)2*x+(T)4/3;
                dweight(j)(3)(k)=((T).5*x2+(T)2*x+(T)2)*one_over_dx;
            } break;
            default: PHYSBAM_FATAL_ERROR();}
    }
};
}

#endif
