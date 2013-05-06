//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __MPLE_POINT__
#define __MPLE_POINT__

namespace PhysBAM{

template<class TV,int w>
class MPLE_POINT
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    enum WORKAROUND{ghost=3};

public:

    TV X;                     // position

    TV_INT base_node;         // corresponds to the min most grid point the particle has influence on
    VECTOR<TV,w> weights;     // 'w' weights per dimension 

    MPLE_POINT(){}

    ~MPLE_POINT(){}

    void Update_Base_And_Weights(const GRID<TV>& grid)
    {
        base_node=grid.Cell(X-(T).5*(w-2)*grid.DX(),ghost*2);
        const TV X_eval=X-grid.Node(base_node);
        for(int k=0;k<TV::m;k++){
            const T x=X_eval(k);
            const T one_over_h=(T)1/grid.DX()(k);
            Compute_Weights(x*one_over_h,one_over_h,k);}
    }

private:

    void Compute_Weights(T x,const T one_over_h,const int k)
    {
        switch(w){
            case 2:{
                weights(0)(k)=(T)1-x;
                x-=(T)1;
                weights(1)(k)=(T)1+x;
            } break;
            case 3:{
                weights(0)(k)=(T).5*x*x-(T)1.5*x+(T)1.125;
                x-=(T)1;
                weights(1)(k)=-x*x+(T).75;
                x-=(T)1;
                weights(2)(k)=(T).5*x*x+(T)1.5*x+(T)1.125;
            } break;
            case 4:{
                T x2,x3;
                x2=x*x;x3=x2*x;
                weights(0)(k)=-((T)1/6)*x3+x2-(T)2*x+(T)4/3;
                x-=(T)1;
                x2=x*x;x3=x2*x;
                weights(1)(k)=(T).5*x3-x2+(T)2/3;
                x-=(T)1;
                x2=x*x;x3=x2*x;
                weights(2)(k)=-(T).5*x3-x2+(T)2/3;
                x-=(T)1;
                x2=x*x;x3=x2*x;
                weights(3)(k)=((T)1/6)*x3+x2+(T)2*x+(T)4/3;
            } break;
            default: PHYSBAM_FATAL_ERROR();}
    }
};
}
#endif
