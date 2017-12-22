//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include<Core/Log/DEBUG_UTILITIES.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Optimization/IPOPT.h>
#include <Tools/Parsing/PARSE_ARGS.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

class TEST_FUNC:public IPOPT<T>
{
public:
    T restitution,friction,m;
    TV v,n;

    TEST_FUNC(const TV& j0,const T& r,const T& f,const T& mass,const TV& vel,const TV& normal):
        restitution(r),friction(f),m(mass),v(vel),n(normal)
    {
        this->initial_guess=j0;
        tolerance=1e-7;
        dof_range.Resize(2,true,true,INTERVAL<T>::Full_Box());
        constraint_range=VECTOR<INTERVAL<T>,2>({0,FLT_MAX},{0,FLT_MAX});
        Compute_Hessian(initial_guess);
        Compute_Constraint_Jacobian(initial_guess);
    }

    virtual T Compute_Objective(const ARRAY<T>& x) override
    {
        TV j(x);
//        LOG::printf("E %P %P\n",j,(T)0.5*m*(v+j/m).Magnitude_Squared());
        return (T)0.5*m*(v+j/m).Magnitude_Squared();
    }

    virtual void Compute_Gradient(ARRAY<T>& g,const ARRAY<T>& x) override
    {
        TV j(x);
        g=v+j/m;
//        LOG::printf("dE %P %P\n",j,g);
    }

    virtual void Compute_Hessian(const ARRAY<T>& x) override
    {
        TV j(x);
        for(int i=0;i<j.m;i++) H.Set({i,i},1/m);
    }

    virtual void Compute_Constraints(ARRAY<T>& f,const ARRAY<T>& x) override
    {
        TV j(x);
        f(0)=n.Dot(v+j/m)-max((T)0,-restitution*n.Dot(v));
        f(1)=(sqr(friction)+1)*sqr(j.Dot(n))-j.Magnitude_Squared();
//        LOG::printf("C %P %P\n",x,f);
//        f(2)=j.Dot(n);
    }

    virtual void Compute_Constraint_Jacobian(const ARRAY<T>& x) override
    {
        TV j(x),J0=n/m,J1=2*(sqr(friction)+1)*j.Dot(n)*n-2*j;
//        LOG::printf("dC %P %P %P\n",j,J0,J1);
        for(int i=0;i<j.m;i++){
            J.Set({0,i},J0(i));
            J.Set({1,i},J1(i));
//            J.Set({2,i},n(i));
        }
    }
};

int main(int argc, char* argv[])
{
    TV initial_guess,vel;
    T restitution=0,friction=0;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-s",&initial_guess,"float","initial guess");
    parse_args.Add("-v",&vel,"float","velocity");
    parse_args.Add("-restitution",&restitution,"float","restitution");
    parse_args.Add("-friction",&friction,"float","friction");
    parse_args.Parse();

    PHYSBAM_ASSERT(vel.Dot(TV(0,1))<=0);

    LOG::printf("initial_guess: %P\n",initial_guess);
    LOG::printf("velocity: %P\n",vel);
    LOG::printf("restitution: %f\n",restitution);
    LOG::printf("friction: %f\n",friction);

    TV sol;
    sol.y=-restitution*vel.y;
    if(fabs(vel.x)<friction*(sol.y-vel.y))
        sol.x=0;
    else sol.x=vel.x-sign(vel.x)*friction*(sol.y-vel.y);

    TEST_FUNC tf(initial_guess,restitution,friction,1,vel,TV(0,1));
    tf.Solve();
    LOG::printf("sol: %P\n",tf.solution);
    LOG::printf("new velocity: %P\n",vel+tf.solution/tf.m);
    LOG::printf("analytic velocity: %P\n",sol);
    LOG::printf("ERROR: %P\n",(vel+tf.solution/tf.m-sol).Magnitude());

    
    return 0;
}
