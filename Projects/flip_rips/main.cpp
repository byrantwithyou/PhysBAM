//#####################################################################
// Copyright 2014, Alexey Stomakhin
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include "FLIP_DRIVER.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    const int dimension=2;
    const int order=1;

    typedef float T;
    typedef VECTOR<T,dimension> TV;
    typedef VECTOR<int,dimension> TV_INT;
    typedef FLIP_DRIVER<TV,order> T_DRIVER;
    typedef T_DRIVER::T_PARTICLE T_PARTICLE;

    LOG::Initialize_Logging();

    int threads=15;
    int resolution=32;
    int last_frame=100;
    T max_dt=1e-3;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-threads",&threads,"threads","Number of threads to use");
    parse_args.Add("-resolution",&resolution,"resolution","Grid resolution");
    parse_args.Add("-last_frame",&last_frame,"last frame","Last simulation frame");
    parse_args.Add("-max_dt",&max_dt,"max dt","Maximum timestep allowed");
    parse_args.Parse();

    omp_set_num_threads(threads);
#pragma omp parallel
#pragma omp master
    if(omp_get_num_threads()!=threads){
        LOG::cerr<<"OpenMP error"<<std::endl;
        exit(-1);}
    else LOG::cout<<"Running on "<<threads<<" threads"<<std::endl;

    RANGE<TV> domain(TV(),TV()+1);
    TV_INT counts=TV_INT()+resolution;
    GRID<TV> grid(counts,domain,true);
    T_DRIVER driver(grid,threads);
    driver.max_dt=max_dt;
    driver.frames=last_frame;
    
    RANGE<TV> seed_domain1(TV(.25,0.25),TV(.75,.65));
    TV_INT seed_counts1=TV_INT()+resolution;
    GRID<TV> seed_grid1(seed_counts1,seed_domain1,true);

    for(CELL_ITERATOR<TV> it(seed_grid1);it.Valid();it.Next()){
        T_PARTICLE p;
        p.X=it.Location();
        driver.particles.Append(p);}

    RANGE<TV> seed_domain2(TV(.35,0.7),TV(.65,.9));
    TV_INT seed_counts2=TV_INT()+resolution/2;
    GRID<TV> seed_grid2(seed_counts2,seed_domain2,true);

    for(CELL_ITERATOR<TV> it(seed_grid2);it.Valid();it.Next()){
        T_PARTICLE p;
        p.X=it.Location();
        driver.particles.Append(p);}

    FLIP_COLLIDABLE_OBJECT_STATIC<TV>* bottom_wall=
        new FLIP_COLLIDABLE_OBJECT_STATIC<TV>(
            new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(
                RANGE<TV>(TV(-1,-1),TV(2,0.01))),TV(),ROTATION<TV>());
    driver.collidable_objects.Append(bottom_wall);

    FLIP_COLLIDABLE_OBJECT_STATIC<TV>* left_wall=
        new FLIP_COLLIDABLE_OBJECT_STATIC<TV>(
            new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(
                RANGE<TV>(TV(-1,-1),TV(0.01,2))),TV(),ROTATION<TV>());
    driver.collidable_objects.Append(left_wall);

    FLIP_COLLIDABLE_OBJECT_STATIC<TV>* right_wall=
        new FLIP_COLLIDABLE_OBJECT_STATIC<TV>(
            new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(
                RANGE<TV>(TV(0.99,-1),TV(2,2))),TV(),ROTATION<TV>());
    driver.collidable_objects.Append(right_wall);

    FLIP_COLLIDABLE_OBJECT_STATIC<TV>* top_wall=
        new FLIP_COLLIDABLE_OBJECT_STATIC<TV>(
            new ANALYTIC_IMPLICIT_OBJECT<RANGE<TV> >(
                RANGE<TV>(TV(-1,0.99),TV(2,2))),TV(),ROTATION<TV>());
    driver.collidable_objects.Append(top_wall);
    
    driver.Run();

    return 0;
}
