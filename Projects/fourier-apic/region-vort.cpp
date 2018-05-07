//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Tools/Parsing/PARSE_ARGS.h>

using namespace PhysBAM;

typedef double RW;
typedef double T;
typedef VECTOR<T,2> TV;
typedef VECTOR<int,TV::m> TV_INT;

template<class TV>
void Vorticity_Analysis(const GRID<TV>& grid,const RANGE<TV>& domain,T time,
    const ARRAY<T,FACE_INDEX<TV::m> >& mass,const ARRAY<T,FACE_INDEX<TV::m> >& velocity,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N)
{
    RANGE<TV_INT> cell_domain=grid.Clamp_To_Cell(domain);
    auto valid=[&](FACE_INDEX<TV::m> face){return mass(face) && (psi_N.domain_indices.Empty() || !psi_N(face));};
    decltype(MATRIX<T,TV::m>().Contract_Permutation_Tensor()) vort;
    for(CELL_ITERATOR<TV> it(grid,cell_domain);it.Valid();it.Next()){
        bool ok=true;
        MATRIX<T,TV::m> dV;
        for(int i=0;i<TV::m;i++)
            for(int j=0;j<TV::m;j++){
                if(i==j) continue;
                FACE_INDEX<TV::m> A(i,it.index),B(A),C(A),D(A);
                C.index(i)++;
                D.index(i)++;
                A.index(j)++;
                B.index(j)--;
                C.index(j)++;
                D.index(j)--;
                if(valid(A) && valid(B) && valid(C) && valid(D))
                    dV(i,j)=(T).25*grid.one_over_dX(j)*(
                        velocity(A)-velocity(B)+
                        velocity(C)-velocity(D));
                else ok=false;}
        if(!ok) continue;
        vort+=dV.Contract_Permutation_Tensor();}
    LOG::printf("window vorticity: %P %P\n",time,vort);
}

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    std::string directory;
    parse_args.Extra(&directory,"dir","simulation result directory");
    T x0,x1,y0,y1;
    parse_args.Extra(&x0,"float","x0");
    parse_args.Extra(&y0,"float","y0");
    parse_args.Extra(&x1,"float","x1");
    parse_args.Extra(&y1,"float","y1");
    parse_args.Parse();

    int last_step;
    RANGE<TV> domain(TV(x0,y0),TV(x1,y1));
    Read_From_Text_File(directory+"/common/last_grid_data",last_step);
    GRID<TV> grid;
    Read_From_File(directory+"/common/grid",grid);
    ARRAY<T,FACE_INDEX<TV::m> > mass, velocity;
    ARRAY<bool,FACE_INDEX<TV::m> > psi_N;
    T time;
    for(int i=0;i<=last_step;i++){
        std::string d=LOG::sprintf("%s/v%d",directory,i);
        Read_From_File(d+"/mass",mass);
        Read_From_File(d+"/mac_velocities",velocity);
        Read_From_File(d+"/psi_N",psi_N);
        Read_From_Text_File(d+"/time",time);
        Vorticity_Analysis(grid,domain,time,mass,velocity,psi_N);}
    return 0;
}

