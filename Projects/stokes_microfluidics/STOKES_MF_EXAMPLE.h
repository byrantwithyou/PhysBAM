//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __STOKES_MF_EXAMPLE__
#define __STOKES_MF_EXAMPLE__
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Grid_Tools/Grids/GRID.h>
#include <functional>
namespace PhysBAM{

template<class TV> class DEBUG_PARTICLES;
class PARSE_ARGS;

template<class TV>
class STOKES_MF_EXAMPLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    STREAM_TYPE stream_type;

    GRID<TV> grid;
    ARRAY<TV_INT> vertices;
    HASHTABLE<PAIR<int,int> > edges;
    
    T density=1;
    T viscosity=0.1;

    std::string output_directory;
    std::string data_directory;
    int test_number=0;
    int resolution=8;
    ARRAY<TV_INT> iverts;
    ARRAY<PAIR<int,int> > iedges;

    // debugging
    std::unique_ptr<DEBUG_PARTICLES<TV> > debug_particles;

    STOKES_MF_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args);
    STOKES_MF_EXAMPLE(const STOKES_MF_EXAMPLE&) = delete;
    void operator=(const STOKES_MF_EXAMPLE&) = delete;
    virtual ~STOKES_MF_EXAMPLE();
    
    void Write_Output_Files(const int frame);
    void Read_Output_Files(const int frame);
    virtual void Initialize()=0;

    int Add_Vertex(const TV_INT& x);
    int Add_Edge(int v,int dir,int len);
    void Add_Edge(int v0,int v1);
    void Build_Grid(T len_scale,int cross_section_radius=1);
//#####################################################################
};
}
#endif
