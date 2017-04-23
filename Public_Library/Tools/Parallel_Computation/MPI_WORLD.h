//#####################################################################
// Copyright 2005-2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_WORLD
//#####################################################################
#ifndef __MPI_WORLD__
#define __MPI_WORLD__

namespace PhysBAM{

class PARSE_ARGS;
class MPI_WORLD
{
public:
    bool initialized;
    int rank;

    MPI_WORLD();
    MPI_WORLD(PARSE_ARGS& parse_args);
    MPI_WORLD(const MPI_WORLD&) = delete;
    void operator=(const MPI_WORLD&) = delete;
    ~MPI_WORLD();

//#####################################################################
    static bool Initialized();
private:
    void Initialize(bool force_mpi);
//#####################################################################
};
}
#endif
