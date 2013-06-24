//#####################################################################
// Copyright 2002-2011, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GNUPLOT_OUTPUT 
//#####################################################################
#ifndef __GNUPLOT_OUTPUT__
#define __GNUPLOT_OUTPUT__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <fstream>
namespace PhysBAM{
template<class TV> class GRID;
template<class T,int d> class VECTOR;
class GNUPLOT_OUTPUT
{
public: 
    GNUPLOT_OUTPUT();
    ~GNUPLOT_OUTPUT();

    template<class T> void Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& output,const int stepnumber);
    template<class T> void Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,2> >& grid,const ARRAY<T,VECTOR<int,2> >& output,const int stepnumber);
    template<class T> void Write_Output_File(const std::string& file_name,const GRID<VECTOR<T,3> >& grid,const ARRAY<T,VECTOR<int,3> >& output,const int stepnumber);
    template<class T,int d> void Write_Output_File(const std::string& file_name,const ARRAY_VIEW<VECTOR<T,d> >& X,const int stepnumber);
};
}
#endif
