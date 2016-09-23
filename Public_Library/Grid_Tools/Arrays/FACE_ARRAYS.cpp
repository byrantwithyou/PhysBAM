//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
namespace PhysBAM{
template ARRAY<bool,FACE_INDEX<1> >::~ARRAY();
template ARRAY<bool,FACE_INDEX<2> >::~ARRAY();
template ARRAY<bool,FACE_INDEX<3> >::~ARRAY();
template ARRAY<int,FACE_INDEX<1> >::~ARRAY();
template ARRAY<int,FACE_INDEX<2> >::~ARRAY();
template ARRAY<int,FACE_INDEX<3> >::~ARRAY();
template ARRAY<float,FACE_INDEX<1> >::~ARRAY();
template ARRAY<float,FACE_INDEX<2> >::~ARRAY();
template ARRAY<float,FACE_INDEX<3> >::~ARRAY();
template ARRAY<double,FACE_INDEX<1> >::~ARRAY();
template ARRAY<double,FACE_INDEX<2> >::~ARRAY();
template ARRAY<double,FACE_INDEX<3> >::~ARRAY();

template ARRAY<VECTOR<double,1>, FACE_INDEX<1> >::~ARRAY();
template ARRAY<VECTOR<double,2>, FACE_INDEX<2> >::~ARRAY();
template ARRAY<VECTOR<double,3>, FACE_INDEX<3> >::~ARRAY();
template ARRAY<VECTOR<float,1>, FACE_INDEX<1> >::~ARRAY();
template ARRAY<VECTOR<float,2>, FACE_INDEX<2> >::~ARRAY();
template ARRAY<VECTOR<float,3>, FACE_INDEX<3> >::~ARRAY();

template ARRAY<VECTOR<double,3>, FACE_INDEX<1> >::~ARRAY();
template ARRAY<VECTOR<double,4>, FACE_INDEX<2> >::~ARRAY();
template ARRAY<VECTOR<double,5>, FACE_INDEX<3> >::~ARRAY();
template ARRAY<VECTOR<float,3>, FACE_INDEX<1> >::~ARRAY();
template ARRAY<VECTOR<float,4>, FACE_INDEX<2> >::~ARRAY();
template ARRAY<VECTOR<float,5>, FACE_INDEX<3> >::~ARRAY();

template void ARRAY<bool, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const bool&);
template void ARRAY<bool, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const bool&);
template void ARRAY<bool, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const bool&);
template void ARRAY<int, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const int&);
template void ARRAY<int, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const int&);
template void ARRAY<int, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const int&);
template void ARRAY<float, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const float&);
template void ARRAY<float, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const float&);
template void ARRAY<float, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const float&);
template void ARRAY<double, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const double&);
template void ARRAY<double, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const double&);
template void ARRAY<double, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const double&);

template void ARRAY<VECTOR<double,1>, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const VECTOR<double,1>&);
template void ARRAY<VECTOR<double,2>, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const VECTOR<double,2>&);
template void ARRAY<VECTOR<double,3>, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const VECTOR<double,3>&);
template void ARRAY<VECTOR<float,1>, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const VECTOR<float,1>&);
template void ARRAY<VECTOR<float,2>, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const VECTOR<float,2>&);
template void ARRAY<VECTOR<float,3>, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const VECTOR<float,3>&);

template void ARRAY<VECTOR<double,3>, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const VECTOR<double,3>&);
template void ARRAY<VECTOR<double,4>, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const VECTOR<double,4>&);
template void ARRAY<VECTOR<double,5>, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const VECTOR<double,5>&);
template void ARRAY<VECTOR<float,3>, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const VECTOR<float,3>&);
template void ARRAY<VECTOR<float,4>, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const VECTOR<float,4>&);
template void ARRAY<VECTOR<float,5>, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const VECTOR<float,5>&);

template void ARRAY<int, FACE_INDEX<1> >::Read<float>(std::istream&);
template void ARRAY<int, FACE_INDEX<2> >::Read<float>(std::istream&);
template void ARRAY<int, FACE_INDEX<3> >::Read<float>(std::istream&);
template void ARRAY<float, FACE_INDEX<1> >::Read<float>(std::istream&);
template void ARRAY<float, FACE_INDEX<2> >::Read<float>(std::istream&);
template void ARRAY<float, FACE_INDEX<3> >::Read<float>(std::istream&);
template void ARRAY<double, FACE_INDEX<1> >::Read<float>(std::istream&);
template void ARRAY<double, FACE_INDEX<2> >::Read<float>(std::istream&);
template void ARRAY<double, FACE_INDEX<3> >::Read<float>(std::istream&);
template void ARRAY<int, FACE_INDEX<1> >::Read<double>(std::istream&);
template void ARRAY<int, FACE_INDEX<2> >::Read<double>(std::istream&);
template void ARRAY<int, FACE_INDEX<3> >::Read<double>(std::istream&);
template void ARRAY<float, FACE_INDEX<1> >::Read<double>(std::istream&);
template void ARRAY<float, FACE_INDEX<2> >::Read<double>(std::istream&);
template void ARRAY<float, FACE_INDEX<3> >::Read<double>(std::istream&);
template void ARRAY<double, FACE_INDEX<1> >::Read<double>(std::istream&);
template void ARRAY<double, FACE_INDEX<2> >::Read<double>(std::istream&);
template void ARRAY<double, FACE_INDEX<3> >::Read<double>(std::istream&);

template void ARRAY<int, FACE_INDEX<1> >::Write<float>(std::ostream&) const;
template void ARRAY<int, FACE_INDEX<2> >::Write<float>(std::ostream&) const;
template void ARRAY<int, FACE_INDEX<3> >::Write<float>(std::ostream&) const;
template void ARRAY<float, FACE_INDEX<1> >::Write<float>(std::ostream&) const;
template void ARRAY<float, FACE_INDEX<2> >::Write<float>(std::ostream&) const;
template void ARRAY<float, FACE_INDEX<3> >::Write<float>(std::ostream&) const;
template void ARRAY<double, FACE_INDEX<1> >::Write<float>(std::ostream&) const;
template void ARRAY<double, FACE_INDEX<2> >::Write<float>(std::ostream&) const;
template void ARRAY<double, FACE_INDEX<3> >::Write<float>(std::ostream&) const;
template void ARRAY<int, FACE_INDEX<1> >::Write<double>(std::ostream&) const;
template void ARRAY<int, FACE_INDEX<2> >::Write<double>(std::ostream&) const;
template void ARRAY<int, FACE_INDEX<3> >::Write<double>(std::ostream&) const;
template void ARRAY<float, FACE_INDEX<1> >::Write<double>(std::ostream&) const;
template void ARRAY<float, FACE_INDEX<2> >::Write<double>(std::ostream&) const;
template void ARRAY<float, FACE_INDEX<3> >::Write<double>(std::ostream&) const;
template void ARRAY<double, FACE_INDEX<1> >::Write<double>(std::ostream&) const;
template void ARRAY<double, FACE_INDEX<2> >::Write<double>(std::ostream&) const;
template void ARRAY<double, FACE_INDEX<3> >::Write<double>(std::ostream&) const;
}
