//#####################################################################
// Copyright 2005, Jerry Talton
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

template<class T> static int compare_array_min_m(const void *e1, const void *e2){if((*(T*)(e1)).m<(*(T*)(e2)).m)return -1;else return 0;}
template<class T> static int compare_array_max_m(const void *e1, const void *e2){if((*(T*)(e1)).m>(*(T*)(e2)).m)return -1;else return 0;}
template<class T> static int compare_vector_min_x(const void *e1, const void *e2){if((*(T*)(e1)).x<(*(T*)(e2)).x)return -1;else return 0;}
template<class T> static int compare_vector_max_x(const void *e1, const void *e2){if((*(T*)(e1)).x>(*(T*)(e2)).x)return -1;else return 0;}
