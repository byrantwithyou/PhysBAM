//#####################################################################
// Copyright 2014, Yuting Wang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_2D
//##################################################################### 
#ifndef __CUTTING_2D__
#define __CUTTING_2D__

#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class CONSISTENT_INTERSECTIONS;
template<class T> class TRIANGULATED_AREA;
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class TV> class SEGMENTED_CURVE;
template<class T> class TRIANGULATED_SURFACE;
template<class TV> class CUTTING;

template<class T>
class CUTTING<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> I2;
    typedef VECTOR<int,3> I3;
    typedef VECTOR<int,4> I4;
    typedef VECTOR<T,3> T3;
    enum WORKAROUND {d=TV::m};
    
public:
    TRIANGULATED_AREA<T>& ta;
    SEGMENTED_CURVE<TV>& sc;
    
    struct TRI_CUTTING
    {
        ARRAY<int> material;
        VECTOR<bool,12> turned_on;
        VECTOR<TV,3> edge_centers;
        T3 face_center;
        
        TRI_CUTTING()
        {
            for(int i=0;i<6;++i)
                material.Append(i);
            for(int i=0;i<3;++i)
                edge_centers(i)=TV(0.5,0.5);
            face_center=T3(1./3,1./3,1./3);
        }
    };
    ARRAY<TRI_CUTTING> tri_cuttings;
    
    CUTTING(TRIANGULATED_AREA<T>& ta,SEGMENTED_CURVE<TV>& sc);
    ~CUTTING(){}

    void Run(T tol);
};
}
#endif
