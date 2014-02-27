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
    
    template<class TV>
    struct POTENTIAL_CENTER
    {
        TV sum;
        int n;
        bool set;
        
        POTENTIAL_CENTER():n(0),set(0)
        {
        }
        
        TV Value()
        {
            if(!n){
                TV c;
                for(int i=0;i<TV::m;++i)
                    c(i)=1./TV::m;
                return c;
            }
            return sum/n;
        }
        
        void Add(const TV& c){
            if(!set){
                sum+=c;
                ++n;
            }
        }
    };
    
    struct TRI_CUTTING
    {
        ARRAY<int> material;
        VECTOR<bool,12> turned_on;
        VECTOR<POTENTIAL_CENTER<TV>,3> edge_centers;
        POTENTIAL_CENTER<T3> face_center;
        
        TRI_CUTTING()
        {
            for(int i=0;i<6;++i)
                material.Append(i);
        }
    };
    ARRAY<TRI_CUTTING> tri_cuttings;
    
    CUTTING(TRIANGULATED_AREA<T>& ta,SEGMENTED_CURVE<TV>& sc);
    ~CUTTING(){}

    void Run(T tol);
};
}
#endif
