//#####################################################################
// Copyright 2014, Yuting Wang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_2D
//##################################################################### 
#ifndef __CUTTING_2D__
#define __CUTTING_2D__

#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/PAIR.h>
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
    typedef VECTOR<int,5> I5;
    typedef VECTOR<T,3> T3;
    typedef PAIR<int,T3> PS;
    enum WORKAROUND {d=TV::m};
    
public:
    TRIANGULATED_AREA<T>* sim_ta;
    TRIANGULATED_AREA<T>* ta;
    SEGMENTED_CURVE<TV>* sc;
    
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
            if(!n||!set) return TV()+(T)1/TV::m;
            return sum/n;
        }
        
        void Add(const TV& c){
            if(!set){
                sum+=c;
                n++;
            }
        }
    };
    
    struct TRI_CUTTING
    {
        VECTOR<bool,6> materials;
        VECTOR<bool,12> turned_on;
        VECTOR<POTENTIAL_CENTER<TV>,3> edge_centers;
        POTENTIAL_CENTER<T3> face_center;
        
        TRI_CUTTING()
        {
            materials.Fill(1);
        }
    };
    
    ARRAY<TRI_CUTTING> tri_cuttings;
    ARRAY<int> tri_in_sim;
    ARRAY<PS> particle_in_sim;
    
    CUTTING(TRIANGULATED_AREA<T>* sim_ta_,SEGMENTED_CURVE<TV>* sc_);
    ~CUTTING(){}

    T3 Weight_In_Tri(const I3& tri,const I3& tri1,const T3 w1)
    {
        T3 w;
        for(int l=0;l<3;++l)
            for(int m=0;m<3;++m)
                if(tri(l)==tri1(m))
                    w(l)=w1(m);
        return w;
    }
    
    T3 Weight_In_Sim(int tri_id,const T3& w)
    {
        int pid=tri_in_sim(tri_id);
        I3 ptri=sim_ta->mesh.elements(pid);
        I3 tri=ta->mesh.elements(tri_id);
        T3 weight;
        for(int i=0;i<3;++i){
            int pid1=particle_in_sim(tri(i)).x;
            I3 ptri1=sim_ta->mesh.elements(pid1);
            T3 w1=particle_in_sim(tri(i)).y;
            T3 ww=Weight_In_Tri(ptri,ptri1,w1);
            weight+=(ww*w[i]);
        }
        return weight;
    }
    
    void Run(T tol);
    void Update_Material_Particles();
};
}
#endif
