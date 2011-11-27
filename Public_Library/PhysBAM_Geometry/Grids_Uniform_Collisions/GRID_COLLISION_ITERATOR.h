//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_COLLISION_ITERATOR
//#####################################################################
#ifndef __GRID_COLLISION_ITERATOR__
#define __GRID_COLLISION_ITERATOR__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class TV> class COLLISION_GEOMETRY;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

template<class TV>
class GRID_COLLISION_ITERATOR
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    const GRID<TV>& grid;
    T thickness;
    int first_z_face;

    struct ENTRY
    {
        FACE_INDEX<TV::m> face;
        VECTOR<bool,2> inside;
        VECTOR<int,2> closest;

        bool operator<(const ENTRY&e) const
        {
            if(face.axis!=e.face.axis) return face.axis<e.face.axis;
            for(int a=1;a<=TV::m;a++) if(a!=face.axis) if(face.index(a)!=e.face.index(a)) return face.index(a)<e.face.index(a);
            return face.index(face.axis)<e.face.index(face.axis);
        }
    };
    ARRAY<ENTRY> faces;

    GRID_COLLISION_ITERATOR(const GRID<TV>& grid_input);
    GRID_COLLISION_ITERATOR(const GRID<TV>& grid_input,const TRIANGULATED_SURFACE<T>& surface);
    ~GRID_COLLISION_ITERATOR();

    void Initialize(const TRIANGULATED_SURFACE<T>& ts);
    void Initialize(const ARRAY<TRIANGLE_3D<T> >& elements);
    void Initialize(const COLLISION_GEOMETRY<TV>& cg);

    struct FACE_ITERATOR
    {
        const ARRAY<ENTRY>& faces;
        const GRID<TV>& grid;
        int cur;
        FACE_INDEX<TV::m> index;
        VECTOR<RANGE<TV_INT>,TV::dimension> face_domain;

        FACE_ITERATOR(GRID_COLLISION_ITERATOR<TV>& gci,int ghost_input=0);
        ~FACE_ITERATOR();

        const FACE_INDEX<TV::m>& Index() const PHYSBAM_ALWAYS_INLINE
        {return index;}

        const ENTRY& Entry() const PHYSBAM_ALWAYS_INLINE
        {return faces(cur);}

        bool Valid() const PHYSBAM_ALWAYS_INLINE
        {return cur<=faces.m;}

        void Next()
        {cur++;if(cur<=faces.m){index=faces(cur).face;if(!face_domain(index.axis).Lazy_Inside(index.index)) Next();}}
    };

    struct INTERIOR_FACE_ITERATOR
    {
        const ARRAY<ENTRY>& faces;
        const GRID<TV>& grid;
        int cur;
        int last;
        int ghost;
        FACE_INDEX<TV::m> index;
        VECTOR<RANGE<TV_INT>,TV::dimension> face_domain;

        INTERIOR_FACE_ITERATOR(GRID_COLLISION_ITERATOR<TV>& gci,int ghost_input=0);
        ~INTERIOR_FACE_ITERATOR();

        const FACE_INDEX<TV::m>& Index() const PHYSBAM_ALWAYS_INLINE
        {return index;}

        const ENTRY& Entry() const PHYSBAM_ALWAYS_INLINE
        {return faces(cur);}

        bool Valid() const PHYSBAM_ALWAYS_INLINE
        {return index.index(index.axis)<=last;}

        void Next() PHYSBAM_ALWAYS_INLINE
        {if(++index.index(index.axis)>last) Next_Helper();}

        void Next_Helper();
    };

    struct CELL_ITERATOR
    {
        const ARRAY<ENTRY>& faces;
        const GRID<TV>& grid;
        int cur;
        int last;
        TV_INT index;
        int ghost;

        CELL_ITERATOR(GRID_COLLISION_ITERATOR<TV>& gci,int ghost_input=0);
        ~CELL_ITERATOR();

        const TV_INT& Index() const PHYSBAM_ALWAYS_INLINE
        {return index;}

        bool Valid() const PHYSBAM_ALWAYS_INLINE
        {return index(TV::m)<=last;}

        void Next() PHYSBAM_ALWAYS_INLINE
        {if(++index(TV::m)>last) Next_Helper();}

        void Next_Helper();
    };
};
}
#endif
