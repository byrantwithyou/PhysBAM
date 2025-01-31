//#####################################################################
// Copyright 2008, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ARRAYS_BINARY_UNIFORM
//#####################################################################
#ifndef __FACE_ARRAYS_BINARY_UNIFORM__
#define __FACE_ARRAYS_BINARY_UNIFORM__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/SIDED_FACE_INDEX.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T,int d>
class ARRAY<T,SIDED_FACE_INDEX<d> >:public ARRAY<T,FACE_INDEX<d> >
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    struct UNUSABLE{};
public:
    enum {dimension=d};
    typedef ARRAY<T,FACE_INDEX<d> > BASE;
    typedef ARRAY_VIEW<T,TV_INT> T_ARRAY_VIEW;
    typedef int HAS_UNTYPED_READ_WRITE;

    template<class T2> struct REBIND{typedef ARRAY<T2,SIDED_FACE_INDEX<dimension> > TYPE;};
    typedef T ELEMENT;

    using BASE::operator();using BASE::Component;
    using BASE::data;using BASE::domain_indices;

    BASE u2;
    VECTOR<VECTOR<T_ARRAY_VIEW*,dimension>,2> sided_data;

    ARRAY()
    {Initialize();}

    ARRAY(const ARRAY& old_array,const bool initialize_with_old_array)
//        :u(old_array.u,initialize_with_old_array),v(old_array.v,initialize_with_old_array),w(old_array.w,initialize_with_old_array)
        :BASE(old_array.u,initialize_with_old_array),u2(old_array.u2) // TODO: this always initializes with the old_array
    {
        Initialize();
    }

    virtual ~ARRAY()
    {}

    void Initialize()
    {for(int i=0;i<dimension;i++){sided_data[0](i)=&BASE::Component(i);sided_data[1](i)=&u2.Component(i);}}

    template<class T2_GRID>
    ARRAY(const T2_GRID& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
    {Initialize();Resize(grid,ghost_cells,initialize_using_default_constructor);} // TODO: make T2_GRID only allow valid things

    template<class T2>
    ARRAY& operator*=(const T2 a)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i)*=a;
    return *this;}

    ARRAY& operator+=(const ARRAY& a)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i)+=a.Component(side,i);
    return *this;}

    ARRAY& operator-=(const ARRAY& a)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i)-=a.Component(side,i);
    return *this;}

    void Resize(const RANGE<TV_INT>& domain,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {BASE::Resize(domain,initialize_new_elements,copy_existing_elements,initialization_value);u2.Resize(domain,initialize_new_elements,copy_existing_elements,initialization_value);}

    template<class T2_GRID>
    void Resize(const T2_GRID& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(grid.Domain_Indices(ghost_cells),initialize_new_elements,copy_existing_elements,initialization_value);}

    void Clean_Memory()
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i).Clean_Memory();}

    void Delete_Pointers_And_Clean_Memory()
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i).Delete_Pointers_And_Clean_Memory();}

    int Number_Of_Ghost_Cells() const
    {return -domain_indices.min_corner.x;}

    T& operator()(const int side,const int axis,const TV_INT& face)
    {return Component(side,axis)(face);}

    const T& operator()(const int side,const int axis,const TV_INT& face) const
    {return Component(side,axis)(face);}

    T& operator()(const SIDED_FACE_INDEX<dimension>& index)
    {return Component(index.side,index.axis)(index.index);}

    const T& operator()(const SIDED_FACE_INDEX<dimension>& index) const
    {return Component(index.side,index.axis)(index.index);}

    T_ARRAY_VIEW& Component(const int side,const int axis)
    {assert((unsigned)axis<dimension);return *sided_data[side](axis);}

    const T_ARRAY_VIEW& Component(const int side,const int axis) const
    {assert((unsigned)axis<dimension);return *sided_data[side](axis);}

    TV Cell_Centered_Average(const TV_INT& cell_index) const
    {TV average;
    for(int axis=0;axis<dimension;axis++) average(axis)=(T).5*(Component(2,axis)(cell_index)+Component(1,axis)(cell_index+TV_INT::Axis_Vector(axis)));
    return average;}

    void Set_All_Faces(const T& value,const TV_INT& cell_index)
    {for(int axis=0;axis<dimension;axis++) Component(2,axis)(cell_index)=Component(1,axis)(cell_index+TV_INT::Axis_Vector(axis))=value;}

    static void Extract_Dimension(const ARRAY& old_array,ARRAY<typename ELEMENT_OF_VECTOR<T>::TYPE,SIDED_FACE_INDEX<dimension> >& extracted_array,const TV_INT& dimensions_to_extract)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) T_ARRAY_VIEW::Extract_Dimension(old_array.Component(side,i),extracted_array.Component(side,i),dimensions_to_extract(i));}

    void Fill(const T& constant)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i).Fill(constant);}

    void Copy(const ARRAY& old_copy,ARRAY& new_copy)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i).Copy(old_copy.Component(side,i));}

    template<class T2>
    void Copy(const T2 c,const ARRAY& old)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i).Copy(c,old.Component(side,i));}

    void Fill(const TV& value)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i).Fill(value(i));}

    template<class T2>
    void Copy(const T2 c1,const ARRAY& v1,const T2 c2,const ARRAY& v2)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) Component(side,i).Copy(c1,v1.Component(side,i),c2,v2.Component(side,i));}

    static void Put(const ARRAY& old_copy,ARRAY& new_copy)
    {for(int side=0;side<2;side++)for(int i=0;i<dimension;i++) T_ARRAY_VIEW::Put(old_copy.Component(side,i),new_copy.Component(side,i));}

    TV Max_Abs() const
    {TV maxabs_values;
    for(int i=0;i<dimension;i++) maxabs_values(i)=max(Component(1,i).Max_Abs(),Component(2,i).Max_Abs());
    return maxabs_values;}

    static void Exchange(ARRAY& a,ARRAY& b)
    {BASE::Exchange(a,b);BASE::Exchange(a.u2,b.u2);a.Initialize();b.Initialize();}
    
    template<class RW> void Read(std::istream& input)
    {BASE::Read(input);Read_Binary<RW>(input,u2);}

    template<class RW> void Write(std::ostream& output) const
    {BASE::Write(output);Write_Binary<RW>(output,u2);}

//#####################################################################
};
}
#endif
