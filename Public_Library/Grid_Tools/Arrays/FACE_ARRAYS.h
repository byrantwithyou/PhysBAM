//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ARRAY
//#####################################################################
#ifndef __FACE_ARRAY__
#define __FACE_ARRAY__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Arrays_Nd/ARRAYS_ND_VIEW.h>
#include <Core/Vectors/SCALAR_POLICY.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/GRID.h>
namespace PhysBAM{

template<class T,int d>
class ARRAY<T,FACE_INDEX<d> >
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    struct UNUSABLE{};
public:
    ARRAY_VIEW<T> array; // one-dimensional data storage

    enum {dimension=TV::m};
    template<class T2> struct REBIND{typedef ARRAY<T2,FACE_INDEX<dimension> > TYPE;};
    typedef ARRAY_VIEW<T,TV_INT> T_ARRAY_VIEW;
    typedef T ELEMENT;typedef FACE_INDEX<d> INDEX;typedef T& RESULT_TYPE;typedef const T& CONST_RESULT_TYPE;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;
    typedef int HAS_UNTYPED_READ_WRITE;

    RANGE<TV_INT> domain_indices;
    VECTOR<T_ARRAY_VIEW,dimension> data;
   
    ARRAY()
    {}

    template<class T2>
    ARRAY(const GRID<VECTOR<T2,dimension> >& grid,const int ghost_cells=0,const bool initialize_using_default_constructor=true)
    {Resize(grid,ghost_cells,initialize_using_default_constructor);}

    ARRAY(const RANGE<TV_INT>& domain_indices_input,const bool initialize_using_default_constructor=true)
    {Resize(domain_indices_input,initialize_using_default_constructor);}

    ARRAY(const ARRAY& old_array)
    {Resize(old_array.Domain_Indices(),false,false);Copy(old_array);}

    ARRAY(ARRAY&& old_array) = default;

    virtual ~ARRAY()
    {delete[] array.base_pointer;}

    template<class T2>
    ARRAY& operator*=(const T2 a)
    {for(int i=0;i<dimension;i++) data(i)*=a;
    return *this;}

    ARRAY& operator+=(const ARRAY& a)
    {data+=a.data;return *this;}

    ARRAY& operator-=(const ARRAY& a)
    {data-=a.data;return *this;}

    const RANGE<TV_INT>& Domain_Indices() const
    {return domain_indices;}

    ARRAY& operator=(const ARRAY& source)
    {if(source.Domain_Indices()!=Domain_Indices()) Resize(source.Domain_Indices(),false,false);Copy(source);return *this;}

    ARRAY& operator=(ARRAY&&) = default;

    template<class T2>
    void Resize(const GRID<VECTOR<T2,dimension> >& grid,const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {Resize(grid.Domain_Indices(ghost_cells),initialize_new_elements,copy_existing_elements,initialization_value);}

    virtual void Resize(const RANGE<TV_INT>& domain,const bool initialize_new_elements=true,const bool copy_existing_elements=true,const T& initialization_value=T())
    {TV_INT sizes_new;bool early_return=true;
    VECTOR<RANGE<TV_INT>,dimension> domains;
    domain_indices=domain;
    for(int i=0;i<dimension;i++){domains(i)=domain;domains(i).max_corner(i)++;sizes_new(i)=(domains(i).Edge_Lengths()).Product();if(domains(i)!=data(i).Domain_Indices()) early_return=false;}
    if(early_return) return;
    int new_length=sizes_new.Sum();
    T* p=new T[new_length];
    if(initialize_new_elements){T* p_start=p;
        for(int i=0;i<dimension;i++){
            T_ARRAY_VIEW array_new(domains(i),p_start);
            array_new.Fill(initialization_value);
            p_start+=sizes_new(i);}}
    if(copy_existing_elements){T* p_start=p;
        for(int i=0;i<dimension;i++){
            T_ARRAY_VIEW array_new(domains(i),p_start);
            T_ARRAY_VIEW::Limited_Shifted_Get(array_new,data(i),TV_INT());
            p_start+=sizes_new(i);}}
    T* p_start=p;
    for(int i=0;i<dimension;i++){
        T_ARRAY_VIEW array_new(domains(i),p_start);
        T_ARRAY_VIEW::Exchange(array_new,data(i));
        p_start+=sizes_new(i);}
    delete[] array.base_pointer;
    array.Set(new_length,p);}

    void Clean_Memory()
    {Resize(RANGE<TV_INT>(TV_INT::All_Ones_Vector(),TV_INT()),false,false);}

    void Delete_Pointers_And_Clean_Memory()
    {for(int i=0;i<array.m;i++) delete array.base_pointer[i];Clean_Memory();}

    int Number_Of_Ghost_Cells() const
    {return -domain_indices.min_corner.x;}

    T& operator()(const int axis,const TV_INT& face)
    {return Component(axis)(face);}

    const T& operator()(const int axis,const TV_INT& face) const
    {return Component(axis)(face);}

    T& operator()(const FACE_INDEX<dimension>& index)
    {return Component(index.axis)(index.index);}

    const T& operator()(const FACE_INDEX<dimension>& index) const
    {return Component(index.axis)(index.index);}

    bool Valid_Index(const FACE_INDEX<dimension>& index) const
    {return ((unsigned)index.axis<(unsigned)dimension) && Component(index.axis).Valid_Index(index.index);}

    int Standard_Index(const FACE_INDEX<TV::m>& index) const
    {assert((unsigned)index.axis<(unsigned)dimension);
    return data(index.axis).array.base_pointer-array.base_pointer+data(index.axis).Standard_Index(index.index);}

    T_ARRAY_VIEW& Component(const int axis)
    {assert((unsigned)axis<dimension);return data(axis);}

    const T_ARRAY_VIEW& Component(const int axis) const
    {assert((unsigned)axis<dimension);return data(axis);}

    TV Cell_Centered_Average(const TV_INT& cell_index) const
    {TV average;
    for(int i=0;i<dimension;i++) average(i)=(T).5*(data(i)(cell_index)+data(i)(cell_index+TV_INT::Axis_Vector(i)));
    return average;}

    void Set_All_Faces(const T& value,const TV_INT& cell_index)
    {for(int i=0;i<dimension;i++) data(i)(cell_index)=data(i)(cell_index+TV_INT::Axis_Vector(i))=value;}

    static void Extract_Dimension(const ARRAY& old_array,ARRAY<typename ELEMENT_OF_VECTOR<T>::TYPE,FACE_INDEX<dimension> >& extracted_array,const TV_INT& dimensions_to_extract)
    {for(int i=0;i<dimension;i++) T_ARRAY_VIEW::Extract_Dimension(old_array.data(i),extracted_array.data(i),dimensions_to_extract(i));}

    void Fill(const T& constant)
    {for(int i=0;i<dimension;i++) data(i).Fill(constant);}

    void Copy(const ARRAY& old_copy)
    {for(int i=0;i<dimension;i++) data(i).Copy(old_copy.data(i));}

    template<class T2>
    void Copy(const T2 c,const ARRAY& old)
    {for(int i=0;i<dimension;i++) data(i).Copy(c,old.data(i));}

    void Fill(const TV& value)
    {for(int i=0;i<dimension;i++) data(i).Fill(value.x);}

    template<class T2>
    void Copy(const T2 c1,const ARRAY& v1,const T2 c2,const ARRAY& v2)
    {for(int i=0;i<dimension;i++) data(i).Copy(c1,v1.data(i),c2,v2.data(i));}

    static void Put(const ARRAY& old_copy,ARRAY& new_copy)
    {for(int i=0;i<dimension;i++) T_ARRAY_VIEW::Put(old_copy.data(i),new_copy.data(i));}

    TV Max_Abs() const
    {TV maxabs_vector;
    for(int i=0;i<dimension;i++) maxabs_vector(i)=data(i).Max_Abs();
    return maxabs_vector;}

    static void Exchange(ARRAY& a,ARRAY& b)
    {a.Exchange(b);}

    void Exchange(ARRAY& a)
    {exchange(domain_indices,a.domain_indices);array.Exchange(a.array);
    for(int i=0;i<dimension;i++) T_ARRAY_VIEW::Exchange(data(i),a.data(i));}

    template<class T_INDICES>
    INDIRECT_ARRAY<ARRAY,T_INDICES&> Subset(const T_INDICES& indices)
    {return INDIRECT_ARRAY<ARRAY,T_INDICES&>(*this,indices);}

    template<class T_INDICES>
    INDIRECT_ARRAY<const ARRAY,T_INDICES&> Subset(const T_INDICES& indices) const
    {return INDIRECT_ARRAY<const ARRAY,T_INDICES&>(*this,indices);}
    
    template<class RW> void Read(std::istream& input)
    {Clean_Memory();Read_Binary<RW>(input,domain_indices);Read_Binary<RW>(input,array.m);
    if(array.m<0){
        char buff[100];
        sprintf(buff,"Invalid negative array size %d",array.m);
        throw READ_ERROR(buff);}
    if(!array.m) return;
    array.base_pointer=new T[array.m];
    Read_Binary_Array<RW>(input,array.base_pointer,array.m);
    T* p_start=array.base_pointer;
    for(int i=0;i<d;i++){
        RANGE<TV_INT> domain;
        Read_Binary<RW>(input,domain);
        T_ARRAY_VIEW array_new(domain,p_start);
        T_ARRAY_VIEW::Exchange(array_new,data(i));
        p_start+=(domain.Edge_Lengths()).Product();}}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,domain_indices);
    Write_Binary<RW>(output,array.m);
    Write_Binary_Array<RW>(output,array.base_pointer,array.m);
    for(int i=0;i<d;i++) Write_Binary<RW>(output,data(i).domain);}
//#####################################################################
};

template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const ARRAY<T,FACE_INDEX<1> >& a)
{
    for(int i=a.domain_indices.min_corner.x;i<a.domain_indices.max_corner.x+1;i++)
        output_stream<<a.Component(0)(i)<<" ";
    output_stream<<std::endl;
    return output_stream;
}

template<class T,int d> inline std::ostream& operator<<(std::ostream& output_stream,const ARRAY<T,FACE_INDEX<d> >& a)
{
    for(int i=0;i<d;i++) output_stream<<a.data(i)<<std::endl;
    return output_stream;
}

extern template ARRAY<bool,FACE_INDEX<1> >::~ARRAY();
extern template ARRAY<bool,FACE_INDEX<2> >::~ARRAY();
extern template ARRAY<bool,FACE_INDEX<3> >::~ARRAY();
extern template ARRAY<int,FACE_INDEX<1> >::~ARRAY();
extern template ARRAY<int,FACE_INDEX<2> >::~ARRAY();
extern template ARRAY<int,FACE_INDEX<3> >::~ARRAY();
extern template ARRAY<float,FACE_INDEX<1> >::~ARRAY();
extern template ARRAY<float,FACE_INDEX<2> >::~ARRAY();
extern template ARRAY<float,FACE_INDEX<3> >::~ARRAY();
extern template ARRAY<double,FACE_INDEX<1> >::~ARRAY();
extern template ARRAY<double,FACE_INDEX<2> >::~ARRAY();
extern template ARRAY<double,FACE_INDEX<3> >::~ARRAY();

extern template ARRAY<VECTOR<double,1>, FACE_INDEX<1> >::~ARRAY();
extern template ARRAY<VECTOR<double,2>, FACE_INDEX<2> >::~ARRAY();
extern template ARRAY<VECTOR<double,3>, FACE_INDEX<3> >::~ARRAY();
extern template ARRAY<VECTOR<float,1>, FACE_INDEX<1> >::~ARRAY();
extern template ARRAY<VECTOR<float,2>, FACE_INDEX<2> >::~ARRAY();
extern template ARRAY<VECTOR<float,3>, FACE_INDEX<3> >::~ARRAY();

extern template ARRAY<VECTOR<double,3>, FACE_INDEX<1> >::~ARRAY();
extern template ARRAY<VECTOR<double,4>, FACE_INDEX<2> >::~ARRAY();
extern template ARRAY<VECTOR<double,5>, FACE_INDEX<3> >::~ARRAY();
extern template ARRAY<VECTOR<float,3>, FACE_INDEX<1> >::~ARRAY();
extern template ARRAY<VECTOR<float,4>, FACE_INDEX<2> >::~ARRAY();
extern template ARRAY<VECTOR<float,5>, FACE_INDEX<3> >::~ARRAY();

extern template void ARRAY<bool, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const bool&);
extern template void ARRAY<bool, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const bool&);
extern template void ARRAY<bool, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const bool&);
extern template void ARRAY<int, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const int&);
extern template void ARRAY<int, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const int&);
extern template void ARRAY<int, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const int&);
extern template void ARRAY<float, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const float&);
extern template void ARRAY<float, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const float&);
extern template void ARRAY<float, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const float&);
extern template void ARRAY<double, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const double&);
extern template void ARRAY<double, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const double&);
extern template void ARRAY<double, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const double&);

extern template void ARRAY<VECTOR<double,1>, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const VECTOR<double,1>&);
extern template void ARRAY<VECTOR<double,2>, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const VECTOR<double,2>&);
extern template void ARRAY<VECTOR<double,3>, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const VECTOR<double,3>&);
extern template void ARRAY<VECTOR<float,1>, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const VECTOR<float,1>&);
extern template void ARRAY<VECTOR<float,2>, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const VECTOR<float,2>&);
extern template void ARRAY<VECTOR<float,3>, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const VECTOR<float,3>&);

extern template void ARRAY<VECTOR<double,3>, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const VECTOR<double,3>&);
extern template void ARRAY<VECTOR<double,4>, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const VECTOR<double,4>&);
extern template void ARRAY<VECTOR<double,5>, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const VECTOR<double,5>&);
extern template void ARRAY<VECTOR<float,3>, FACE_INDEX<1> >::Resize(const RANGE<VECTOR<int, 1> >&, bool, bool, const VECTOR<float,3>&);
extern template void ARRAY<VECTOR<float,4>, FACE_INDEX<2> >::Resize(const RANGE<VECTOR<int, 2> >&, bool, bool, const VECTOR<float,4>&);
extern template void ARRAY<VECTOR<float,5>, FACE_INDEX<3> >::Resize(const RANGE<VECTOR<int, 3> >&, bool, bool, const VECTOR<float,5>&);

extern template void ARRAY<int, FACE_INDEX<1> >::Read<float>(std::istream&);
extern template void ARRAY<int, FACE_INDEX<2> >::Read<float>(std::istream&);
extern template void ARRAY<int, FACE_INDEX<3> >::Read<float>(std::istream&);
extern template void ARRAY<float, FACE_INDEX<1> >::Read<float>(std::istream&);
extern template void ARRAY<float, FACE_INDEX<2> >::Read<float>(std::istream&);
extern template void ARRAY<float, FACE_INDEX<3> >::Read<float>(std::istream&);
extern template void ARRAY<double, FACE_INDEX<1> >::Read<float>(std::istream&);
extern template void ARRAY<double, FACE_INDEX<2> >::Read<float>(std::istream&);
extern template void ARRAY<double, FACE_INDEX<3> >::Read<float>(std::istream&);
extern template void ARRAY<int, FACE_INDEX<1> >::Read<double>(std::istream&);
extern template void ARRAY<int, FACE_INDEX<2> >::Read<double>(std::istream&);
extern template void ARRAY<int, FACE_INDEX<3> >::Read<double>(std::istream&);
extern template void ARRAY<float, FACE_INDEX<1> >::Read<double>(std::istream&);
extern template void ARRAY<float, FACE_INDEX<2> >::Read<double>(std::istream&);
extern template void ARRAY<float, FACE_INDEX<3> >::Read<double>(std::istream&);
extern template void ARRAY<double, FACE_INDEX<1> >::Read<double>(std::istream&);
extern template void ARRAY<double, FACE_INDEX<2> >::Read<double>(std::istream&);
extern template void ARRAY<double, FACE_INDEX<3> >::Read<double>(std::istream&);

extern template void ARRAY<int, FACE_INDEX<1> >::Write<float>(std::ostream&) const;
extern template void ARRAY<int, FACE_INDEX<2> >::Write<float>(std::ostream&) const;
extern template void ARRAY<int, FACE_INDEX<3> >::Write<float>(std::ostream&) const;
extern template void ARRAY<float, FACE_INDEX<1> >::Write<float>(std::ostream&) const;
extern template void ARRAY<float, FACE_INDEX<2> >::Write<float>(std::ostream&) const;
extern template void ARRAY<float, FACE_INDEX<3> >::Write<float>(std::ostream&) const;
extern template void ARRAY<double, FACE_INDEX<1> >::Write<float>(std::ostream&) const;
extern template void ARRAY<double, FACE_INDEX<2> >::Write<float>(std::ostream&) const;
extern template void ARRAY<double, FACE_INDEX<3> >::Write<float>(std::ostream&) const;
extern template void ARRAY<int, FACE_INDEX<1> >::Write<double>(std::ostream&) const;
extern template void ARRAY<int, FACE_INDEX<2> >::Write<double>(std::ostream&) const;
extern template void ARRAY<int, FACE_INDEX<3> >::Write<double>(std::ostream&) const;
extern template void ARRAY<float, FACE_INDEX<1> >::Write<double>(std::ostream&) const;
extern template void ARRAY<float, FACE_INDEX<2> >::Write<double>(std::ostream&) const;
extern template void ARRAY<float, FACE_INDEX<3> >::Write<double>(std::ostream&) const;
extern template void ARRAY<double, FACE_INDEX<1> >::Write<double>(std::ostream&) const;
extern template void ARRAY<double, FACE_INDEX<2> >::Write<double>(std::ostream&) const;
extern template void ARRAY<double, FACE_INDEX<3> >::Write<double>(std::ostream&) const;
}
#endif
