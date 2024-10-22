//#####################################################################
// Copyright 2009-2010, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __Relaxation_Interior_Helper__
#define __Relaxation_Interior_Helper__

namespace PhysBAM{

template<class T>
class Relaxation_Interior_Size_Specific_Helper_Base
{
protected:
    T* const u;
    const T* const b;
    T* const delta;
    const unsigned char* const bit_writemask;

public:
    explicit Relaxation_Interior_Size_Specific_Helper_Base(T* const u_input,const T* const b_input,T* const delta_input,
        const unsigned char* const bit_writemask_input)
        :u(u_input),b(b_input),delta(delta_input),bit_writemask(bit_writemask_input)
    {}

    virtual ~Relaxation_Interior_Size_Specific_Helper_Base() {}

//#####################################################################
    virtual void Run()=0;
    virtual void Run_Parallel(const int number_of_partitions)=0;
//#####################################################################
};

template<class T>
class Relaxation_Interior_Helper{
    const int x_size,y_size,z_size;
    Relaxation_Interior_Size_Specific_Helper_Base<T>* derived;

public:
    ~Relaxation_Interior_Helper() {delete derived;}
    void Run() {derived->Run();}
    void Run_Parallel(const int number_of_partitions) {derived->Run_Parallel(number_of_partitions);}

//#####################################################################
    Relaxation_Interior_Helper(const int x_size_input,const int y_size_input,const int z_size_input,T* const u,const T* const b,T* const delta,
        const unsigned char* const bit_writemask);
//#####################################################################
};

template<class T,int y_size,int z_size=y_size>
class Relaxation_Interior_Size_Specific_Helper:public Relaxation_Interior_Size_Specific_Helper_Base<T>
{
    typedef Relaxation_Interior_Size_Specific_Helper_Base<T> Base;
    using Base::u;using Base::b;using Base::delta;using Base::bit_writemask;

    const int x_size,padded_x_size;
    const int coarse_x_size,coarse_padded_x_size;
    
    enum WORKAROUND{
        x_block_size=4,
        y_block_size=4,
        z_block_size=4,
        padded_y_size=y_size+2,
        padded_z_size=z_size+2,
        coarse_y_size=y_size/2,
        coarse_z_size=z_size/2,
        coarse_padded_y_size=coarse_y_size+2,
        coarse_padded_z_size=coarse_z_size+2,
        x_shift=padded_y_size*padded_z_size,
        y_shift=padded_z_size,
        z_shift=1,
        coarse_x_shift=coarse_padded_y_size*coarse_padded_z_size,
        coarse_y_shift=coarse_padded_z_size,
        coarse_z_shift=1,
        x_plus_one_z_minus_one_shift=x_shift-z_shift,
        x_plus_one_z_plus_two_shift=x_shift+2*z_shift,
        x_plus_one_y_plus_two_z_plus_one_shift=x_shift+2*y_shift+z_shift,
        x_plus_one_y_plus_one_z_minus_one_shift=x_shift+y_shift-z_shift,
        x_plus_one_y_plus_one_z_plus_two_shift=x_shift+y_shift+2*z_shift,
        x_plus_one_y_plus_two_shift=x_shift+2*y_shift,
        x_plus_two_y_plus_one_z_plus_one_shift=2*x_shift+y_shift+z_shift,
        x_plus_two_y_plus_one_shift=2*x_shift+y_shift,
        x_plus_one_y_minus_one_z_plus_one_shift=x_shift-y_shift+z_shift,
        x_plus_one_y_minus_one_shift=x_shift-y_shift,
        z_plus_two_shift=2*z_shift,
        x_plus_one_y_plus_one_z_plus_one_shift=x_shift+y_shift+z_shift,
        x_plus_one_y_plus_one_shift=x_shift+y_shift,
        z_minus_one_shift=-z_shift,
        y_minus_one_z_plus_one_shift=-y_shift+z_shift,
        x_plus_one_z_plus_one_shift=x_shift+z_shift,
        x_plus_one_shift=x_shift,
        y_minus_one_shift=-y_shift,
        y_plus_two_z_plus_one_shift=2*y_shift+z_shift,        
        y_plus_two_shift=2*y_shift,
        x_plus_two_shift=2*x_shift,
        x_plus_two_z_plus_one_shift=2*x_shift+z_shift,
        y_plus_one_z_plus_two_shift=y_shift+2*z_shift,
        y_plus_one_z_minus_one_shift=y_shift-z_shift,
        x_minus_one_z_plus_one_shift=-x_shift+z_shift,
        y_plus_one_shift=y_shift,
        x_minus_one_y_plus_one_shift=-x_shift+y_shift,
        x_minus_one_shift=-x_shift,
        z_plus_one_shift=z_shift,
        x_minus_one_y_plus_one_z_plus_one_shift=-x_shift+y_shift+z_shift,
        y_plus_one_z_plus_one_shift=y_shift+z_shift
    };

public:
    explicit Relaxation_Interior_Size_Specific_Helper(const int x_size_input,T* const u_input,const T* const b_input,T* const delta_input,const unsigned char* const bit_writemask_input)
        :Base(u_input,b_input,delta_input,bit_writemask_input)
        ,x_size(x_size_input),padded_x_size(x_size_input+2)
        ,coarse_x_size(x_size/2),coarse_padded_x_size(coarse_x_size+2)
    {}

    void Run()
    {Run_Major_X_Range(1,coarse_x_size);Run_Minor_X_Range(1,coarse_x_size);}
    
    // For debugging purposes only

//     static void Allocate_Data(T*& u,T*& b,T*& r,unsigned char*& bit_writemask)
//     {u=new T[padded_length];b=new T[padded_length];r=new T[padded_length];bit_writemask=new unsigned char[coarse_padded_length];}

//     static void Initialize_Data(T* const u,T* const b,T* const r,unsigned char* const bit_writemask)
//     {for(int i=0;i<padded_length;i++) u[i]=b[i]=r[i]=(T)i;for(int i=0;i<coarse_padded_length;i++) bit_writemask[i]=0xff;}

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    // compute deltas for all indices in range. apply deltas to xmin to xmax-1
    void Run_Major_X_Range(const int xmin,const int xmax);
    // apply deltas to xmax
    void Run_Minor_X_Range(const int xmin,const int xmax);
//#####################################################################
};
}
#endif
