#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <fstream>
#include <iostream>

using namespace PhysBAM;

bool verbose=false;
bool opt_double=false,opt_float=false,opt_int=false,opt_bool=false,opt_opengl_color=false,opt_array=false;
bool opt_arrays=false,opt_1d=false,opt_2d=false,opt_3d=false,opt_vec2d=false,opt_vec3d=false,opt_vec3d=false;
int num_arrays=1;

template<template<class> class ARRAY,class T> void Print_Min_Value(const ARRAY<T> &array)
{ if(array.m>0) std::cout << "Min = " << ARRAY<T>::min(array) << std::endl; }

template<template<class> class ARRAY,class T> void Print_Max_Value(const ARRAY<T> &array)
{ if(array.m>0) std::cout << "Max = " << ARRAY<T>::max(array) << std::endl; }

template<template<class> class ARRAY,class T> void Print_Max_Value(const ARRAY<VECTOR<T,2> > &array)
{
    T max_value = ARRAY<VECTOR<T,2> >::template Maximum_Vector_Magnitude<T>(array);
    std::cout << "Max = " << max_value << std::endl;
}

template<template<class> class ARRAY,class T> void Print_Max_Value(const ARRAY<VECTOR<T,3> > &array)
{ 
    T max_value = ARRAY<VECTOR<T,3> >::template Maximum_Vector_Magnitude<T>(array);
    std::cout << "Max = " << max_value << std::endl;
}

template<> void Print_Max_Value(const ARRAY<VECTOR<float,3> > &array)
{ 
    float max_magnitude_squared=-1;int max_index=0;
    for(int i=0;i<array.m;i++){float magnitude_squared=array(i).Magnitude_Squared();
        if(magnitude_squared>max_magnitude_squared){max_magnitude_squared=magnitude_squared;max_index=i;}}
    if(max_index) std::cout << "Max = " << sqrt(max_magnitude_squared) << " (at " << max_index << ")" << std::endl;
}

template<template<class> class ARRAY,class T> void Print_Min_Value(const ARRAY<VECTOR<T,2> > &array) {}
template<template<class> class ARRAY,class T> void Print_Min_Value(const ARRAY<VECTOR<T,3> > &array) {}
template<template<class> class ARRAY> void Print_Min_Value(const ARRAY<OPENGL_COLOR> &array) {}
template<template<class> class ARRAY> void Print_Max_Value(const ARRAY<VECTOR<int,2> > &array) {}
template<template<class> class ARRAY> void Print_Max_Value(const ARRAY<VECTOR<bool,2> > &array) {}
template<template<class> class ARRAY> void Print_Max_Value(const ARRAY<VECTOR<int,3> > &array) {}
template<template<class> class ARRAY> void Print_Max_Value(const ARRAY<VECTOR<bool,3> > &array) {}
template<template<class> class ARRAY> void Print_Max_Value(const ARRAY<OPENGL_COLOR> &array) {}

template<class T,class RW> void Print_Array(std::istream &input, int num_arrays)
{
    ARRAY<T> array;

    for (int t=1;t<=num_arrays;t++)
    {
        array.template Read<RW>(input);
        printf("Array %d:\n", t);
        printf("m=%d\n", array.m);
        Print_Min_Value(array); Print_Max_Value(array);
        if (verbose) {
            for (int i=1;i<=array.m;i++)
                std::cout << "array(" << i << ") = " << array(i) << std::endl;
        }
        std::cout << std::endl;
    }
}

template<class T,class RW> void Print_Arrays(std::istream &input, int num_arrays)
{
    ARRAYS<T> array;

    for (int t=1;t<=num_arrays;t++)
    {
        array.template Read<RW>(input);
        printf("Array %d:\n", t);
        printf("length=%d m=%d\n", array.length, array.m);
        //Print_Min_Value(array); Print_Max_Value(array);
        if (verbose) {
            for (int i=1;i<=array.m;i++)
                for(int t=0;t<array.length;t++)
                    std::cout << "array(" << t << "," << i << ") = " << array(t,i) << std::endl;
        }
        std::cout << std::endl;
    }
}


template<class T,class RW> void Print_Array_1D(std::istream &input, int num_arrays)
{
    ARRAYS<VECTOR<T,1> > array;

    for (int t=1;t<=num_arrays;t++)
    {
        array.template Read<RW>(input);
        printf("Array %d:\n", t);
        printf("length=%d m_start=%d m_end=%d\n", array.length, array.m_start, array.m_end);
        Print_Min_Value(array); Print_Max_Value(array);
        if (verbose) {
            if (array.length == 1) {
                for (int i=array.m_start;i<=array.m_end;i++)
                    std::cout << "array(" << i << ") = " << array(i) << std::endl;
            } else {
                for (int i=array.m_start;i<=array.m_end;i++)
                    for (int t=1;t<=array.length;t++)
                        std::cout << "array(" << t << "," << i << ") = " << array(t,i) << std::endl;
            }
        }
        std::cout << std::endl;
    }
}

template<class T,class RW> void Print_Array_2D(std::istream &input, int num_arrays)
{
    ARRAYS<VECTOR<T,2> > array;

    for (int t=1;t<=num_arrays;t++)
    {
        array.template Read<RW>(input);
        printf("Array %d:\n", t);
        printf("length=%d m_start=%d m_end=%d n_start=%d n_end=%d\n", array.length, array.m_start, array.m_end, array.n_start, array.n_end);
        Print_Min_Value(array); Print_Max_Value(array);
        if (verbose) {
            if (array.length == 1) {
                for (int i=array.m_start;i<=array.m_end;i++)
                    for (int j=array.n_start;j<=array.n_end;j++)
                        std::cout << "array(" << i << "," << j << ") = " << array(i,j) << std::endl;
            } else {
                for (int i=array.m_start;i<=array.m_end;i++)
                    for (int j=array.n_start;j<=array.n_end;j++)
                        for (int t=1;t<=array.length;t++)
                            std::cout << "array(" << t << "," << i << "," << j << ") = " << array(t,i,j) << std::endl;
            }
        }
        std::cout << std::endl;
    }
}

template<class T,class RW> void Print_Array_3D(std::istream &input, int num_arrays)
{
    ARRAYS<VECTOR<T,3> > array;

    for (int t=1;t<=num_arrays;t++)
    {
        array.template Read<RW>(input);
        printf("Array %d:\n", t);
        printf("length=%d m_start=%d m_end=%d n_start=%d n_end=%d mn_start=%d mn_end=%d\n",
               array.length, array.m_start, array.m_end, array.n_start, array.n_end, array.mn_start, array.mn_end);
        Print_Min_Value(array); Print_Max_Value(array);
        if (verbose) {
            if (array.length == 1) {
                for (int i=array.m_start;i<=array.m_end;i++)
                    for (int j=array.n_start;j<=array.n_end;j++)
                        for (int k=array.mn_start;k<=array.mn_end;k++)
                            std::cout << "array(" << i << "," << j << "," << k << ") = " << array(i,j,k) << std::endl;
            } else {
                for (int i=array.m_start;i<=array.m_end;i++)
                    for (int j=array.n_start;j<=array.n_end;j++)
                        for (int k=array.mn_start;k<=array.mn_end;k++)
                            for (int t=1;t<=array.length;t++)
                                std::cout << "array(" << t << "," << i << "," << j << "," << k << ") = " << array(t,i,j,k) << std::endl;
            }
        }
        std::cout << std::endl;
    }
}

template<class T,class RW> void By_Composite_Type(std::istream &input, PARSE_ARGS &parse_args)
{
    if (opt_3d) Print_Array_3D<T,RW>(input,num_arrays);
    else if (opt_2d) Print_Array_2D<T,RW>(input,num_arrays);
    else if (opt_1d) Print_Array_1D<T,RW>(input,num_arrays);
    else if (opt_arrays) Print_Arrays<T,RW>(input,num_arrays);
    else Print_Array<T,RW>(input,num_arrays);
}

template<class T> void By_Base_Type(std::istream &input, PARSE_ARGS &parse_args)
{
    if (opt_vec3d) By_Composite_Type<VECTOR<T,3>,T>(input, parse_args);
    else if (opt_vec2d) By_Composite_Type<VECTOR<T,2>,T>(input, parse_args);
    else By_Composite_Type<T,T>(input, parse_args);
}

int main(int argc, char *argv[])
{
    int skip=0,opt_n=1;

    std::string filename;
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-double",&opt_double,"double");
    parse_args.Add("-float",&opt_float,"float");
    parse_args.Add("-int",&opt_int,"int");
    parse_args.Add("-bool",&opt_bool,"bool");
    parse_args.Add("-opengl_color",&opt_opengl_color,"opengl_color");
    parse_args.Add("-array",&opt_array,"array");   // ARRAY
    parse_args.Add("-arrays",&opt_arrays,"arrays");   // ARRAYS
    parse_args.Add("-1d",&opt_1d,"1d");      // ARRAYS_1D
    parse_args.Add("-2d",&opt_2d,"2d");      // ARRAYS_2D
    parse_args.Add("-3d",&opt_3d,"3d");      // ARRAYS_3D
    parse_args.Add("-vec2d",&opt_vec2d,"vec2d");
    parse_args.Add("-vec3d",&opt_vec3d,"vec3d");
    parse_args.Add("-vec3d",&opt_vec3d,"vec3d");
    parse_args.Add("-v",&verbose,"verbose");
    parse_args.Add("-skip", &skip, "bytes", "skip header bytes");
    parse_args.Add("-n", &num_arrays, "num-arrays", "number of consecutive arrays in the file");
    parse_args.Extra(&filename,"filename","filename");
    parse_args.Parse();

    verbose=opt_v;

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    if(!(*input)) {
        std::cerr << "Cannot open " << filename << std::endl;
        return 1;
    }

    if (skip > 0) {
        std::cout << "Skipping " << skip << " bytes" << std::endl;
        input->ignore(skip);
    }

    if (opt_double) {
        std::cout << "Base type = double" << std::endl;
        By_Base_Type<double>(*input, parse_args);
    }
    else if (opt_float) {
        std::cout << "Base type = float" << std::endl;
        By_Base_Type<float>(*input, parse_args);
    }
    else if (opt_int) {
        std::cout << "Base type = int" << std::endl;
        By_Base_Type<int>(*input, parse_args);
    }
    else if (opt_opengl_color) {
        std::cout << "Base type = OPENGL_COLOR" << std::endl;
        By_Composite_Type<OPENGL_COLOR,float>(*input, parse_args);
    } else {
        std::cout << "Base type = bool" << std::endl;
        By_Base_Type<bool>(*input, parse_args);
    }

    delete input;
}
