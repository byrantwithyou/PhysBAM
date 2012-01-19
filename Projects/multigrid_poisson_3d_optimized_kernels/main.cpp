//#####################################################################
// Copyright 2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include "FINITE_DIFFERENCE_DISCRETIZATION.h"
#include "STENCIL_ITERATOR.h"
using namespace PhysBAM;

struct Always_True_Functor{
    bool operator()(const VECTOR<int,3>& index){return true;}
    bool operator()(const VECTOR<int,2>& index){return true;}
} always_true;

template<class T,int d>
class Stencil_Implementation_Helper{
    typedef VECTOR<int,d> T_INDEX;
    typedef STENCIL<T,d> T_STENCIL;
    typedef STENCIL_ITERATOR<const T,d> T_CONST_STENCIL_ITERATOR;

public:

    static void Print_Enumerators(std::ostream& output,const T_STENCIL& stencil)
    {
        bool on_first_line_of_output=true;
        for(T_CONST_STENCIL_ITERATOR iterator(stencil);iterator.Valid();iterator.Next()){
            const T_INDEX& index=iterator.Key();
            if(index==T_INDEX()) continue;
            if(!on_first_line_of_output) output<<","<<std::endl;
            else on_first_line_of_output=false;
            output<<"    "<<Shift_Name_From_Index(index)<<"="<<Shift_Value_From_Index(index);}
        if(!on_first_line_of_output) output<<std::endl;
    }

    static void Print_Operation(std::ostream& output,const T_STENCIL& stencil)
    {
        HASHTABLE<T,ARRAY<T_INDEX> > hash;
        for(T_CONST_STENCIL_ITERATOR iterator(stencil);iterator.Valid();iterator.Next())
            hash.Get_Or_Insert(iterator.Data(),ARRAY<T_INDEX>()).Append(iterator.Key());
        bool on_first_coefficient=true;
        for(HASHTABLE_ITERATOR<T,ARRAY<T_INDEX> > iterator(hash);iterator.Valid();iterator.Next()){
            output<<"    result";
            if(!on_first_coefficient) output<<"+";else on_first_coefficient=false;
            output<<"=(T)"<<iterator.Key()<<"*(";
            const ARRAY<T_INDEX>& indices=iterator.Data();
            for(int i=0;i<indices.m;i++){
                output<<std::endl<<"        ";
                if(i>1) output<<"+input[index"; else output<<"input[index";
                if(indices(i)!=T_INDEX()) output<<"+"<<Shift_Name_From_Index(indices(i));
                output<<"]";}
            output<<");"<<std::endl;}
    }

private:
    static std::string Coordinate_Name_From_Index(const int index)
    {
        switch(index){
            case 1: PHYSBAM_ASSERT(d>=1);return std::string("x");
            case 2: PHYSBAM_ASSERT(d>=2);return std::string("y");
            case 3: PHYSBAM_ASSERT(d==3);return std::string("z");}
        PHYSBAM_FATAL_ERROR();
    }

    static std::string Offset_Name_From_Index(const int index)
    {
        switch(index){
            case -4: return std::string("minus_four");
            case -3: return std::string("minus_three");
            case -2: return std::string("minus_two");
            case -1: return std::string("minus_one");
            case 1: return std::string("plus_one");
            case 2: return std::string("plus_two");
            case 3: return std::string("plus_three");
            case 4: return std::string("plus_four");}
        PHYSBAM_FATAL_ERROR();
    }

    static std::string Shift_Name_From_Index(const T_INDEX& index)
    {
        std::string shift_name;
        for(int i=0;i<d;i++){
            int j=index(i);
            if(j==0) continue;
            if(shift_name!="") shift_name+="_";
            shift_name+=Coordinate_Name_From_Index(i)+"_"+Offset_Name_From_Index(j);}
        if(index!=T_INDEX()) shift_name+="_shift";
        return shift_name;
    }

    static std::string Shift_Value_From_Index(const T_INDEX& index)
    {
        std::string shift_value;
        for(int i=0;i<d;i++){
            int j=index(i);
            if(j==0) continue;
            if(shift_value!="" && j>0) shift_value+="+";
            if(j==-1) shift_value+="-";
            else if(j!=1) shift_value+=STRING_UTILITIES::Value_To_String(j)+"*";
            shift_value+=Coordinate_Name_From_Index(i)+"_shift";}
        return shift_value;
    }
};

int main(int argc,char *argv[])
{
    typedef float T;
    static const int d=3;
    typedef VECTOR<int,d> T_INDEX;
    typedef STENCIL<T,d> T_STENCIL;
    typedef STENCIL_ITERATOR<T,d> T_STENCIL_ITERATOR;

    T_STENCIL residual;
    for(int v=0;v<d;v++){
        residual.Insert(T_INDEX::Axis_Vector(v),1);
        residual.Insert(-T_INDEX::Axis_Vector(v),1);}
    residual.Insert(T_INDEX(),3);
//     std::cout<<"== Residual stencil =="<<std::endl<<std::endl;
//     Print_Stencil(std::cout,residual);
//     std::cout<<std::endl;

    FINITE_DIFFERENCE_DISCRETIZATION<T,d> discretization((T)1);
    T_STENCIL interpolation;
    discretization.Set_Offset_Multilinear_Interpolation_Stencil(T_INDEX(),interpolation,(T)64,always_true,always_true);
//     std::cout<<"== Interpolation stencil =="<<std::endl<<std::endl;
//     Print_Stencil(std::cout,interpolation);
//     std::cout<<std::endl;

    T_STENCIL product=residual.Convolve(interpolation);
//     std::cout<<"== Combined stencil =="<<std::endl<<std::endl;
//     Print_Stencil(std::cout,product);
//     std::cout<<std::endl;

    std::cout<<std::endl<<"Enumerators:"<<std::endl<<std::endl;
    Stencil_Implementation_Helper<T,d>::Print_Enumerators(std::cout,product);
    std::cout<<std::endl<<"Operation:"<<std::endl<<std::endl;
    Stencil_Implementation_Helper<T,d>::Print_Operation(std::cout,product);

    return 0;
}
