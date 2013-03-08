//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <PhysBAM_Tools/Log/LOG.h>
#include "DELAUNAY_TRIANGULATION_2D.h"
namespace PhysBAM{
//#####################################################################
// Function Triangulate
//#####################################################################
template<class T> void DELAUNAY_TRIANGULATION_2D<T>::
Triangulate(const ARRAY_VIEW<TV>& X,ARRAY<E>& elements,const T threshold)
{
    std::ofstream ptsf("pts");
    if(ptsf.is_open()){
        ptsf<<"2\n";
        ptsf<<X.m<<"\n";
        for(int i=0;i<X.m;i++) ptsf<<X(i).x<<" "<<X(i).y<<"\n";
        ptsf.close();}
    if(system("qdelaunay i < pts > tris")==0){
        LOG::cout<<"Delaunay triangulation done."<<std::endl;
        elements.Clean_Memory();
        std::string line;
        std::istringstream ss;
        std::ifstream trisf("tris");
        if(trisf.is_open()){
            int line_number=1;
            while(getline(trisf,line)){
                ss.clear();
                ss.str(line);
                int a,b,c;
                if(line_number++!=1){
                    ss>>a>>b>>c;
                    elements.Append(E(a,b,c));}}
            trisf.close();}}
    else LOG::cout<<"Delaunay triangulation failed."<<std::endl;
    ARRAY<E> filted_elements;
    T threshold2=sqr(threshold);
    for(int i=0;i<elements.m;i++){
        TV A=X(elements(i)(0)),B=X(elements(i)(1)),C=X(elements(i)(2));
        if((A-B).Magnitude_Squared()>threshold2 || (B-C).Magnitude_Squared()>threshold2 || (C-A).Magnitude_Squared()>threshold2){}
        else filted_elements.Append(elements(i));}
    elements=filted_elements;
}
//#####################################################################
template class DELAUNAY_TRIANGULATION_2D<float>;
template class DELAUNAY_TRIANGULATION_2D<double>;
}
