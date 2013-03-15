//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "DELAUNAY_TRIANGULATION_2D.h"
namespace PhysBAM{
//#####################################################################
// Function Triangulate
//#####################################################################
template<class T> void DELAUNAY_TRIANGULATION_2D<T>::
Triangulate(const ARRAY_VIEW<TV>& X,TRIANGULATED_AREA<T>& ta,const T length_max,const T theta_min)
{
    std::ofstream ptsf("pts");
    if(ptsf.is_open()){
        ptsf<<"2\n";
        ptsf<<X.m<<"\n";
        for(int i=0;i<X.m;i++) ptsf<<X(i).x<<" "<<X(i).y<<"\n";
        ptsf.close();}
    ARRAY<E> elements;
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
    if(system("rm pts")) PHYSBAM_FATAL_ERROR();
    if(system("rm tris")) PHYSBAM_FATAL_ERROR();
    ARRAY<E> filted_elements;
    T length_max2=sqr(length_max);
    T theta_min_converted=theta_min/(T)360*(T)6.28318530717959;
    T cos_max=cos(theta_min_converted);LOG::cout<<cos_max<<std::endl;
    for(int i=0;i<elements.m;i++){
        TV A=X(elements(i)(0)),B=X(elements(i)(1)),C=X(elements(i)(2));
        T a=(B-C).Magnitude(),b=(A-C).Magnitude(),c=(A-B).Magnitude();
        T cosA=(b*b+c*c-a*a)/((T)2*b*c),cosB=(a*a+c*c-b*b)/((T)2*a*c),cosC=(b*b+a*a-c*c)/((T)2*b*a);
        bool good=true;
        if((A-B).Magnitude_Squared()>length_max2 || (B-C).Magnitude_Squared()>length_max2 || (C-A).Magnitude_Squared()>length_max2) good=false;
        if(abs(cosA)>cos_max || abs(cosB)>cos_max || abs(cosC)>cos_max) good=false;
        if(good) filted_elements.Append(elements(i));}
    ta.Clean_Memory();
    ta.particles.Delete_All_Elements();
    ta.particles.Add_Elements(X.m);
    for(int i=0;i<X.m;i++) ta.particles.X(i)=X(i);
    ta.mesh.number_nodes=X.m;
    ta.mesh.elements.Exact_Resize(filted_elements.m);
    for(int i=0;i<filted_elements.m;i++) ta.mesh.elements(i).Set(filted_elements(i).x,filted_elements(i).y,filted_elements(i).z);
    Test(ta);
}
//#####################################################################
// Function Test
//#####################################################################
template<class T> void DELAUNAY_TRIANGULATION_2D<T>::
Test(const TRIANGULATED_AREA<T>& ta)
{
    ARRAY<bool> covered(ta.particles.number);
    covered.Fill(false);
    for(int i=0;i<ta.mesh.elements.m;i++)
        for(int k=0;k<3;k++)
            covered(ta.mesh.elements(i)(k))=true;
    for(int i=0;i<covered.m;i++)
        if(!covered(i)){
            LOG::cout<<"Error: Particle "<<i<<" is not covered by Delaunay Triangulation. (Deleted by threshold?)";
            PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
template class DELAUNAY_TRIANGULATION_2D<float>;
template class DELAUNAY_TRIANGULATION_2D<double>;
}
