//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ORIENTED_BOX
//##################################################################### 
#include <Core/Log/LOG.h>
#include <Tools/Auto_Diff/AUTO_DIFF.h>
#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
namespace PhysBAM{
//#####################################################################
// Function Intersection_Helper
//#####################################################################
template<class T,class T_BOX> static bool Intersection_Helper(const ORIENTED_BOX<VECTOR<T,1> >& self,const T_BOX& box,const VECTOR<T,1>& u,const VECTOR<T,1>& v){PHYSBAM_NOT_IMPLEMENTED();}
template<class T,class T_BOX> static bool Intersection_Helper(const ORIENTED_BOX<VECTOR<T,2> >& self,const T_BOX& box,const VECTOR<T,2>& u,const VECTOR<T,2>& v){return false;}
template<class T,class T_BOX> static bool Intersection_Helper(const ORIENTED_BOX<VECTOR<T,3> >& self,const T_BOX& box,const VECTOR<T,3>& u,const VECTOR<T,3>& v)
{
    return self.Separating_Test(box,VECTOR<T,3>::Cross_Product(u,v));
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class TV> bool ORIENTED_BOX<TV>::
Intersection(const ORIENTED_BOX& box) const
{
    for(int i=0;i<d;i++) if(Separating_Test(box,edges.Column(i)) || Separating_Test(box,box.edges.Column(i))) return false;
    for(int i=0;i<d;i++) for(int j=0;j<d;j++) if(Intersection_Helper(*this,box,edges.Column(i),box.edges.Column(j))) return false;
    return true; // otherwise
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class TV> bool ORIENTED_BOX<TV>::
Intersection(const RANGE<TV>& box) const
{
    for(int j=0;j<d;j++){
        T line_min=corner(j),line_max=corner(j);
        for(int i=0;i<d;i++) if(edges.Column(i)(j)>0) line_max+=edges.Column(i)(j);else line_min+=edges.Column(i)(j);
        if(line_max<box.min_corner(j) || line_min>box.max_corner(j)) return false;}
    for(int i=0;i<d;i++) if(Separating_Test(box,edges.Column(i))) return false;
    if(d==3) for(int j=0;j<d;j++){
        TV box_edge;box_edge(j)=box.max_corner(j)-box.min_corner(j);
        for(int i=0;i<d;i++) if(Intersection_Helper(*this,box,edges.Column(i),box_edge)) return false;}
    return true; // otherwise
}
//#####################################################################
// Function Separating_Test
//#####################################################################
template<class TV> inline  bool ORIENTED_BOX<TV>::
Separating_Test(const ORIENTED_BOX& box,const TV& direction) const
{
    T min1,max1;Project_Points_Onto_Line(direction,min1,max1);
    T min2,max2;box.Project_Points_Onto_Line(direction,min2,max2);
    if(max2<min1 || min2>max1) return true;else return false;
}
//#####################################################################
// Function Separating_Test
//#####################################################################
template<class TV> inline bool ORIENTED_BOX<TV>::
Separating_Test(const RANGE<TV>& box,const TV& direction) const
{
    T min1,max1;Project_Points_Onto_Line(direction,min1,max1);
    T min2,max2;box.Project_Points_Onto_Line(direction,min2,max2);
    if(max2<min1 || min2>max1) return true;else return false;
}
//#####################################################################
// Function Project_Points_Onto_Line
//#####################################################################
template<class TV> inline void ORIENTED_BOX<TV>::
Project_Points_Onto_Line(const TV& direction,T& line_min,T& line_max) const
{
    line_min=line_max=TV::Dot_Product(direction,corner);
    TV e=edges.Transpose_Times(direction);
    for(int i=0;i<d;i++) if(e(i)>0) line_max+=e(i);else line_min+=e(i);
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class TV> typename TV::SCALAR ORIENTED_BOX<TV>::
Signed_Distance(const TV& X) const
{
    TV lengths=edges.Column_Magnitudes();
    if(lengths.Contains(0)) PHYSBAM_NOT_IMPLEMENTED();
    TV phi=abs(edges.Transpose_Times(X-Center()))/lengths-(T).5*lengths;
    if(!phi.All_Less_Equal(TV())) return TV::Componentwise_Max(phi,TV()).Magnitude();
    return phi.Max();
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV ORIENTED_BOX<TV>::
Normal(const TV& X) const
{
    TV lengths=edges.Column_Magnitudes();
    if(lengths.Contains(0)) PHYSBAM_NOT_IMPLEMENTED();
    MATRIX<T,TV::m> M=DIAGONAL_MATRIX<T,TV::m>(lengths).Inverse()*edges.Transposed();
    AUTO_DIFF<TV,TV> t(M*X,M),phi=abs(t-M*Center())-(T).5*lengths;
    if(phi.x.All_Less_Equal(TV())) return max(phi).dx;
    for(int i=0;i<TV::m;i++) if(phi.x(i)<0) phi.Set_Entry(i,AUTO_DIFF<T,TV>());
    return phi.Magnitude().dx;
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> ORIENTED_BOX<TV>::
Hessian(const TV& X) const
{
    TV lengths=edges.Column_Magnitudes();
    if(lengths.Contains(0)) PHYSBAM_NOT_IMPLEMENTED();
    MATRIX<T,TV::m> M=DIAGONAL_MATRIX<T,TV::m>(lengths).Inverse()*edges.Transposed();
    AUTO_HESS<TV,TV> t(M*X,M),phi=abs(t-M*Center())-(T).5*lengths;
    if(phi.x.All_Less_Equal(TV())) return max(phi).ddx;
    for(int i=0;i<TV::m;i++) if(phi.x(i)<0) phi.Set_Entry(i,AUTO_HESS<T,TV>());
    return phi.Magnitude().ddx;
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string ORIENTED_BOX<TV>::
Name()
{   
    return LOG::sprintf("ORIENTED_BOX<VECTOR<T,%d>",d);
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template std::string ORIENTED_BOX<VECTOR<T,d> >::Name(); \
    template bool ORIENTED_BOX<VECTOR<T,d> >::Intersection(const ORIENTED_BOX<VECTOR<T,d> >& box) const; \
    template bool ORIENTED_BOX<VECTOR<T,d> >::Intersection(const RANGE<VECTOR<T,d> >& box) const; \
    template bool ORIENTED_BOX<VECTOR<T,d> >::Separating_Test(const ORIENTED_BOX& box,const VECTOR<T,d> & direction) const; \
    template bool ORIENTED_BOX<VECTOR<T,d> >::Separating_Test(const RANGE<VECTOR<T,d> >& box,const VECTOR<T,d> & direction) const; \
    template void ORIENTED_BOX<VECTOR<T,d> >::Project_Points_Onto_Line(const VECTOR<T,d> & direction,T& line_min,T& line_max) const; \
    template T ORIENTED_BOX<VECTOR<T,d> >::Signed_Distance(const VECTOR<T,d>& X) const; \
    template SYMMETRIC_MATRIX<T,d> ORIENTED_BOX<VECTOR<T,d> >::Hessian(const VECTOR<T,d>&) const; \
    template VECTOR<T,d> ORIENTED_BOX<VECTOR<T,d> >::Normal(const VECTOR<T,d>& X) const;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
}
