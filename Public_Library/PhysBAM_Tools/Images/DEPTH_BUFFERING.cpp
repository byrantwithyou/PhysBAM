//#####################################################################
// Copyright 2012, Alexey, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Images/DEPTH_BUFFERING.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> DISPLAY_PRIMITIVE<T>::
DISPLAY_PRIMITIVE(const TV &a,const int style):style(style)
{ 
    type=POINT;
    vertices(0)=a;
    vertices(1)=a;
    vertices(2)=a;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> DISPLAY_PRIMITIVE<T>::
DISPLAY_PRIMITIVE(const TV &a,const TV &b,const int style):style(style)
{ 
    type=SEGMENT;
    vertices(0)=a;
    vertices(1)=b;
    vertices(2)=b;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> DISPLAY_PRIMITIVE<T>::
DISPLAY_PRIMITIVE(const TV &a,const TV &b,const TV &c,const int style):style(style)
{ 
    type=TRIANGLE;
    vertices(0)=a;
    vertices(1)=b;
    vertices(2)=c;
}
//#####################################################################
// Function Initialize_Bounding_Box
//#####################################################################
template<class T> void DISPLAY_PRIMITIVE<T>::
Initialize_Bounding_Box()
{
    bounding_box=RANGE<TV2>::Bounding_Box(
        DB::Project(vertices(0)),
        DB::Project(vertices(1)),
        DB::Project(vertices(2)));
}
//#####################################################################
// Function Initialize_Elements
//#####################################################################
template<class T> void DISPLAY_PRIMITIVE_CUTTING<T>::
Initialize_Elements()
{
    elements.Clean_Memory();
    elements.Append(vertices);
}
//#####################################################################
// Function Cut_By_Primitive
//#####################################################################
template<class T> void DISPLAY_PRIMITIVE_CUTTING<T>::
Cut_By_Primitive(const DISPLAY_PRIMITIVE_CUTTING<T> &p)
{
    if(p.type==TRIANGLE && type!=POINT && !RANGE<TV2>::Intersect(bounding_box,p.bounding_box).Empty()){
        Cut_By_Plane(DB::Get_Cutting_Plane(p.vertices(0),p.vertices(1)));
        Cut_By_Plane(DB::Get_Cutting_Plane(p.vertices(1),p.vertices(2)));
        Cut_By_Plane(DB::Get_Cutting_Plane(p.vertices(2),p.vertices(0)));
        Cut_By_Plane(PLANE<T>(p.vertices(0),p.vertices(1),p.vertices(2)));
    }
}
//#####################################################################
// Function Cut_By_Plane
//#####################################################################
template<class T> void DISPLAY_PRIMITIVE_CUTTING<T>::
Cut_By_Plane(const PLANE<T> &p,T tol)
{
    ARRAY<VECTOR<TV, 3> > cut_elements;
    if(type==SEGMENT){
        for(int i=0;i<elements.m;i++){
            VECTOR<T,2> phi_nodes;
            for(int i=0;i<2;i++){phi_nodes(i)=p.Signed_Distance(vertices(i));}
            if(phi_nodes(0)*phi_nodes(1)<0){
                T theta=LEVELSET_UTILITIES<T>::Theta(phi_nodes(0),phi_nodes(1));
                TV interface_location=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(vertices(0),vertices(1),theta);
                if(theta>tol) cut_elements.Append(VECTOR<TV,3>(vertices(0),interface_location,interface_location));
                if(1-theta>tol) cut_elements.Append(VECTOR<TV,3>(interface_location,vertices(1),vertices(1)));}}}
    if(type==TRIANGLE){
        for(int i=0;i<elements.m;i++){
            VECTOR<T,3> phi_nodes;
            for(int i=0;i<3;i++){phi_nodes(i)=p.Signed_Distance(vertices(i));}
            int positive_count=0,single_node_sign;
            for(int i=0;i<3;i++) if(phi_nodes(i)>0) positive_count++;
            switch(positive_count){
                case 0:case 3:
                    cut_elements.Append(vertices);break;
                case 1:case 2:
                    single_node_sign=positive_count==1?1:-1;
                        for(int i=0;i<3;i++)if(LEVELSET_UTILITIES<T>::Sign(phi_nodes(i))==single_node_sign){
                                VECTOR<TV,2> interface_locations;int index=(i+1)%3;
                                VECTOR<T,2> theta;
                                VECTOR<int,2> other_locations;
                                for(int j=0;j<2;j++,index=(index+1)%3){
                                    other_locations(j)=index;
                                    theta(j)=LEVELSET_UTILITIES<T>::Theta(phi_nodes(i),phi_nodes(index));
                                    interface_locations(j)=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(vertices(i),vertices(index),theta(j));}
                                if(1-theta(0)>tol) cut_elements.Append(VECTOR<TV,3>(interface_locations(0),vertices(other_locations(0)),vertices(other_locations(1))));
                                if(1-theta(1)>tol) cut_elements.Append(VECTOR<TV,3>(interface_locations(0),vertices(other_locations(1)),interface_locations(1)));
                                if(theta(0)>tol && theta(1)>tol) cut_elements.Append(VECTOR<TV,3>(vertices(i),interface_locations(0),interface_locations(1)));}}}}
    elements.Exchange(cut_elements);    
}
//#####################################################################
// Function Get_Cutting_Plane
//#####################################################################
template<class T> PLANE<T> DEPTH_BUFFERING<T>::
Get_Cutting_Plane(const TV &a,const TV &b)
{ 
    return PLANE<T>(Embed(Project(a-b).Rotate_Clockwise_90().Normalized()),a);
}
//#####################################################################
// Function Add_Element
//#####################################################################
template<class T> int DEPTH_BUFFERING<T>::
Add_Element(const TV &a,const int style)
{ 
    primitives_cutting.Append(DISPLAY_PRIMITIVE_CUTTING<T>(a,style));
    return primitives_cutting.Size()-1;
}
//#####################################################################
// Function Add_Element
//#####################################################################
template<class T> int DEPTH_BUFFERING<T>::
Add_Element(const TV &a,const TV &b,const int style)
{
    primitives_cutting.Append(DISPLAY_PRIMITIVE_CUTTING<T>(a,b,style));
    return primitives_cutting.Size()-1;
}
//#####################################################################
// Function Add_Element
//#####################################################################
template<class T> int DEPTH_BUFFERING<T>::
Add_Element(const TV &a,const TV &b,const TV &c,const int style)
{
    primitives_cutting.Append(DISPLAY_PRIMITIVE_CUTTING<T>(a,b,c,style));
    return primitives_cutting.Size()-1;
}
//#####################################################################
// Function Process_Primitives
//#####################################################################
template<class T> void DEPTH_BUFFERING<T>::
Process_Primitives()
{
    for(int i=0;i<primitives_cutting.m;i++) {
        primitives_cutting(i).Initialize_Elements();
        primitives_cutting(i).Initialize_Bounding_Box();
    };
    for(int i=0;i<primitives_cutting.m;i++) for(int j=0;j<primitives_cutting.m;j++) if(i!=j) primitives_cutting(i).Cut_By_Primitive(primitives_cutting(j));
    
    
}
//#####################################################################
template class DISPLAY_PRIMITIVE<float>;
template class DEPTH_BUFFERING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DISPLAY_PRIMITIVE<double>;
template class DEPTH_BUFFERING<double>;
#endif
