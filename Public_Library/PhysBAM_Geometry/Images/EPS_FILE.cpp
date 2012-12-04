#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Images/EPS_FILE.h>
using namespace PhysBAM;
template<class T> EPS_FILE<T>::
//#####################################################################
// Constructor
//#####################################################################
EPS_FILE(const std::string& filename,const RANGE<TV>& box)
    :VECTOR_IMAGE<T>(filename,box),head_offset(0),effective_line_width(-1),effective_point_radius(-1),effective_line_opacity(-1),effective_fill_opacity(-1)
{
    stream<<"%!PS-Adobe-3.0 EPSF-3.0"<<std::endl;
    stream<<"%%BoundingBox: ";
    Emit(output_box.min_corner);
    Emit(output_box.max_corner);
    stream<<std::endl;
    head_offset=stream.tellp();
    stream<<"                                                                                                    "<<std::endl;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> EPS_FILE<T>::
~EPS_FILE()
{
    if(!stream) return;
    T scale;
    TV shift;
    Compute_Transform(scale,shift);
    stream.seekp(head_offset,std::ios::beg);
    Emit(shift);
    stream<<"translate"<<std::endl;
    stream<<scale<<" "<<scale<<" scale"<<std::endl;
    stream<<3/scale<<" setlinewidth"<<std::endl;
}
//#####################################################################
// Function Emit
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit(const std::string& str)
{
    stream<<str<<" ";
}
//#####################################################################
// Function Emit
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit(const TV &pt)
{
    stream<<pt.x<<" "<<pt.y<<" ";
}
//#####################################################################
// Function Compute_Transform
//#####################################################################
template<class T> void EPS_FILE<T>::
Compute_Transform(T& scale,TV& shift)
{
    if(output_box.Empty() || bounding_box.Empty()){
        scale=1;
        shift=TV();}
    else{
        scale=(output_box.Edge_Lengths()/bounding_box.Edge_Lengths()).Min();
        shift=output_box.min_corner-bounding_box.min_corner*scale;}
}
//#####################################################################
// Function Update_Effective_Formatting
//#####################################################################
template<class T> void EPS_FILE<T>::
Stroke()
{
    if(!cur_format.line_style) return;
    if(effective_color!=cur_format.line_color){
        effective_color=cur_format.line_color;
        stream<<effective_color.x<<" "<<effective_color.y<<" "<<effective_color.z<<" setrgbcolor ";}
    Emit("stroke");
}
//#####################################################################
// Function Update_Effective_Formatting
//#####################################################################
template<class T> void EPS_FILE<T>::
Fill()
{
    if(!cur_format.fill_style) return;
    if(cur_format.line_style) Emit("gsave");
    if(effective_color!=cur_format.fill_color){
        effective_color=cur_format.fill_color;
        stream<<effective_color.x<<" "<<effective_color.y<<" "<<effective_color.z<<" setrgbcolor ";}
    Emit("fill");
    if(cur_format.line_style) Emit("grestore");
}
//#####################################################################
// Function Update_Effective_Formatting
//#####################################################################
template<class T> void EPS_FILE<T>::
Update_Effective_Formatting()
{
    // if(effective_color!=cur_format.line_color){
    //     effective_color=cur_format.line_color;
    //     stream<<effective_color.x<<" "<<effective_color.y<<" "<<effective_color.z<<" setrgbcolor"<<std::endl;}
//    T effective_line_width,effective_point_radius,effective_line_opacity,effective_fill_opacity;
}
//#####################################################################
// Function Draw_Point
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit_Object(const TV &pt,T radius)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();
    Emit("newpath");
    Emit(pt);
    stream<<radius<<" 0 360 arc closepath ";
    Fill();
    Stroke();
    stream<<std::endl;
}
//#####################################################################
// Function Draw_Line
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit_Object(const TV &a,const TV &b)
{
    if(!cur_format.line_style) return;
    Update_Effective_Formatting();
    Mt(a);
    Lt(b);
    Stroke();
    stream<<std::endl;
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit_Object(const TV &a,const TV &b,const TV &c)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();
    Emit("newpath");
    Mt(a);
    Lt(b);
    Lt(c);
    Emit("closepath");
    Fill();
    Stroke();
    stream<<std::endl;
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit_Object(const RANGE<TV>& box)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();
    TV a(box.min_corner.x,box.max_corner.y),b(box.max_corner.x,box.min_corner.y);
    Emit("newpath");
    Mt(box.min_corner);
    Lt(a);
    Lt(box.max_corner);
    Lt(b);
    Emit("closepath");
    Fill();
    Stroke();
    stream<<std::endl;
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit_Object(ARRAY_VIEW<TV> pts)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();
    if(!pts.Size()) return;
    Emit("newpath");
    Mt(pts(0));
    for(int i=1;i<pts.Size();i++) Lt(pts(i));
    Emit("closepath");
    Fill();
    Stroke();
    stream<<std::endl;
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void EPS_FILE<T>::
Emit_Object(ARRAY_VIEW<TV> pts,ARRAY_VIEW<ARRAY_VIEW<TV> > holes)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();
    if(!pts.Size()) return;
    Emit("newpath");
    Mt(pts(0));
    for(int i=1;i<pts.Size();i++) Lt(pts(i));
    Emit("closepath");
    for(int h=0;h<holes.Size();h++){
        Mt(holes(h)(0));
        for(int i=1;i<holes(h).Size();i++) Lt(holes(h)(i));
        Emit("closepath");}
    Fill();
    Stroke();
    stream<<std::endl;
}
//#####################################################################
// Function Mt
//#####################################################################
template<class T> void EPS_FILE<T>::
Mt(const TV &pt)
{
    Emit(pt);
    Emit("moveto");
}
//#####################################################################
// Function Lt
//#####################################################################
template<class T> void EPS_FILE<T>::
Lt(const TV &pt)
{
    Emit(pt);
    Emit("lineto");
}
namespace PhysBAM{
template class EPS_FILE<float>;
template class EPS_FILE<double>;
}
