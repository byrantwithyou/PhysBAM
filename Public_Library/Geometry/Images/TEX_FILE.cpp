#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Images/TEX_FILE.h>
#include <iomanip>
using namespace PhysBAM;
template<class T> TEX_FILE<T>::
//#####################################################################
// Constructor
//#####################################################################
TEX_FILE(const std::string& filename,const RANGE<TV>& box)
    :VECTOR_IMAGE<T>(filename,box),unit_offset(0),frame_offset(0)
{
    TV size=output_box.Edge_Lengths();
    stream<<"\\documentclass{article}\n";
    stream<<"\\usepackage[margin=0cm,papersize={"<<size.x<<"px,"<<size.y<<"px}]{geometry}\n";
    stream<<"\\usepackage{pstricks}\n";
    stream<<"\\usepackage{color}\n";
    stream<<"\\usepackage{pst-plot}\n";
    stream<<"\\definecolor{fc}{rgb}{0,0,0}\n";
    stream<<"\\definecolor{lc}{rgb}{0,0,0}\n";
    stream<<"\\begin{document}\n";
    stream<<"\\psset{unit=";
    unit_offset=stream.tellp();
    stream<<"                         px}\n";
    stream<<"\\noindent\\begin{pspicture}";
    frame_offset=stream.tellp();
    stream<<"                                                      "<<std::endl;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TEX_FILE<T>::
~TEX_FILE()
{
    if(!stream) return;
    stream<<"\\end{pspicture}\n";
    stream<<"\\end{document}\n";
    stream.seekp(frame_offset,std::ios::beg);
    Emit(bounding_box.min_corner);
    Emit(bounding_box.max_corner);
    stream.seekp(unit_offset,std::ios::beg);
    TV bb=bounding_box.Edge_Lengths(),ob=output_box.Edge_Lengths();
    T sc=1e5;
    if(bb.x>0) sc=min(sc,ob.x/bb.x);
    if(bb.y>0) sc=min(sc,ob.y/bb.y);
    stream<<sc;
}
//#####################################################################
// Function Emit
//#####################################################################
template<class T> void TEX_FILE<T>::
Emit(const TV &pt)
{
    stream.setf(std::ios_base::fixed);
    stream<<"("<<pt.x<<","<<pt.y<<")";
}
//#####################################################################
// Function Update_Effective_Formatting
//#####################################################################
template<class T> void TEX_FILE<T>::
Update_Effective_Formatting()
{
    if(effective_line_color!=cur_format.line_color){
        effective_line_color=cur_format.line_color;
        stream<<"\\definecolor{lc}{rgb}{"<<effective_line_color.x<<","<<effective_line_color.y<<","<<effective_line_color.z<<"}\n";}
    if(effective_fill_color!=cur_format.fill_color){
        effective_fill_color=cur_format.fill_color;
        stream<<"\\definecolor{fc}{rgb}{"<<effective_fill_color.x<<","<<effective_fill_color.y<<","<<effective_fill_color.z<<"}\n";}
}
//#####################################################################
// Function Emit_Options
//#####################################################################
template<class T> bool TEX_FILE<T>::
Emit_Options(bool line,bool fill)
{
    static const char* line_style_str[4]={"none","solid","dotted","dashed"};
    stream<<"[";
    bool split=line && fill && cur_format.line_style && cur_format.line_opacity && cur_format.fill_style && cur_format.fill_opacity && cur_format.line_opacity!=cur_format.fill_opacity;
    if(!split && line && cur_format.line_style && cur_format.line_opacity){
        stream<<"linewidth="<<cur_format.line_width;
        stream<<",linecolor=lc";
        stream<<",linejoin="<<cur_format.line_join;
        if(cur_format.line_style>1) stream<<",linestyle="<<line_style_str[cur_format.line_style];
        if(cur_format.line_opacity!=1) stream<<",opacity="<<cur_format.line_opacity;
        if(cur_format.arrow_style) stream<<",arrows="<<cur_format.arrow_style;}
    else stream<<"linestyle=none";

    if(fill && cur_format.fill_style && cur_format.fill_opacity){
        stream<<",fillcolor=fc";
        if(cur_format.fill_style==1) stream<<",fillstyle=solid";
        if(cur_format.fill_opacity!=1) stream<<",opacity="<<cur_format.fill_opacity;}

    if(cur_format.misc.length()>0) stream<<","<<cur_format.misc;
    stream<<"]";
    return split;
}
//#####################################################################
// Function Draw_Point
//#####################################################################
template<class T> void TEX_FILE<T>::
Emit_Object(const TV &pt,T radius)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();

    stream<<"\\pscircle";
    bool split=Emit_Options(true,true);
    Emit(pt);
    stream<<"{"<<radius<<"}"<<std::endl;

    if(split){
        stream<<"\\pscircle";
        Emit_Options(true,false);
        Emit(pt);
        stream<<"{"<<radius<<"}"<<std::endl;}
}
//#####################################################################
// Function Draw_Line
//#####################################################################
template<class T> void TEX_FILE<T>::
Emit_Object(const TV &a,const TV &b)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();

    stream<<"\\psline";
    bool split=Emit_Options(true,true);
    Emit(a);
    Emit(b);
    stream<<std::endl;

    if(split){
        stream<<"\\psline";
        Emit_Options(true,false);
        Emit(a);
        Emit(b);
        stream<<std::endl;}
}
//#####################################################################
// Function Draw_Line
//#####################################################################
template<class T> void TEX_FILE<T>::
Emit_Object(const TV &a,const TV &b,const TV &c)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();

    stream<<"\\pspolygon";
    bool split=Emit_Options(true,true);
    Emit(a);
    Emit(b);
    Emit(c);
    stream<<std::endl;

    if(split){
        stream<<"\\pspolygon";
        Emit_Options(true,false);
        Emit(a);
        Emit(b);
        Emit(c);
        stream<<std::endl;}
}
//#####################################################################
// Function Draw_Object
//#####################################################################
template<class T> void TEX_FILE<T>::
Emit_Object(const RANGE<TV>& box)
{
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();
    TV a(box.min_corner.x,box.max_corner.y),b(box.max_corner.x,box.min_corner.y);

    stream<<"\\pspolygon";
    bool split=Emit_Options(true,true);
    Emit(box.min_corner);
    Emit(a);
    Emit(box.max_corner);
    Emit(b);
    stream<<std::endl;

    if(split){
        stream<<"\\pspolygon";
        Emit_Options(true,false);
        Emit(box.min_corner);
        Emit(a);
        Emit(box.max_corner);
        Emit(b);
        stream<<std::endl;}
}
//#####################################################################
// Function Emit_Object
//#####################################################################
template<class T> void TEX_FILE<T>::
Emit_Object(ARRAY_VIEW<TV> pts)
{
    if(!pts.Size()) return;
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();

    stream<<"\\pspolygon";
    bool split=Emit_Options(true,true);
    for(int i=0;i<pts.Size();i++) Emit(pts(i));
    stream<<std::endl;

    if(split){
        stream<<"\\pspolygon";
        Emit_Options(true,false);
        for(int i=0;i<pts.Size();i++) Emit(pts(i));
        stream<<std::endl;}
}
//#####################################################################
// Function Emit_Object
//#####################################################################
template<class T> void TEX_FILE<T>::
Emit_Object(ARRAY_VIEW<TV> pts,ARRAY_VIEW<ARRAY_VIEW<TV> > holes)
{
    // TODO: fix
    if(!pts.Size()) return;
    if(!cur_format.line_style && !cur_format.fill_style) return;
    Update_Effective_Formatting();

    stream<<"\\pspolygon";
    bool split=Emit_Options(true,true);
    for(int i=0;i<pts.Size();i++) Emit(pts(i));
    stream<<std::endl;

    if(split){
        stream<<"\\pspolygon";
        Emit_Options(true,false);
        for(int i=0;i<pts.Size();i++) Emit(pts(i));
        stream<<std::endl;}
}
namespace PhysBAM{
template class TEX_FILE<float>;
template class TEX_FILE<double>;
}
