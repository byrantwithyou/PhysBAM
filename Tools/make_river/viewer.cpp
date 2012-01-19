#include "BEZIER_SPLINE.h"
#include <Fl/Fl.h>
#include <Fl/Fl_Button.h>
#include <Fl/Fl_Double_Window.h>
#include <FL/Fl_Gl_Window.h>
#include <GL/gl.h>
#include <GL/glu.h>

using namespace PhysBAM;

template<class T>
void OpenGL_Vertex(VECTOR_2D<T>& X)
{glVertex2f(X.x,X.y);}

BEZIER_SPLINE<float,VECTOR_2D<float> > spline;
// make a Gl window that draws something:
template<class T>
class MyGlWindow : public Fl_Gl_Window {
public:
    int serial;
    GRID_2D<T> grid;
    int drag_point,drag_point_2;
    
    MyGlWindow(int x, int y, int w, int h) 
        :Fl_Gl_Window(x,y,w,h,"My GL Window")
    {
        grid.Initialize(w,h,-.3,1.3,-.3,1.3);
        serial=0;drag_point_2=drag_point=-1;
    }

private:
    void draw() 
    {glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(grid.xmin,grid.xmax,grid.ymin,grid.ymax);
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glColor3f(.2,.2,.2);
    glRectf(0,0,1,1);
    glColor3f(0,0,1);
    glLineWidth(1.0);
    glBegin(GL_LINE_STRIP);
    if(spline.control_points.m>2) for(int i=0;i<spline.segments;i++) for(int k=0;k<spline.samples_per_segment;k++){
            VECTOR_2D<T> position=spline.f(i,k*spline.one_over_samples_per_segment);glVertex2f(position.x,position.y);}
    glEnd();
    
    glColor3f(1,0,0);
    glPointSize(5.0);
    
    for(int s=0;s<spline.segments;s++){
        glBegin(GL_POINTS);
        int base_index=4*(s-1)+1;
        glColor3f(1,1,0);OpenGL_Vertex(spline.control_points(base_index));OpenGL_Vertex(spline.control_points(base_index+3));
        glColor3f(0,1,0);OpenGL_Vertex(spline.control_points(base_index+1));OpenGL_Vertex(spline.control_points(base_index+2));
        glEnd();
        glLineWidth(2.0);
        glBegin(GL_LINES);
        glColor3f(0,0,0.5);OpenGL_Vertex(spline.control_points(base_index));OpenGL_Vertex(spline.control_points(base_index+1));
        OpenGL_Vertex(spline.control_points(base_index+3));OpenGL_Vertex(spline.control_points(base_index+2));
        glEnd();}
    
    if(drag_point>0){
        glBegin(GL_LINES);glColor3f(0,0,0.5);
        OpenGL_Vertex(spline.control_points(drag_point-1));OpenGL_Vertex(spline.control_points(drag_point));
        glEnd();glBegin(GL_POINTS);glColor3f(1,0,0);OpenGL_Vertex(spline.control_points(drag_point-1));glEnd();}}

    int handle(int event)
    {VECTOR_2D<int> index(Fl::event_x()+1,h()-Fl::event_y());
    if(Fl::event_button()==1){
        if(event==FL_PUSH){
            if(spline.control_points.m>0){drag_point_2=spline.Add_Point(grid.X(index));spline.Add_Point(grid.X(index));}
            spline.Add_Point(grid.X(index));drag_point=spline.Add_Point(grid.X(index));}
        else if(event==FL_DRAG){
            if(drag_point>0) spline.control_points(drag_point)=grid.X(index);
            if(drag_point_2>0) spline.control_points(drag_point_2)=(T)2*spline.control_points(drag_point_2+1)-grid.X(index);}
        if(event==FL_RELEASE){drag_point=drag_point_2=0;}
        damage(1);
        return 1;}
    else if(Fl::event_button()==3){
        FILE* fp=fopen("cps.txt","w");
        for(int i=0;i<spline.control_points.m;i++){
            VECTOR_2D<T> cp=spline.control_points(i);
            fprintf(fp,"spline.Add_Point(VECTOR_2D<T>(%f,%f));\n",cp.x,cp.y);}
        fclose(fp);
        return 1;}
    return 0;}
};

main()
{
    MyGlWindow<float> wind(10,10,600,600);wind.show();Fl::run();
    return 0;
}
