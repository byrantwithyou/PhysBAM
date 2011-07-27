//#####################################################################
// Class MAIN_CLASS for Paramerizing the triangle mesh
// Modified by Zhaosheng Bao.
//##################################################################### 
//
//#####################################################################
//#####################################################################
#include "dk.h"
#include "gl/glut.h"
#include "LSCM.h"
#include "SphericalPM.h"
#include <time.h>
//#################################################################
// User parameters
//#################################################################
static double zoom=1.f;
double swingangle=0;
double elevateangle=0;
bool ShowBlue=false;
bool ShowRed=false;
int vertexpointer=1;
int fixedvertex[2];
DK dk;

bool SphereSolving=false;
SphericalPM spm;
//#####################################################################
// Class OPENGL_WRAPPER
//#####################################################################
class OPENGL_WRAPPER
{
public:
    OPENGL_WRAPPER(int *argc,char **argv)
    {
        glutInit(argc, argv);
        glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
        glutInitWindowPosition(100,100);
        glutInitWindowSize(800,600);
        glutCreateWindow ("DK Simplification Implementation");
        glutDisplayFunc(Handle_Display);
        glutReshapeFunc(Handle_Reshape);
        glutKeyboardFunc(Handle_Keypress);
        glutSpecialFunc(Handle_SpecialKeypress);
        glutMainLoop();
    }
    static void Handle_Display()
    {
        glClearColor(0,0,0,0);
        glClearDepth(1);
        glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
        glColor3f(1,1,1);
        glLoadIdentity();

        if( SphereSolving == true)
        {
            spm.Draw(swingangle, elevateangle);
            return;
        }


        gluLookAt(0, 0, 4, 0, 0, 0, 0, 1, 0);
        glRotated(swingangle, 1, 0, 0);
        glRotated(elevateangle, 0, 1, 0);


        for(int t=1; t<=dk.currentmesh->triangles.m; t++)
        {
            int ti,tj,tk;
            dk.currentmesh->triangles.Get(t,ti,tj,tk);
            VECTOR_3D<double> xi=dk.emb_boundary_particles.X(ti);
            VECTOR_3D<double> xj=dk.emb_boundary_particles.X(tj);
            VECTOR_3D<double> xk=dk.emb_boundary_particles.X(tk);
            glColor3f(1,1,1);
            glBegin(GL_LINE_LOOP);
    //        glNormal3dv(normal);
            glVertex3f(zoom*(xi.x-dk.aveX), zoom*(xi.y-dk.aveY), zoom*(xi.z-dk.aveZ));
            glVertex3f(zoom*(xj.x-dk.aveX), zoom*(xj.y-dk.aveY), zoom*(xj.z-dk.aveZ));
            glVertex3f(zoom*(xk.x-dk.aveX), zoom*(xk.y-dk.aveY), zoom*(xk.z-dk.aveZ));
            glEnd();

            if(ShowBlue==true || ShowRed==true)
            {
                for(int i=1; i<=dk.vertices_info(t).m; i++)
                {
                    int index=dk.vertices_info(t)(i);
                    VECTOR_3D<double> pi=dk.emb_boundary_particles.X(dk.vertices_info(t)(i));
                    VECTOR_3D<double> ppi=xi*(*dk.vlist)(index)->u + xj*(*dk.vlist)(index)->v
                            +xk*(1-(*dk.vlist)(index)->u-(*dk.vlist)(index)->v);
                    if(ShowBlue)
                    {
                        glColor3f(0, 0, 1);
                        glBegin(GL_POINTS);
                        glVertex3f(zoom*(pi.x-dk.aveX), zoom*(pi.y-dk.aveY), zoom*(pi.z-dk.aveZ));
                        glEnd();
                    }
                    if(ShowRed)
                    {
                        glColor3f(1, 0, 0);
                        glBegin(GL_POINTS);
                        glVertex3f(zoom*(ppi.x-dk.aveX), zoom*(ppi.y-dk.aveY), zoom*(ppi.z-dk.aveZ));
                        glEnd();
                    }
                }
            }
        }
         ARRAY <DK_VERTEX <double> *> *vertices=dk.vlist->vertex_sequence;

        for(int i=1; i<=vertices->m; i++)
        {
            if( vertexpointer==i)
            {
                VECTOR_3D<double> p=dk.emb_boundary_particles.X((*vertices)(i)->index);
                glColor3f(1, 0, 0);
                glPushMatrix();
                glTranslatef(zoom*(p.x-dk.aveX), zoom*(p.y-dk.aveY), zoom*(p.z-dk.aveZ));
                glutWireSphere(zoom*0.005, 20, 20);
                glPopMatrix();
            }
            //    else if( (*vertices)(i)->tag==true)
            else if( i==fixedvertex[0] || i==fixedvertex[1] )
            {
                VECTOR_3D<double> p=dk.emb_boundary_particles.X((*vertices)(i)->index);
                glColor3f(0, 1, 0);
                glPushMatrix();
                glTranslatef(zoom*(p.x-dk.aveX), zoom*(p.y-dk.aveY), zoom*(p.z-dk.aveZ));
                glutWireSphere(zoom*0.005, 20, 20);
                glPopMatrix();
            }
        }
        
        glutSwapBuffers();
    }

    static void Handle_Reshape(int w,int h)
    {
        glViewport(0,0,800, 600);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective (10, 4.0/3, 1, 10000);
        glMatrixMode(GL_MODELVIEW);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_DEPTH_TEST);
        glutPostRedisplay();
    }
    static void Handle_SpecialKeypress(int key, int x, int y)
    {
        if( key==101)    swingangle-=2;
        if( key==103)    swingangle+=2;
        if( key==100)    elevateangle+=2;
        if( key==102)    elevateangle-=2;
        glutPostRedisplay();
    }
    static void Handle_Keypress(unsigned char key,int mousex,int mousey)
    {
        switch(key)
        {
            case 'q':
                //exit(0);
            case 'Z':
                zoom+=0.05f;
                break;
           case 'z':
               zoom-=0.05f;
                break;
           case 'b':case'B':
               ShowBlue=!ShowBlue;
               break;
           case 'r':case'R':
               ShowRed=!ShowRed;
               break;
           case '.':case'>':
                {
                    vertexpointer++;
                    ARRAY <DK_VERTEX <double> *> *vertices=dk.vlist->vertex_sequence;
                    vertexpointer=(vertexpointer-1+vertices->m)%vertices->m+1;

                    int index=(*vertices)(vertexpointer)->index;
                    printf("Moved  vertex %d, ", vertexpointer);
                    if( (*dk.vlist)(index)->real_u<0 || (*dk.vlist)(index)->real_v<0)
                        printf("U, V: undefined.\n");
                    else printf("U: %f, V: %f.\n", (*dk.vlist)(index)->real_u, 
                        (*dk.vlist)(index)->real_v);
                }
               break;
           case ',':case'<':
               {
                    vertexpointer--;
                    ARRAY <DK_VERTEX <double> *> *vertices=dk.vlist->vertex_sequence;
                    vertexpointer=(vertexpointer-1+vertices->m)%vertices->m+1;
                    int index=(*vertices)(vertexpointer)->index;
                    printf("Moved to the vertex %d, ", vertexpointer);
                    if( (*dk.vlist)(index)->real_u<0 || (*dk.vlist)(index)->real_v<0)
                        printf("U, V: undefined.\n");
                    else printf("U: %f, V: %f.\n", (*dk.vlist)(index)->real_u, 
                        (*dk.vlist)(index)->real_v);
               }
               break;
           case 's':case'S':
            {
            /*    float u, v;
                printf("input U:");
                scanf("%f", &u);
                printf("input V:");
                scanf("%f", &v);

                ARRAY <DK_VERTEX <double> *> *vertices=dk.vlist->vertex_sequence;
                int index=(*vertices)(vertexpointer)->index;
                (*dk.vlist)(index)->tag=true;
                (*dk.vlist)(index)->real_u=u;
                (*dk.vlist)(index)->real_v=v;*/
                if( fixedvertex[0]==0)    fixedvertex[0]=vertexpointer;
                else if( fixedvertex[1]==0)    fixedvertex[1]=vertexpointer;
                else
                {
                    fixedvertex[0]=fixedvertex[1];
                    fixedvertex[1]=vertexpointer;
                }
            }
            break;
            case 'c':case'C':
                //LSCM(dk);
                break;
            case 'n':case 'N':
                if( SphereSolving==true)
                {
                    //Here you can set the iterations in each step.
                    //More steps needed (>2000 or more) if you use two fixed vertices.
                    spm.Interations(2000);
                    break;
                }
                else
                {
                    dk.Run();
                    dk.Init_Run();
                    printf("Triangles left after simplification: %d\n", dk.currentmesh->triangles.m);
                    break;
                }
            case 'o':
                if( dk.Output("param.txt")==true)
                   printf("Output file saved successfully.\n");
                break;
            case 'm':case 'M':
                SphereSolving=!SphereSolving;
                if( SphereSolving == true)
                {
                    swingangle=0;
                    elevateangle=0;
                    spm.Initialze(dk, fixedvertex[0], fixedvertex[1]);
                }
               break;
           case 'p':case 'P':
               if( SphereSolving == true )
               {
                   spm.Output(dk);
                   if( dk.Output("param.txt")==true)
                          printf("Output file saved successfully.\n");
               }
               break;
    
           default:
            break; 
        }
        glutPostRedisplay();
    }
};

//#######################################################################
// Function main
//#######################################################################
int main(int argc, char **argv)
{
    //You can use one fixed vertex, just set fixedvertex[1] to be 0.
    fixedvertex[0]=15732;
    fixedvertex[1]=21;
    char filename[100];
    sprintf(filename, "../../Public_Data/VH_Skin/arm_right.tri");
    dk.Open(filename);
    dk.currentmesh->Initialize_Incident_Triangles();
    dk.currentmesh->Initialize_Neighbor_Nodes();
    for(int i=1; i<0; i++)
    {
        dk.Run();
        printf("Triangles left after simplification %d: %d\n", i, dk.currentmesh->triangles.m);
    }

    dk.Init_Run();
    OPENGL_WRAPPER(&argc, argv);

    return 0;
}
