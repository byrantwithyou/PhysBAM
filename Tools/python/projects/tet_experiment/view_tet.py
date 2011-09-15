#!/usr/bin/python
import physbam
import sys
import OpenGL.GLUT as GLUT
import OpenGL.GL as GL
import OpenGL.GLU as GLU
import math

from  random_tet import RANDOM_TET

class WORLD:
    def __init__(self):
        self.alpha=0
        self.eye=physbam.Vf3(0,.5,-5)
        self.lookat=physbam.Vf3(.01,.01,.01)
        self.up=physbam.Vf3(0,1,0)
        self.width,self.height=640,480
        self.edge_show,self.point_show=None,None
        self.argv=GLUT.glutInit(sys.argv)
        GLUT.glutInitDisplayMode(GLUT.GLUT_DOUBLE|GLUT.GLUT_RGBA|GLUT.GLUT_DEPTH)
        #OpenGL.GLUT.glutDisplayMode(OpenGL.GLUT.GLUT_DOUBLE|OpenGL.GLUT.GLUT_RGBA|OpenGL.GLUT.GLUT_DEPTH)
        GLUT.glutInitWindowSize(self.width,self.height)

        self.main_window=GLUT.glutCreateWindow('tetview')
        GLUT.glutDisplayFunc(self.Display)


        self.rand=RANDOM_TET(True,physbam.Vf3(.5*math.cos(0),0,.5*math.sin(0)),
                             physbam.Vf3(.5*math.cos(-2*math.pi/3),0,.5*math.sin(-2*math.pi/3)),
                             physbam.Vf3(.5*math.cos(-4*math.pi/3),0,.5*math.sin(-4*math.pi/3)),
                             physbam.Vf3(0,.5,0))
        self.tet=self.rand.tet
        #self.tet=None
        #self.tet=physbam.TETRAHEDRON_f()
        #print self.tet.Signed_Volume()
        self.Center()

        self.initialized=False
        GLUT.glutTimerFunc(100,self.Timer,0)
        GLUT.glutReshapeFunc(self.Reshape)
        GLUT.glutKeyboardUpFunc(self.KeyboardUp)
        GLUT.glutMainLoop()


    def renderText(self,lines,highlight):
        count=0
        GLU.gluOrtho2D(0,self.width,0,self.height)
        font=GLUT.GLUT_BITMAP_9_BY_15
        height=20
        offset=8
        width=10
        for i in lines:
            if count==highlight:
                ymin=self.height-20-(count-1)*height-4
                ymax=self.height-20-(count)*height-4
                width=180 # GLUT.glutBitmapWidth(font,i)
                GL.glColor3f(1,0,0)
                GL.glRectf(offset,ymin,width+offset,ymax)
                GL.glColor3f(1,1,1)
            else:
                GL.glColor3f(0,0,0)
            GL.glRasterPos(offset,self.height-20-count*height)
            for j in i:
                GLUT.glutBitmapCharacter(font,ord(j))
            count+=1

    def conditionText(self,condition,text):
        if condition: return text
        else: return ' '*len(text)
    def formatStuff(self):
        lines=[]
        highlight=-1
        for i in range(1,5):
            rank=self.rand.ranked_lookup[i][0]
            lines.append("%d PF%2d %s %s %s"%(rank,i,self.conditionText(i==self.rand.pf_minimum_length,"min"),
                                              self.conditionText(self.rand.weightsValid(self.rand.pf_items[i][5]),"good"),
                                              self.conditionText(self.rand.gmin_distance_primitive==i,"gmin")))
            if i==self.point_show:
                highlight=i-1
        for i in range(1,4):
            rank=self.rand.ranked_lookup[i+4][0]
            lines.append("%d EE%2d %s %s %s"%(rank,i,self.conditionText(i==self.rand.ee_minimum_length,"min"),
                                              self.conditionText(self.rand.weightsValid(self.rand.ee_items[i][6]),"good"),
                                              self.conditionText(self.rand.gmin_distance_primitive==i+4,"gmin")))
            if i==self.edge_show:
                highlight=i-1+4
            
        lines.append("Tet volume: %f"%self.tet.Signed_Volume())

        return lines,highlight
        
        
    def KeyboardUp(self,key,x,y):
        if key=="n":
            self.rand=RANDOM_TET()
            self.tet=self.rand.tet
            self.Center()
        #points
        if key=="1": self.edge_show=None;self.point_show=1
        if key=="2": self.edge_show=None;self.point_show=2
        if key=="3": self.edge_show=None;self.point_show=3
        if key=="4": self.edge_show=None;self.point_show=4
        #edges
        if key=="5": self.edge_show=1;self.point_show=None
        if key=="6": self.edge_show=2;self.point_show=None
        if key=="7": self.edge_show=3;self.point_show=None
            

    def Center(self):
        self.lookat=self.tet.Center()

    def Timer(self,x):
        self.alpha+=.1
        self.eye=2*physbam.Vf3(math.cos(self.alpha),0,math.sin(self.alpha))+self.lookat
        GLUT.glutPostRedisplay();
        GLUT.glutTimerFunc(50,self.Timer,0)

    def Reshape(self,w,h):
        self.width,self.height=w,h
        GL.glViewport(0,0,self.width,self.height)
        GLUT.glutPostRedisplay();

    def Display(self):
        GL.glClearColor(1,1,1,0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT|GL.GL_DEPTH_BUFFER_BIT);

        if not self.initialized:
            self.initialized=True
            GL.glBlendFunc(GL.GL_SRC_ALPHA,GL.GL_ONE_MINUS_SRC_ALPHA)
            GL.glEnable(GL.GL_DEPTH_TEST)
            GL.glFrontFace(GL.GL_CCW);
            GL.glPolygonMode(GL.GL_FRONT_AND_BACK,GL.GL_FILL);

            self.quad=GLU.gluNewQuadric()

            GL.glLineWidth(3.0)
            GL.glDepthMask(GL.GL_TRUE)
            GL.glEnable(GL.GL_LINE_SMOOTH)
            GL.glDisable(GL.GL_CULL_FACE)
            GL.glEnable(GL.GL_BLEND)

        # setup text transform
        GL.glMatrixMode(GL.GL_PROJECTION)
        GL.glLoadIdentity()
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glLoadIdentity()
        lines,highlight=self.formatStuff()
        self.renderText(lines,highlight)

        GL.glMatrixMode(GL.GL_PROJECTION)
        GLU.gluPerspective(60,float(self.width)/self.height,.1,100.)


        # lighting
        #GL.glEnable(GL.GL_LIGHTING)
        #GL.glEnable(GL.GL_LIGHT0)
        #GL.glLightfv(GL.GL_LIGHT0,GL.GL_DIFFUSE,(.3,.3,.3,1))
        #GL.glLightfv(GL.GL_LIGHT0,GL.GL_SPECULAR,(.3,.3,.3,1))
        #GL.glLightfv(GL.GL_LIGHT0,GL.GL_POSITION,(0,.3,1,0))
            
        # setup camera
        GL.glMatrixMode(GL.GL_MODELVIEW)
        GL.glLoadIdentity()
        GLU.gluLookAt(self.eye.x,self.eye.y,self.eye.z,self.lookat.x,self.lookat.y,self.lookat.z,self.up.x,self.up.y,self.up.z)

        # material
        GL.glEnable(GL.GL_COLOR_MATERIAL)
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK,GL.GL_AMBIENT,(.1,.1,.1,1))
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK,GL.GL_DIFFUSE,(1,1,1,1))
        GL.glMaterialfv(GL.GL_FRONT_AND_BACK,GL.GL_SPECULAR,(1,1,1,1))

        # draw axes
        GL.glLineWidth(1.)
        GL.glBegin(GL.GL_LINES);
        GL.glColor3f(1,0,0);GL.glVertex3f(0,0,0);GL.glVertex3f(1,0,0)
        GL.glColor3f(0,1,0);GL.glVertex3f(0,0,0);GL.glVertex3f(0,1,0)
        GL.glColor3f(0,0,1);GL.glVertex3f(0,0,0);GL.glVertex3f(0,0,1)
        GL.glEnd()
        GL.glLineWidth(3.)

        # sphere point helper function
        def drawPoint(center):
            GL.glPushMatrix()
            GL.glTranslatef(*center)
            GLU.gluSphere(self.quad,.025,10,10)
            GL.glPopMatrix()

        # draw point face pair
        if self.point_show:
            tri,nodes,area,normal,distance,bary=self.rand.pf_items[self.point_show]
            GL.glBegin(GL.GL_TRIANGLES)
            GL.glColor3f(.5,.5,.5)
            GL.glVertex3f(*tri.x1)
            GL.glVertex3f(*tri.x2)
            GL.glVertex3f(*tri.x3)
            GL.glEnd()
            GL.glBegin(GL.GL_LINES)
            if self.rand.weightsValid(bary): GL.glColor3f(0,0,1)
            else: GL.glColor3f(1,0,0)
            face_point=tri.Point_From_Barycentric_Coordinates(bary)
            GL.glVertex3f(*face_point)
            GL.glVertex3f(*self.rand.Xs[nodes[0]])
            GL.glColor4f(1,0,1,.5)
            GL.glVertex3f(*tri.Point_From_Barycentric_Coordinates(bary))
            GL.glVertex3f(*tri.x1)
            GL.glVertex3f(*tri.Point_From_Barycentric_Coordinates(bary))
            GL.glVertex3f(*tri.x2)
            GL.glVertex3f(*tri.Point_From_Barycentric_Coordinates(bary))
            GL.glVertex3f(*tri.x3)
            GL.glEnd()

            GL.glColor3f(.5,0,.5)
            drawPoint(self.rand.Xs[nodes[0]])
            drawPoint(face_point)

        # draw edge edge pair
        if self.edge_show:
            edge1,edge2,join_edge,nodes,s,distance,weights=self.rand.ee_items[self.edge_show]
            GL.glLineWidth(6.0)

            GL.glColor3f(.5,0,.5)
            drawPoint(join_edge.x1)
            drawPoint(join_edge.x2)

            GL.glBegin(GL.GL_LINES)
            GL.glColor3f(1,0,0)
            GL.glVertex3f(*edge1.x1)
            GL.glVertex3f(*edge1.x2)
            GL.glVertex3f(*edge2.x1)
            GL.glVertex3f(*edge2.x2)
            if self.rand.weightsValid(weights): GL.glColor3f(0,0,1)
            else: GL.glColor3f(1,0,0)
            GL.glVertex3f(*join_edge.x1)
            GL.glVertex3f(*join_edge.x2)
            GL.glEnd()
            GL.glLineWidth(3.0)

        # draw the tet and its edges
        if self.tet != None:
            triangles=[self.tet.triangle1,self.tet.triangle2,self.tet.triangle3,self.tet.triangle4]
            #colors=[(1,.2,.2,.1),(.2,1,.2,.1),(.2,.2,1,.1),(1,.2,1,.1)]
            trans=.5
            colors=[(1,0,0,trans),(0,1,0,trans),(0,0,1,trans),(.5,.5,.5,trans)]
            for i in range(len(triangles)):
                tri=triangles[i]
                #TODO: putback
                #GL.glDepthMask(0)
                GL.glBegin(GL.GL_TRIANGLES) # GL_LINE_LOOP)
                #GL.glBegin(GL.GL_LINE_LOOP)
                GL.glColor4f(*colors[i])
                GL.glVertex3f(*tri.x1)
                GL.glVertex3f(*tri.x2)
                GL.glVertex3f(*tri.x3)
                GL.glEnd()
                GL.glDepthMask(1)
                
                GL.glBegin(GL.GL_LINES)
                GL.glColor4f(0,0,0,1)
                GL.glVertex3f(*tri.Center())
                GL.glVertex3f(*(tri.Center()+.1*tri.Normal()))
                GL.glEnd()

                GL.glBegin(GL.GL_LINE_LOOP)
                GL.glColor4f(0,0,0,1)
                GL.glVertex3f(*tri.x1)
                GL.glVertex3f(*tri.x2)
                GL.glVertex3f(*tri.x3)
                GL.glEnd()


        GLUT.glutSwapBuffers()

    
WORLD()
