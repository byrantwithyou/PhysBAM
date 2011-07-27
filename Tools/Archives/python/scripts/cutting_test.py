#!/usr/bin/python
from physbam import *
import cairo


def pt(x):  return tuple(map(lambda x: int(100*x),x))
def subset(X,nodes): return map(lambda x: X[x],nodes)

def draw_embed(ctx,area):
    for nodes in area.mesh.elements:
        X=subset(area.particles.X,nodes)
        ctx.set_source_rgba (0.3, 0.2, 0.5,.5)
        ctx.move_to(*pt(X[0]))
        ctx.line_to(*pt(X[1]))
        ctx.line_to(*pt(X[2]))
        ctx.close_path()
        ctx.fill()
        
def draw_cut(ctx,cut):
    ctx.set_source_rgb (0.0, 0.0, 0.0)
    for nodes in cut.mesh.elements:
        X=subset(cut.particles.X,nodes)
        ctx.new_path()
        ctx.move_to(*pt(X[0]))
        ctx.line_to(*pt(X[1]))
        ctx.stroke()
        
        

# Embedding
area=TRIANGULATED_AREA_f.Create()
area.particles.Add_Particles(3)
area.particles.X[1]=(.7,1)
area.particles.X[2]=(1.4,1.7)
area.particles.X[3]=(.3,1.8)
area.mesh.elements.Append((1,2,3))
area.Update_Number_Nodes()

# Cutting
cutting=CUTTING_GEOMETRY_2D_Vf2_2(True)
cutting.Initialize_Original_Embedding(area)

# Cutting Surface
cut=cutting.cutting
cut.particles.Add_Particles(2)
cut.particles.X[1]=(0,1.2)
cut.particles.X[2]=(2,1.2)
cut.mesh.elements.Append((1,2))
cut.Update_Number_Nodes()


# Do cut
area2=TRIANGULATED_AREA_f.Create()
cutting.Cut_Material(area2)
#cutting.Cut_Material(area2)

# Draw
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, 640,480)
ctx=cairo.Context(surface)
ctx.set_line_width(1)
draw_embed(ctx,area)
draw_cut(ctx,cut)
surface.write_to_png("triangle.png")
