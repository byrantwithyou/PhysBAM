#!/usr/bin/python

import os
import sys
import optparse
from struct import *

from physbam import *

parser = optparse.OptionParser("usage: %prog [options] input.tri output.tri")
parser.add_option('-d','--delete_connected',default=-1,type='int',help='delete connected components')
(options,args)=parser.parse_args()
if len(args)!=2: parser.error("invalid number of arguments")

tri=TRIANGULATED_SURFACE_f.Create()
Read_From_File("float",args[0],tri)

if 0:
    status=-1;
    if(options.delete_connected>=0):
        tri.mesh.Initialize_Neighbor_Nodes()
        components=0
        component=[0 for i in range(len(tri.particles)+1)]
        component_list=[[]]
        #loop through elems
        for i in range(1,len(tri.particles)+1):
            #print "Particle %i"%i
            new_status=(int(float(i)/float(len(tri.particles))*100))*1
            if (new_status is not status): print "%i"%new_status
            status=new_status
            #only deal with if we haven't already
            if (not component[i]):
                current_component=0;
                #loop through neighbors
                for j in range(1,len(tri.mesh.neighbor_nodes[i])):
                    n=tri.mesh.neighbor_nodes[i][j]
                    #if we need to merge components do so
                    if (component[n] and current_component and (component[n] is not current_component)):
                        old_component=component[n]
                        #print "Old %i"%len(component_list[old_component])
                        #print "Current %i"%len(component_list[current_component])
                        if(len(component_list[old_component])>len(component_list[current_component])):
                            tmp=current_component
                            current_component=old_component
                            old_component=tmp
                        for k in range(len(component_list[old_component])):
                            component[component_list[old_component][k]]=current_component
                        component_list[current_component].extend(component_list[old_component]);
                        component_list[old_component]=[]
                    #set the component
                    elif (component[n]):
                        current_component=component[n]
                        component[i]=current_component
                        component_list[component[n]].append(i);
                #if no neighbors have been thouched yet make a new component
                if (not current_component):
                    components+=1
                    component[i]=components
                    component_list.append([i]);

        for i in range(len(component_list)-1,-1,-1):
            if len(component_list[i]) is 0:
                component_list.pop(i);

        print component_list
        print len(component_list);

tri.mesh.Initialize_Neighbor_Nodes()
uf=UNION_FIND_int();
uf.Initialize(len(tri.particles));
total=0
connected_components={}
for i in range(1,len(tri.particles)+1):
    for j in range(1,len(tri.mesh.neighbor_nodes[i])):
        uf.Union(i,tri.mesh.neighbor_nodes[i][j])
for i in range(1,len(tri.particles)+1):
    total+=int(uf.Is_Root(i));
    root=uf.Find(i)
    if (root in connected_components): connected_components[root].append(i)
    else: connected_components[root]=[i]
#print connected_components
print total

surfaces={}
for i in range(1,len(tri.mesh.elements)+1):
    root1=uf.Find(tri.mesh.elements[i][1])
    root2=uf.Find(tri.mesh.elements[i][2])
    root3=uf.Find(tri.mesh.elements[i][3])
    if(root1==root2 and root2==root3):
        if not (root1 in surfaces):
            surfaces[root1]=TRIANGULATED_SURFACE_f.Create(tri.particles)
            surfaces[root1].Update_Number_Nodes()
        surfaces[root1].mesh.elements.Append(tri.mesh.elements[i])
for i,j in surfaces.iteritems():
    Write_To_File("float","%s_%i.tri"%(args[1],i),j)
    print "Wrote %s_%i"%(args[1],i)
