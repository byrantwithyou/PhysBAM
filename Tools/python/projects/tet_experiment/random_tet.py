#!/usr/bin/python

import physbam
import random
import math

rand=random.Random()

def triangle_nodes(tri_index): # (point,tri1,tri2,tri3)
    return {1:(1,2,4,3),2:(2,1,3,4),3:(3,1,4,2),4:(4,1,2,3)}[tri_index]

def edge_nodes(edge_index):
    return {1:(1,2,3,4),2:(2,3,1,4),3:(1,3,4,2)}[edge_index]


class RANDOM_TET:
    """Makes a random tetrahedron and computes the various lengths between skew edges and point-face altitudes"""
    def cyclic_shift(nodes):
        return [nodes[-1]]+nodes[:-1]

    def __init__(self,verbose=True):
        self.verbose=verbose
        self.compute()

    def __init__(self,verbose=True,x1=None,x2=None,x3=None,x4=None):
        self.verbose=verbose
        self.Xs=[0,x1,x2,x3,x4]
        if x1==None:
            self.Xs=[0,physbam.Vf3(rand.uniform(0,1),rand.uniform(0,1),rand.uniform(0,1)),
                     physbam.Vf3(rand.uniform(0,1),rand.uniform(0,1),rand.uniform(0,1)),
                     physbam.Vf3(rand.uniform(0,1),rand.uniform(0,1),rand.uniform(0,1)),
                     physbam.Vf3(rand.uniform(0,1),rand.uniform(0,1),rand.uniform(0,1))]
        self.compute()

    def weightsValid(self,w):
        for i in w:
            if i<0. or i>1.: return False
        return True
        
    def compute(self):
        self.nodes=[1,2,3,4]
        self.tet=physbam.TETRAHEDRON_f(*map(lambda x: self.Xs[x],self.nodes))


        Xs=self.Xs

        if self.verbose:
            print "Tet Volume %f"%self.tet.Signed_Volume()

        # go through triangles
        self.ranked_items=[] # (item#, distance, weights_good)

        self.pf_items=[None] # (nodes,area,normal,distance,bary)
        min_distance,min_distance_tri=1e99,0
        max_area,max_area_tri=0,0
        for tri_index in range(1,5):
            node1,node2,node3,node4=triangle_nodes(tri_index)
            tri=physbam.TRIANGLE_3D_f(Xs[node2],Xs[node3],Xs[node4])
            
            #cross_area=physbam.Vf3.Cross_Product(Xs[node3]-Xs[node2],Xs[node4]-Xs[node2]).Magnitude_Squared()
            area=tri.Area()
            normal=tri.Normal()
            bary=tri.Barycentric_Coordinates(Xs[node1])
            distance=math.fabs(physbam.Vf3.Dot_Product(Xs[node1]-Xs[node2],normal))
            #bary="Fu"
            #print "tri_index=%d area=%f bary=%s distance=%f"%(tri_index,tri.Area(),bary,distance)
            if distance<min_distance: min_distance,min_distance_tri=distance,tri_index
            if area>max_area: max_area,max_area_tri=area,tri_index
            self.ranked_items.append((tri_index,distance,self.weightsValid(bary)))
            self.pf_items.append((tri,(node1,node2,node3,node4),area,normal,distance,bary))
        self.pf_minimum_length=min_distance_tri
        if self.verbose:
            print min_distance_tri,max_area_tri
            print "---POINT FACE-------------------------------"
            print "Idx  Min Dist Max Area      Bary"
        for tri_index in range(1,5):
            tri,nodes,area,normal,distance,bary=self.pf_items[tri_index]
            
            max_area_string=" "
            if max_area_tri==tri_index: max_area_string="*"
            min_distance_string=" "
            if min_distance_tri==tri_index: min_distance_string="*"
            bary_ok="*"
    
            for i in bary:
                #print i
                if i<0. or i>1.: bary_ok=" "
    
            if self.verbose: print "%3d %s%10.5f %s%10.5f %s %-30s"%(tri_index,min_distance_string,distance,max_area_string,area,bary_ok,bary)
    
        if self.verbose: print "---EDGE EDGE-------------------------------"
        self.ee_items=[None] # (nodes,direction,distance,weights)
        min_distance_edge,min_distance_edge_index=1e99,0
        gmin_distance,self.gmin_distance_primitive=min_distance,min_distance_tri
        for edge_index in range(1,4):
            weights=physbam.Vf2()
            node1,node2,node3,node4=edge_nodes(edge_index)
            edge1=physbam.SEGMENT_3D_f(Xs[node1],Xs[node2])
            edge2=physbam.SEGMENT_3D_f(Xs[node3],Xs[node4])
            s=edge1.Shortest_Vector_Between_Lines(edge2,weights)
            join_segment=physbam.SEGMENT_3D_f(Xs[node1]*(1-weights.x)+Xs[node2]*weights.x,Xs[node3]*(1-weights.y)+Xs[node4]*weights.y)
            distance=s.Magnitude()
            if distance<min_distance_edge:
                min_distance_edge,min_distance_edge_index=distance,edge_index
            if distance<gmin_distance:
                gmin_distance,self.gmin_distance_primitive=distance,edge_index+4
            self.ee_items.append((edge1,edge2,join_segment,(node1,node2,node3,node4),s,distance,weights))
            self.ranked_items.append((edge_index+4,distance,self.weightsValid(weights)))
        self.ee_minimum_length=min_distance_edge_index
        for edge_index in range(1,4):
            edge1,edge2,join_segment,nodes,s,distance,weights=self.ee_items[edge_index]
            min_distance_string=" "
            if edge_index==min_distance_edge_index: min_distance_string="*"
            weights_ok="*"
            for i in weights:
                if i<0. or i>1.: weights_ok=" "
            if self.verbose: print "%3d %s%10.5f %s %-30s"%(tri_index,min_distance_string,distance,weights_ok,weights)
        if self.verbose: print "gmin %d"%self.gmin_distance_primitive

        def rank(x,y):
            good=cmp(y[2],x[2])
            if good==0:
                return cmp(x[1],y[1])
            else:
                return good
        self.ranked_items.sort(rank)
        self.ranked_lookup={}
        ranking=0
        for idx,good,dist in self.ranked_items:
            self.ranked_lookup[idx]=(ranking,good,dist)
            ranking+=1
            
    def computeArea(self):
        Xs=self.Xs

        print "Volume %f"%self.tet.Signed_Volume()
        
        for tri_index in range(1,5):
            node1,node2,node3,node4=triangle_nodes(tri_index)
            tri=physbam.TRIANGLE_3D_f(Xs[node2],Xs[node3],Xs[node4])
            u_cross_v=physbam.Vf3.Cross_Product(Xs[node3]-Xs[node2],Xs[node4]-Xs[node2])
            area=u_cross_v.Normalize()
            distance=physbam.Vf3.Dot_Product(u_cross_v,Xs[node1]-Xs[node3])
            print "Tri  %d %5.3f %5.3f %f"%(tri_index,area,distance,area*distance/6)

        for edge_index in range(1,4):
            node1,node2,node3,node4=edge_nodes(edge_index)
            u_cross_v=physbam.Vf3.Cross_Product(Xs[node2]-Xs[node1],Xs[node4]-Xs[node3])
            area=u_cross_v.Normalize()
            distance=physbam.Vf3.Dot_Product(u_cross_v,Xs[node1]-Xs[node3])
            print "Edge %d %5.3f %5.3f %f"%(edge_index,area,distance,area*distance/6)
        
            
            
        
            
if __name__=="__main__":
    #for i in range(1000000):
    #    if i%10000==0: print "%d"%i
    #    tet=RANDOM_TET(False)
    #    if tet.ranked_items[0][2]==False:
    #        print "ERROR %s"%tet.Xs
    tet=RANDOM_TET(False)
    tet.computeArea()
        
    
