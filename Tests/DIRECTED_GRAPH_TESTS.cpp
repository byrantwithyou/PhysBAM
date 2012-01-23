//#####################################################################
// Copyright 2012, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIRECTED_GRAPH_TESTS
//#####################################################################

#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Utilities/TEST_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

class DIRECTED_GRAPH_TESTS:public TEST_BASE
{
    bool Test(const bool test,const char* str_test,int line,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<str_test<<"    line "<<line<<std::endl;ok=false;return false;}

public:
    DIRECTED_GRAPH_TESTS()
        :TEST_BASE("directed_graph")
    {}

    virtual ~DIRECTED_GRAPH_TESTS(){}

#define TEST(x) Test((x),#x,__LINE__,ok)

    TEST_RESULT Run_Test(int n)
    {
        bool ok=true;
        DIRECTED_GRAPH<int> g(6);
        TEST(g.Number_Of_Levels()==0);
        TEST(g.Number_Nodes()==6);
        g.Initialize(3);
        TEST(g.Number_Nodes()==3);
        g.Resize(4);
        TEST(g.Number_Nodes()==4);
        TEST(g.Number_Of_Levels()==0);
        g.Initialize(6);
        TEST(g.Number_Nodes()==6);

        g.Add_Edge(0,1);
        g.Add_Edge(1,2,true);
        g.Add_Edge(4,1);
        g.Add_Edge(4,2);
        g.Add_Edge(4,3);
        g.Add_Edge(5,4);

        ARRAY<int> parents=g.Parents(4);
        TEST(parents.Size()==1);
        TEST(parents(0)==5);
        ARRAY<int> children=g.Children(4);
        TEST(children.Size()==3);

        g.Generate_Levels();
        TEST(g.Number_Of_Levels()==6);
        g.Add_Edge(2,5);
        g.Generate_Levels();
        TEST(g.Number_Of_Levels()!=6);
        g.Remove_Edge(2,5);
        g.Generate_Levels();
        TEST(g.Number_Of_Levels()==6);

        ARRAY<int,int> finish_time,node_index;
        TEST(g.Topological_Sort(finish_time,node_index)==true);
        TEST(node_index(0)==2);
        TEST(node_index(1)==1);
        TEST(node_index(2)==0);
        TEST(node_index(3)==3);
        TEST(node_index(4)==4);
        TEST(node_index(5)==5);
       
        g.Topological_Sort_Assuming_Cycle_Free(finish_time,node_index);
        TEST(node_index(0)==2);
        TEST(node_index(1)==1);
        TEST(node_index(2)==0);
        TEST(node_index(3)==3);
        TEST(node_index(4)==4);
        TEST(node_index(5)==5);

        ARRAY<int,int> depths;
        g.Maximal_Depth_On_Acyclic_Graph(depths);
        TEST(depths(0)==1);
        TEST(depths(1)==3);
        TEST(depths(2)==4);
        TEST(depths(3)==3);
        TEST(depths(4)==2);
        TEST(depths(5)==1);
       
        ARRAY<int,int> component_id;
        g.Add_Edge(2,5);
        g.Remove_Edge(0,1);
        TEST(g.Topological_Sort(finish_time,node_index)==false);

        DIRECTED_GRAPH<int> g2(2);
        g2.Initialize(2);
        g2.Add_Edge(0,1);
        TEST(g2.Topological_Sort(finish_time,node_index)==true);
        TEST(node_index(0)==1);
        TEST(node_index(1)==0);
        TEST(g2.Strongly_Connected_Components(component_id)==2);
        g2.Remove_Edge(0,1);
        TEST(g2.Strongly_Connected_Components(component_id)==2);
        g2.Add_Edge(0,1);
        g2.Add_Edge(1,0);
        TEST(g2.Strongly_Connected_Components(component_id)==1);
                
        return ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static DIRECTED_GRAPH_TESTS directed_graph_tests;
}
