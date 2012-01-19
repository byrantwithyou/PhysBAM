//#####################################################################
// Copyright 2009, Michael Lentine, Bridget Vuong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_DIRECTED_GRAPH.h>
#include <PhysBAM_Dynamics/Solids_Evolution/SEARCH_CONTROLLER.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>

using namespace PhysBAM;

template<class TV>
class PREPROCESS_CONTROLLER_DATA
    {
public:
    int number_nodes;
    typedef typename TV::SCALAR T;
    typedef GRID<TV> T_GRID;
    ARRAY<ENVIRONMENTAL_STATE<T_GRID>*> states;
    DIRECTED_GRAPH<int> states_graph;
    
    ARRAY<ENVIRONMENTAL_STATE<T_GRID>*> graph_index_to_state;
    int next_hash_number;
    
 
    PREPROCESS_CONTROLLER_DATA(std::string graph_file, std::string hashtable_file):states_graph(1)
    { 
        if(FILE_UTILITIES::File_Exists(graph_file)) FILE_UTILITIES::Read_From_File<float>(graph_file,states_graph);
        if(FILE_UTILITIES::File_Exists(hashtable_file)) FILE_UTILITIES::Read_From_File<float>(hashtable_file,graph_index_to_state);
        number_nodes=states_graph.Number_Nodes(); next_hash_number=graph_index_to_state.Size()+1;
    }

    ~PREPROCESS_CONTROLLER_DATA()
    {}

    void Process_File(std::string file)
    {
        FILE_UTILITIES::Read_From_File<float>(file,states);
        //assert(states(1)->angles.Size()!=0); //
        for(int i=0;i<states.Size();i++) {for(int j=i+1; j<=states.Size(); j++) {
            int current_index, next_index;
            if(!(current_index=graph_index_to_state.Find(states(i)))){
                graph_index_to_state.Resize(next_hash_number);
                current_index=next_hash_number;graph_index_to_state(next_hash_number++)=states(i);number_nodes++;}
            if(!(next_index=graph_index_to_state.Find(states(j)))){
                graph_index_to_state.Resize(next_hash_number);
                next_index=next_hash_number;graph_index_to_state(next_hash_number++)=states(j);number_nodes++;}
            if(number_nodes>=states_graph.Number_Nodes()) states_graph.Resize(states_graph.Number_Nodes()*2);
            states_graph.Add_Edge(current_index,next_index);}}
        //for(int i=0;i<graph_index_to_state.Size();i++) assert(graph_index_to_state(i)->angles.Size()!=0); //
    }
};
//#####################################################################
// Function main
//#####################################################################
int main(int argc,char *argv[])
{
    if(argc<4) {
        LOG::cout << std::endl << "Usage: ./preprocess_controller_data graph_output_file graph_index_to_state_output_file input_directory" << std::endl;
        exit(1);
    }
    PREPROCESS_CONTROLLER_DATA<VECTOR<float,3> > data(argv[1], argv[2]);
    //for(int i=0;i<=100;i++) data.Process_File(str(boost::format("%s/search_controller.run_%d")%argv[3]%i));
    data.Process_File("/disk2/bvuong/PhysBAM/Projects/motion_capture/Standard_Tests/Test_1/environmental_state");
    data.states_graph.Generate_Levels();
    FILE_UTILITIES::Write_To_File<float>(argv[1], data.states_graph);
    FILE_UTILITIES::Write_To_File<float>(argv[2], data.graph_index_to_state);
    return 0;
}

