//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
using namespace PhysBAM;
//#####################################################################
// Class MESH_PARTITIONING
//#####################################################################
template<class T>
class MESH_PARTITIONING
{
    SEGMENT_MESH segment_mesh;
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume;
    RANDOM_NUMBERS<T> random_numbers;
    ARRAY<int> active_nodes;
    ARRAY<int> node_labels;
    ARRAY<ARRAY<int> > partition_nodes;
public:
    std::string output_directory;
    int partitions,epochs,iterations,cross_edge_factor,reverse_range;
    double initial_temperature,annealing_factor;

    MESH_PARTITIONING()
        :tetrahedralized_volume(0)
    {
        random_numbers.Set_Seed();
    }

    ~MESH_PARTITIONING()
    {}

    void Initialize_Segment_Mesh(const std::string& segment_mesh_filename)
    {FILE_UTILITIES::Read_From_File<T>(segment_mesh_filename,segment_mesh);segment_mesh.Initialize_Neighbor_Nodes();}

    void Initialize_Segment_Mesh_From_Geometry()
    {if(tetrahedralized_volume){
        tetrahedralized_volume->mesh.Initialize_Segment_Mesh();
        segment_mesh.Initialize_Mesh(*tetrahedralized_volume->mesh.segment_mesh);}
    else PHYSBAM_FATAL_ERROR();
    segment_mesh.Initialize_Neighbor_Nodes();}

    void Initialize_Geometry(const std::string& geometry_filename)
    {switch(FILE_UTILITIES::Get_File_Type(geometry_filename)){
        case FILE_UTILITIES::TET_FILE:
            tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create();
            FILE_UTILITIES::Read_From_File<T>(geometry_filename,*tetrahedralized_volume);
            break;
        default:
            PHYSBAM_FATAL_ERROR();}}

    void Visualize_Partitioning()
    {if(tetrahedralized_volume){
        TETRAHEDRON_MESH visualization_mesh;
        TETRAHEDRALIZED_VOLUME<T> visualization_volume(visualization_mesh,tetrahedralized_volume->particles);
        visualization_volume.Update_Number_Nodes();
        for(int t=1;t<=tetrahedralized_volume->mesh.elements.m;t++){
            VECTOR<int,4> nodes=tetrahedralized_volume->mesh.elements(t);
            VECTOR<int,4> element_labels=VECTOR<int,4>(node_labels.Subset(nodes)).Sorted();
            if(element_labels[1]==element_labels[4]) visualization_mesh.elements.Append(nodes);}
        FILE_UTILITIES::Write_To_File<T>(output_directory+"/partitioning.tet",visualization_volume);}}

//#####################################################################
// Function Flood_Fill_Mesh
//#####################################################################
int Flood_Fill_Mesh(const ARRAY<int>& centers,ARRAY<int>& labels,ARRAY<int>& distances)
{
    int functional=0;
    labels.Resize(segment_mesh.number_nodes);labels.Fill(0);
    distances.Resize(segment_mesh.number_nodes);distances.Fill(-1);
    QUEUE<int> queue(segment_mesh.number_nodes);
    for(int c=1;c<=centers.m;c++){labels(centers(c))=c;distances(centers(c))=0;queue.Enqueue(centers(c));}
    while(!queue.Empty()){
        int node=queue.Dequeue();
        for(int i=1;i<=(*segment_mesh.neighbor_nodes)(node).m;i++){
            int neighbor=(*segment_mesh.neighbor_nodes)(node)(i);
            if(!labels(neighbor)){
                labels(neighbor)=labels(node);
                distances(neighbor)=distances(node)+1;
                queue.Enqueue(neighbor);
                functional+=sqr(distances(neighbor));}}} // TODO: make this robust against overflow
    return functional;
}
//#####################################################################
// Function Compute_Initial_Clustering
//#####################################################################
void Compute_Initial_Clustering()
{
    LOG::SCOPE scope("Initial clustering");
    segment_mesh.elements.Flattened().Get_Unique(active_nodes);
    LOG::cout<<"active nodes = "<<active_nodes.m<<std::endl;
    if(partitions<1 || partitions>active_nodes.m) PHYSBAM_FATAL_ERROR();

    // Initial labeling
    LOG::Time("Initial labeling");
    ARRAY<int> centers,distances;centers.Append(active_nodes(random_numbers.Get_Uniform_Integer(1,active_nodes.m)));
    while(centers.m<partitions){Flood_Fill_Mesh(centers,node_labels,distances);centers.Append(ARRAYS_COMPUTATIONS::Argmax(distances));}
    int functional=Flood_Fill_Mesh(centers,node_labels,distances);
    for(int partition=1;partition<=partitions;partition++) LOG::cout<<"partition "<<partition<<" has "<<node_labels.Count_Matches(partition)<<" particles"<<std::endl;
    LOG::cout<<node_labels.Count_Matches(0)<<" particles do not belong to any partition"<<std::endl;
    LOG::cout<<"sum of squared radii = "<<functional<<std::endl;

    // K-means clustering
    for(int partition=1,last_optimized_partition=0;partition!=last_optimized_partition;partition=partition%partitions+1){
        for(;;){
            int old_center=centers(partition);
            for(int i=1;i<=(*segment_mesh.neighbor_nodes)(old_center).m;i++){
                centers(partition)=(*segment_mesh.neighbor_nodes)(old_center)(i);
                int new_functional=Flood_Fill_Mesh(centers,node_labels,distances);
                if(new_functional<functional){functional=new_functional;last_optimized_partition=partition;break;} else centers(partition)=old_center;}
            if(centers(partition)!=old_center) LOG::cout<<"optimized center of partition "<<partition<<", new sum of squared radii = "<<functional<<std::endl;else break;}}

    FILE_UTILITIES::Write_To_File<T>(output_directory+"/initial_node_labels",node_labels);
    Visualize_Partitioning();
}
//#####################################################################
// Function Update_Boundary_Node
//#####################################################################
void Update_Boundary_Node(ARRAY<int>& boundary_nodes,ARRAY<int>& boundary_node_index,ARRAY<int>& removed_boundary_nodes,const int node)
{
    bool is_boundary_node=false;
    for(int i=1;i<=(*segment_mesh.neighbor_nodes)(node).m;i++)
        if(node_labels((*segment_mesh.neighbor_nodes)(node)(i))!=node_labels(node)){is_boundary_node=true;break;}
    if(is_boundary_node){
        if(boundary_node_index(node)) return;
        int new_index;
        if(removed_boundary_nodes.m){new_index=removed_boundary_nodes.Pop();boundary_nodes(new_index)=node;}else new_index=boundary_nodes.Append(node);
        boundary_node_index(node)=new_index;}
    else{if(!boundary_node_index(node)) return;removed_boundary_nodes.Append(boundary_node_index(node));boundary_node_index(node)=0;}
}
//#####################################################################
// Function Optimize_Clustering
//#####################################################################
void Optimize_Clustering()
{
    LOG::SCOPE scope("Optimization");
    FILE_UTILITIES::Read_From_File<T>(output_directory+"/initial_node_labels",node_labels);partitions=ARRAYS_COMPUTATIONS::Max(node_labels);
    ARRAY<int> cluster_sizes(partitions);for(int partition=1;partition<=partitions;partition++) cluster_sizes(partition)=node_labels.Count_Matches(partition);
    int cluster_size_functional=0;
    for(int i=1;i<partitions;i++) for(int j=i+1;j<=partitions;j++) cluster_size_functional+=sqr(cluster_sizes(i)-cluster_sizes(j));
    int cross_edges=0;
    ARRAY<int> boundary_nodes,boundary_node_index(node_labels.m),removed_boundary_nodes;
    for(int s=1;s<=segment_mesh.elements.m;s++){
        VECTOR<int,2>& segment=segment_mesh.elements(s);
        if(node_labels(segment.x)!=node_labels(segment.y)){
            cross_edges++;
            if(!boundary_node_index(segment.x)) boundary_node_index(segment.x)=boundary_nodes.Append(segment.x);
            if(!boundary_node_index(segment.y)) boundary_node_index(segment.y)=boundary_nodes.Append(segment.y);}}
    int total_functional=cluster_size_functional+cross_edge_factor*cross_edges;

    {LOG::SCOPE scope("Parameters");
    LOG::cout<<"partitions           : "<<partitions<<std::endl;
    LOG::cout<<"total edges          : "<<segment_mesh.elements.m<<std::endl;
    LOG::cout<<"initial temperature  : "<<initial_temperature<<std::endl;
    LOG::cout<<"annealing factor     : "<<annealing_factor<<std::endl;
    LOG::cout<<"annealing epochs     : "<<epochs<<std::endl;
    LOG::cout<<"iterations per epoch : "<<iterations<<std::endl;
    LOG::cout<<"cross-edge factor    : "<<cross_edge_factor<<std::endl;
    LOG::cout<<"class variance       : "<<ARRAYS_COMPUTATIONS::Min(cluster_sizes)<<" - "<<ARRAYS_COMPUTATIONS::Max(cluster_sizes)<<std::endl;
    LOG::cout<<"cross-edges          : "<<cross_edges<<std::endl;
    LOG::cout<<"boundary nodes       : "<<boundary_nodes.m<<std::endl;
    LOG::cout<<"Total functional     : "<<total_functional<<std::endl;}

    double temperature=initial_temperature;
    if(boundary_nodes.m) for(int epoch=1;epoch<=epochs;epoch++){
        for(int iteration=1;iteration<=iterations;iteration++){
            int node;do{node=boundary_nodes(random_numbers.Get_Uniform_Integer(1,boundary_nodes.m));}while(!boundary_node_index(node));
            int old_label=node_labels(node);
            ARRAY<int> candidate_labels;
            for(int i=1;i<=(*segment_mesh.neighbor_nodes)(node).m;i++)
                if(node_labels((*segment_mesh.neighbor_nodes)(node)(i))!=old_label)
                    candidate_labels.Append_Unique(node_labels((*segment_mesh.neighbor_nodes)(node)(i)));
            if(!candidate_labels.m) continue;
            int new_label=candidate_labels(random_numbers.Get_Uniform_Integer(1,candidate_labels.m)),new_total_functional=total_functional,new_cross_edges=cross_edges;
            new_total_functional+=2*cluster_sizes.m*(cluster_sizes(new_label)-cluster_sizes(old_label)+1);
            for(int i=1;i<=(*segment_mesh.neighbor_nodes)(node).m;i++)
                if(node_labels((*segment_mesh.neighbor_nodes)(node)(i))==old_label){new_total_functional+=cross_edge_factor;new_cross_edges++;}
                else if(node_labels((*segment_mesh.neighbor_nodes)(node)(i))==new_label){new_total_functional-=cross_edge_factor;new_cross_edges--;}
            double threshold=exp((total_functional-new_total_functional)/temperature);
            if(new_total_functional<total_functional || threshold>random_numbers.Get_Number()){
                node_labels(node)=new_label;cluster_sizes(old_label)--;cluster_sizes(new_label)++;total_functional=new_total_functional;cross_edges=new_cross_edges;
                Update_Boundary_Node(boundary_nodes,boundary_node_index,removed_boundary_nodes,node);
                for(int i=1;i<=(*segment_mesh.neighbor_nodes)(node).m;i++)
                    Update_Boundary_Node(boundary_nodes,boundary_node_index,removed_boundary_nodes,(*segment_mesh.neighbor_nodes)(node)(i));}}
        temperature*=annealing_factor;

        if(epoch%100==0){
            LOG::SCOPE scope("Progress");
            LOG::cout<<"partitions           : "<<partitions<<std::endl;
            LOG::cout<<"temperature          : "<<temperature<<std::endl;
            LOG::cout<<"epoch                : "<<epoch<<std::endl;
            LOG::cout<<"class variance       : "<<ARRAYS_COMPUTATIONS::Min(cluster_sizes)<<" - "<<ARRAYS_COMPUTATIONS::Max(cluster_sizes)<<std::endl;
            LOG::cout<<"cross-edges          : "<<cross_edges<<std::endl;
            LOG::cout<<"boundary nodes       : "<<boundary_nodes.m-removed_boundary_nodes.m<<std::endl;
            LOG::cout<<"total functional     : "<<total_functional<<std::endl;}}
    else LOG::cout<<"Only one partition.  No need to optimize."<<std::endl;

    partition_nodes.Resize(partitions);for(int i=1;i<=node_labels.m;i++) if(node_labels(i)) partition_nodes(node_labels(i)).Append(i);
    FILE_UTILITIES::Write_To_File<T>(output_directory+"/optimized_node_labels",node_labels);
    FILE_UTILITIES::Write_To_File<T>(output_directory+"/partition_nodes",partition_nodes);
    Visualize_Partitioning();
}
//#####################################################################
// Function Log_Distance
//#####################################################################
typedef int FUNCTIONAL;
static int Log_Distance(const VECTOR<int,2>& pair)
{
    return integer_log(abs(pair.x-pair.y));
}
/*
typedef double FUNCTIONAL;
static double Log_Distance(const VECTOR<int,2>& pair)
{
    return log(abs(pair.x-pair.y))/log(2);
}
*/
//#####################################################################
// Function Optimize_Particle_Order
//#####################################################################
void Optimize_Particle_Order(const int reorder_partition_input)
{
    LOG::SCOPE scope("Reordering");
    int min_reorder_partition,max_reorder_partition;
    if(partitions){
        LOG::Time("Reading partition_nodes");
        FILE_UTILITIES::Read_From_File<T>(output_directory+"/partition_nodes",partition_nodes);
        PHYSBAM_ASSERT(partitions==partition_nodes.m);
        if(reorder_partition_input){
            if(reorder_partition_input<1 || reorder_partition_input>partitions)
                PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Invalid reorder_partition %d",reorder_partition_input));
            min_reorder_partition=max_reorder_partition=reorder_partition_input;}
        else{min_reorder_partition=1;max_reorder_partition=partitions;}}
    else{
        LOG::cout<<"Starting with default ordering."<<std::endl;
        partitions=min_reorder_partition=max_reorder_partition=1;
        partition_nodes.Resize(1);partition_nodes(1)=IDENTITY_ARRAY<>(segment_mesh.number_nodes);}

    for(int reorder_partition=min_reorder_partition;reorder_partition<=max_reorder_partition;reorder_partition++){
        ARRAY<int> local_indices(segment_mesh.number_nodes);
        for(int i=1;i<=partition_nodes(reorder_partition).m;i++) local_indices(partition_nodes(reorder_partition)(i))=i;
        SEGMENT_MESH local_segment_mesh;
        for(int s=1;s<=segment_mesh.elements.m;s++) local_segment_mesh.elements.Append(VECTOR<int,2>::Map(local_indices,segment_mesh.elements(s)));
        local_segment_mesh.Delete_Elements_With_Missing_Nodes();local_segment_mesh.number_nodes=partition_nodes(reorder_partition).m;
        local_segment_mesh.Initialize_Neighbor_Nodes();
        const ARRAY<ARRAY<int> >& neighbor_nodes=*local_segment_mesh.neighbor_nodes;
        ARRAY<int> node_locations(IDENTITY_ARRAY<>(local_segment_mesh.number_nodes));

        FUNCTIONAL functional=0;for(int s=1;s<=local_segment_mesh.elements.m;s++) functional+=Log_Distance(VECTOR<int,2>::Map(node_locations,local_segment_mesh.elements(s)));

        static const bool swap_neighbors=true; // 8.39834
        static const int swap_nearby=0;//20;
        //reverse_range 10: 8.02942, 7.89422
        STATIC_ASSERT(swap_neighbors+!!swap_nearby==1);

        LOG::SCOPE scope("Reorder partition","Reorder partition %d",reorder_partition);
        {LOG::SCOPE scope("Parameters");
        LOG::cout<<"partititions         : "<<partitions<<std::endl;
        LOG::cout<<"reordering partition : "<<reorder_partition<<std::endl;
        LOG::cout<<"nodes in partition   : "<<partition_nodes(reorder_partition).m<<std::endl;    
        LOG::cout<<"local edges          : "<<local_segment_mesh.elements.m<<std::endl;
        LOG::cout<<"initial temperature  : "<<initial_temperature<<std::endl;
        LOG::cout<<"annealing factor     : "<<annealing_factor<<std::endl;
        LOG::cout<<"annealing epochs     : "<<epochs<<std::endl;
        LOG::cout<<"iterations per epoch : "<<iterations<<std::endl;
        LOG::cout<<"functional           : "<<functional<<std::endl;
        LOG::cout<<"average log distance : "<<(double)functional/local_segment_mesh.elements.m<<std::endl;
        if(reverse_range)
            LOG::cout<<"mutations            : reverse_range = "<<reverse_range<<std::endl;
        else if(swap_neighbors)
            LOG::cout<<"mutations            : swap_neighbors"<<std::endl;
        else if(swap_nearby)
            LOG::cout<<"mutations            : swap_nearby = "<<swap_nearby<<std::endl;}

        double temperature=initial_temperature;
        for(int epoch=1;epoch<=epochs;epoch++){
            if(!reverse_range)
                for(int iteration=1;iteration<=iterations;iteration++){
                    FUNCTIONAL new_functional=functional;
                    int node1=random_numbers.Get_Uniform_Integer(1,neighbor_nodes.m),node2;
                    if(swap_neighbors)
                        node2=neighbor_nodes(node1)(random_numbers.Get_Uniform_Integer(1,neighbor_nodes(node1).m));
                    else if(swap_nearby){
                        node2=node1+random_numbers.Get_Uniform_Integer(1,swap_nearby);if(node2>neighbor_nodes.m) continue;}
                    for(int i=1;i<=neighbor_nodes(node1).m;i++){
                        int neighbor=neighbor_nodes(node1)(i);
                        if(neighbor!=node2){
                            new_functional-=Log_Distance(VECTOR<int,2>::Map(node_locations,VECTOR<int,2>(node1,neighbor)));
                            new_functional+=Log_Distance(VECTOR<int,2>::Map(node_locations,VECTOR<int,2>(node2,neighbor)));}}
                    for(int i=1;i<=neighbor_nodes(node2).m;i++){
                        int neighbor=neighbor_nodes(node2)(i);
                        if(neighbor!=node1){
                            new_functional-=Log_Distance(VECTOR<int,2>::Map(node_locations,VECTOR<int,2>(node2,neighbor)));
                            new_functional+=Log_Distance(VECTOR<int,2>::Map(node_locations,VECTOR<int,2>(node1,neighbor)));}}
                    if(new_functional<functional || exp((functional-new_functional)/temperature)>random_numbers.Get_Number()){
                        exchange(node_locations(node1),node_locations(node2));functional=new_functional;}}
            else // reverse
                for(int iteration=1;iteration<=iterations;iteration++){
                    FUNCTIONAL new_functional=functional;
                    int start=random_numbers.Get_Uniform_Integer(1,neighbor_nodes.m);
                    int end=start+random_numbers.Get_Uniform_Integer(1,reverse_range);
                    if(end>neighbor_nodes.m) continue;
                    for(int n1=start;n1<=end;n1++)for(int a=1;a<=neighbor_nodes(n1).m;a++){int n2=neighbor_nodes(n1)(a);
                        if(n2<start || n2>end){
                            new_functional-=Log_Distance(VECTOR<int,2>::Map(node_locations,VECTOR<int,2>(n1,n2)));
                            new_functional+=Log_Distance(VECTOR<int,2>::Map(node_locations,VECTOR<int,2>(start+end-n1,n2)));}
                        else if(n1<n2){
                            new_functional-=Log_Distance(VECTOR<int,2>::Map(node_locations,VECTOR<int,2>(n1,n2)));
                            new_functional+=Log_Distance(VECTOR<int,2>::Map(node_locations,VECTOR<int,2>(start+end-n1,start+end-n2)));}}
                    if(new_functional<functional || exp((functional-new_functional)/temperature)>random_numbers.Get_Number()){
                        for(int i=start,j=end;i<j;i++,j--) exchange(node_locations(i),node_locations(j));
                        functional=new_functional;}}
            temperature*=annealing_factor;
            if(epoch%100==0){
                LOG::SCOPE scope("Progress");
                LOG::cout<<"partitions           : "<<partitions<<std::endl;
                LOG::cout<<"temperature          : "<<temperature<<std::endl;
                LOG::cout<<"epoch                : "<<epoch<<std::endl;
                LOG::cout<<"functional           : "<<functional<<std::endl;
                LOG::cout<<"average log distance : "<<(double)functional/local_segment_mesh.elements.m<<std::endl;}}
        partition_nodes(reorder_partition).Subset(node_locations)=ARRAY<int>(partition_nodes(reorder_partition));}
    LOG::Time("Writing partition_nodes");
    FILE_UTILITIES::Write_To_File<T>(output_directory+"/partition_nodes",partition_nodes);
}
//#####################################################################
};
//#####################################################################
// Function main
//#####################################################################
int main(int argc,char* argv[])
{
    PARSE_ARGS parse_args;
    parse_args.Add_Option_Argument("-all","-cluster, -optimize, and -reorder");
    parse_args.Add_Option_Argument("-cluster","compute initial clustering");
    parse_args.Add_Option_Argument("-optimize","optimize clustering");
    parse_args.Add_Option_Argument("-reorder","optimize particle order");
    parse_args.Add_String_Argument("-geometry","","geometry file");
    parse_args.Add_String_Argument("-mesh","","connectivity segment mesh");
    parse_args.Add_String_Argument("-o","partitioning_data","output directory");
    parse_args.Add_Integer_Argument("-partitions",0,"partitions","partitions to use");
    parse_args.Add_Integer_Argument("-reorder_partition",0,"partition","partition to reorder (defaults to all)");
    parse_args.Add_Integer_Argument("-epochs",10000,"epochs","epochs of simulated annealing");
    parse_args.Add_Integer_Argument("-iterations",1000,"iterations","iterations per annealing epoch");
    parse_args.Add_Double_Argument("-initial_temperature",100,"temperature","initial annealing temperature");
    parse_args.Add_Double_Argument("-annealing_factor",.999,"factor","annealing factor");
    parse_args.Add_Integer_Argument("-cross_edge_factor",1,"factor","cross edge functional multiplier");
    parse_args.Add_Integer_Argument("-reverse",0,"range","use reverse mutations for reordering with given range (e.g. 10) (default is swap neighbors)");
    parse_args.Parse(argc,argv);

    LOG::Initialize_Logging();
    MESH_PARTITIONING<float> mesh_partitioning;

    mesh_partitioning.output_directory=parse_args.Get_String_Value("-o");
    FILE_UTILITIES::Create_Directory(mesh_partitioning.output_directory);
    LOG::Instance()->Copy_Log_To_File(mesh_partitioning.output_directory+"/log.txt",false);

    {LOG::SCOPE scope("Initialization");
    if(parse_args.Is_Value_Set("-geometry")) mesh_partitioning.Initialize_Geometry(parse_args.Get_String_Value("-geometry"));
    if(parse_args.Is_Value_Set("-mesh")) mesh_partitioning.Initialize_Segment_Mesh(parse_args.Get_String_Value("-mesh"));
    else if(parse_args.Is_Value_Set("-geometry")) mesh_partitioning.Initialize_Segment_Mesh_From_Geometry();
    else{LOG::cerr<<"Must specify either -geometry or -mesh."<<std::endl;parse_args.Print_Usage(true);}}

    if(!parse_args.Is_Value_Set("-partitions")){LOG::cerr<<"Must specify -partitions."<<std::endl;parse_args.Print_Usage(true);}
    mesh_partitioning.partitions=parse_args.Get_Integer_Value("-partitions");
    int reorder_partition=parse_args.Get_Integer_Value("-reorder_partition");
    mesh_partitioning.epochs=parse_args.Get_Integer_Value("-epochs");
    mesh_partitioning.iterations=parse_args.Get_Integer_Value("-iterations");
    mesh_partitioning.initial_temperature=parse_args.Get_Double_Value("-initial_temperature");
    mesh_partitioning.annealing_factor=parse_args.Get_Double_Value("-annealing_factor");    
    mesh_partitioning.cross_edge_factor=parse_args.Get_Integer_Value("-cross_edge_factor");    
    mesh_partitioning.reverse_range=parse_args.Get_Integer_Value("-reverse");

    bool all=parse_args.Is_Value_Set("-all");
    bool cluster=parse_args.Is_Value_Set("-cluster"),optimize=parse_args.Is_Value_Set("-optimize"),reorder=parse_args.Is_Value_Set("-reorder");
    if(all+cluster+optimize+reorder<1){LOG::cerr<<"Must specify one of -all, -cluster, -optimize, or -reorder."<<std::endl;parse_args.Print_Usage();}
    if(all) cluster=optimize=reorder=true;
    if(cluster && !optimize && reorder){LOG::cerr<<"Can't specify -cluster and -reorder without -optimize."<<std::endl;parse_args.Print_Usage();}

    if(cluster) mesh_partitioning.Compute_Initial_Clustering();
    if(optimize) mesh_partitioning.Optimize_Clustering();
    if(reorder) mesh_partitioning.Optimize_Particle_Order(reorder_partition);

    LOG::Finish_Logging();
    return 0;
}
//#################################################################
