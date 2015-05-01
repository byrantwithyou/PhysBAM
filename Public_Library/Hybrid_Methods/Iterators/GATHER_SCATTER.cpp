//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Vectors/VECTOR.h>
#include <Hybrid_Methods/Iterators/GATHER_SCATTER.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_ITERATOR.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <map>
#include <queue>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GATHER_SCATTER<TV>::
GATHER_SCATTER(const GRID<TV>& grid,const ARRAY<int>& simulated_particles)
    :simulated_particles(simulated_particles),weights(0),grid(grid),threads(1)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GATHER_SCATTER<TV>::
~GATHER_SCATTER()
{
}
//#####################################################################
// Function Prepare_Scatter
//#####################################################################
template<class TV> void GATHER_SCATTER<TV>::
Prepare_Scatter(const MPM_PARTICLES<TV>& particles)
{
    if(threads<2) return;
    int min_width=(weights->Order()+1)/2;/(min_width*partitions));
    ARRAY<int> counts(grid.numbers_of_cells.x),sum_counts(grid.numbers_of_cells.x);
    for(int k=0;k<simulated_particles.m;++k){
        int p=simulated_particles(k);
        TV_INT cell=grid.Cell(particles.X(p),0);
        counts(cell.x)++;}
    sum_counts(0)=counts(0);
    for(int i=1;i<counts.m;i++) sum_counts(i)=sum_counts(i-1)+counts(i);

    ARRAY<int> bin_ends2,bin_ends3,bin_ends;
    Compute_Optimal_Bins(bin_ends2,counts,sum_counts,2*min_width,2*threads);
    Compute_Optimal_Bins(bin_ends3,counts,sum_counts,min_width,3*threads);

    int max2[2]={0},max3[3]={0};

    max2[0]=sum_counts(bin_ends2(0));
    max3[0]=sum_counts(bin_ends3(0));
    for(int i=1;i<bin_ends2.m;i++)
        max2[i%2]=max(max2[i%2],sum_counts(bin_ends2(i))-sum_counts(bin_ends2(i-1)));
    for(int i=1;i<bin_ends3.m;i++)
        max3[i%3]=max(max3[i%3],sum_counts(bin_ends3(i))-sum_counts(bin_ends3(i-1)));
    int cost2=max2[0]+max2[1],cost3=max3[0]+max3[1]+max3[2];
    LOG::printf("COST %i %i   %i %i %i   (%i %i)\n",max2[0],max2[1],max3[0],max3[1],max3[2],cost2,cost3);
    if(cost2>cost3+10){
        LOG::printf("USING 3 STRIPES\n");
        bin_ends.Exchange(bin_ends3);
        partitions=3;}
    else{
        LOG::printf("USING 2 STRIPES\n");
        bin_ends.Exchange(bin_ends2);
        min_width*=2;
        partitions=2;}

    ARRAY<int> bin_map(counts.m);
    for(int b=1;b<bin_ends.m;b++)
        for(int i=bin_ends(b-1)+1;i<=bin_ends(b);i++)
            bin_map(i)=b;

    bins.Remove_All();
    bins.Resize(partitions*threads);
    for(int k=0;k<simulated_particles.m;++k){
        int p=simulated_particles(k);
        TV_INT cell=grid.Cell(particles.X(p),0);
        // int c=bin_map(cell.x);
        // Add_Debug_Particle(particles.X(p),VECTOR<T,3>(c&1,c/2&1,c/4&1));
        bins(bin_map(cell.x)).Append(p);}
    LOG::printf("bins:");
    for(int i=0;i<bins.m;i++)
        LOG::printf(" %i",bins(i).m);
    LOG::printf("\n");
}
namespace
{
struct COST_PAIR
{
    int value,prev_end;

    COST_PAIR():value(INT_MAX/2),prev_end(INT_MAX/2) {}
    COST_PAIR(int value,int prev_end):value(value),prev_end(prev_end) {}

    bool operator<(const COST_PAIR& c) const {return value<c.value;}
};
struct IVAL
{
    int prev_end,ival_start,ival_end,value;

    IVAL(): prev_end(INT_MAX/2),ival_start(INT_MAX/2),ival_end(INT_MAX/2),value(INT_MAX/2) {}
    IVAL(int prev_end,int ival_start,int ival_end,int value)
        :prev_end(prev_end),ival_start(ival_start),ival_end(ival_end),value(value) {}
};
struct REMOVE_EVENT
{
    int ival_end;
    std::multimap<int,IVAL>::iterator it;

    bool operator<(const REMOVE_EVENT& c) const {return ival_end>c.ival_end;}
};
}
//#####################################################################
// Function Compute_Optimal_Bins
//#####################################################################
template<class TV> int GATHER_SCATTER<TV>::
Compute_Optimal_Bins(ARRAY<int>& bin_ends,ARRAY<int>& counts,ARRAY<int>& sum_counts,int min_width,int max_bins)
{
    ARRAY<ARRAY<COST_PAIR> > table(max_bins);
    for(int i=0;i<max_bins;i++) table(i).Resize(counts.m);
    ARRAY<COST_PAIR> C0(counts.m),C1(counts.m),T0(counts.m);
    ARRAY<IVAL> next_intervals;

    C0(0)=COST_PAIR(counts(0),-1);

    for(int b=0;b<max_bins;b++){
        if(b){
            C1.Fill(COST_PAIR());
            next_intervals.Sort([](const IVAL& a,const IVAL& b){return a.ival_start<b.ival_start;});
            std::priority_queue<REMOVE_EVENT> q;
            std::multimap<int,IVAL> current_intervals;
            current_intervals.emplace(INT_MAX/2,IVAL());

            int next=0;
            for(int i=0;i<sum_counts.m;i++){
                while(next<next_intervals.m && next_intervals(next).ival_start==i){
                    IVAL iv=next_intervals(next++);
                    REMOVE_EVENT re={iv.ival_end,current_intervals.emplace(iv.value,iv)};
                    q.push(re);}

                IVAL iv=current_intervals.begin()->second;
                T0(i)={iv.value,iv.prev_end};

                while(!q.empty() && q.top().ival_end==i){
                    current_intervals.erase(q.top().it);
                    q.pop();}}
            next_intervals.Remove_All();}

        for(int i=0;i<sum_counts.m;i++){
            COST_PAIR pr=min(C0(i),T0(i));
            table(b)(i)=pr;
            if(i+1<counts.m) C0(i+1)=min(C0(i+1),COST_PAIR(C0(i).value+counts(i+1),C0(i).prev_end));
            int k=sum_counts.Binary_Search(pr.value+sum_counts(i)+1);
            IVAL iv(i,i+min_width,k-1,pr.value);
            if(iv.ival_start<=iv.ival_end) next_intervals.Append(iv);
            k=max(k,i+min_width);
            if(k<sum_counts.m) C1(k)={sum_counts(k)-sum_counts(i),i};}
        C0.Exchange(C1);}

    int best_b=-1;
    int best_cost=INT_MAX/2;
    for(int b=0;b<table.m;b++)
        if(table(b).Last().value<=best_cost){
            best_cost=table(b).Last().value;
            best_b=b;}

    bin_ends.Resize(best_b+1);
    int end=counts.m-1;
    for(int b=best_b;b>=0;b--){
        bin_ends(b)=end;
        end=table(b)(end).prev_end;}

    return best_cost;
}
template class GATHER_SCATTER<VECTOR<float,2> >;
template class GATHER_SCATTER<VECTOR<float,3> >;
template class GATHER_SCATTER<VECTOR<double,2> >;
template class GATHER_SCATTER<VECTOR<double,3> >;
}
