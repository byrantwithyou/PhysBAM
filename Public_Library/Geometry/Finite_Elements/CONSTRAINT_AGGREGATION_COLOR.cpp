//#####################################################################
// Copyright 2012, Russell Howes.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTRAINT_AGGREGATION_COLOR
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Finite_Elements/CELL_DOMAIN_INTERFACE_COLOR.h>
#include <Geometry/Finite_Elements/CONSTRAINT_AGGREGATION_COLOR.h>
#include <Geometry/Finite_Elements/INTERFACE_POISSON_SYSTEM_VECTOR_COLOR.h>


using namespace PhysBAM;
//#####################################################################
// Function Get_Neighboring_Cells_From_Padded_Node
//#####################################################################
template<class TV> void CONSTRAINT_AGGREGATION_COLOR<TV>::
    Get_Neighboring_Cells_From_Padded_Node(int node,int cube_radius,ARRAY<int>& neighbors)
{   
    assert(cube_radius==2||cube_radius==1);
    assert(neighbors.m==1<<(cube_radius*TV::m));
    
    TV_INT padded_cells;for(int i=0;i<TV::m;i++)padded_cells(i)=2*cdi->padding;padded_cells+=grid.numbers_of_cells;
    TV_INT node_index=Index_Of_Uncompressed_Node(node),a;
    
    a(TV::m-1)=1;for(int i=TV::m-2;i>=0;i--)a(i)=a(i+1)*grid.numbers_of_cells(i+1);
    
    int this_neighbor=0;
    TV_INT counts=TV_INT()+2*cube_radius;
    
    GRID<TV> grid2(counts,RANGE<TV>(TV(),TV()+1),true);
    
    for(CELL_ITERATOR<TV> it(grid2);it.Valid();it.Next()){
        TV_INT this_neighbors_address=node_index+it.index-cube_radius;
        for(int i=0;i<TV::m;i++){
            if(this_neighbors_address(i)<0)this_neighbors_address(i)+=grid.numbers_of_cells(i);
            if(this_neighbors_address(i)>=grid.numbers_of_cells(i))this_neighbors_address(i)-=grid.numbers_of_cells(i);
        }
        neighbors(this_neighbor++)=this_neighbors_address.Dot(a);
    }
}
//#####################################################################
// Function Distance_Between_Node_And_Cell
//#####################################################################
template<class TV> int CONSTRAINT_AGGREGATION_COLOR<TV>::
Distance_Between_Node_And_Cell(int node,int cell)
{   
    
    TV_INT padded_cells;for(int i=0;i<TV::m;i++)padded_cells(i)=2*cdi->padding;padded_cells+=grid.numbers_of_cells;
    TV_INT node_index=Index_Of_Uncompressed_Node(node),a;
    
    a(TV::m-1)=1;for(int i=TV::m-2;i>=0;i--)a(i)=a(i+1)*grid.numbers_of_cells(i+1);
    
    int cell_temp=cell;
    TV_INT cell_index;for(int i=TV::m-1;i>=0;i--){cell_index(i)=cell_temp % grid.numbers_of_cells(i);cell_temp-=cell_index(i);cell_temp/=grid.numbers_of_cells(i);}
    
    TV_INT distance=cell_index*2-node_index*2+1;
    for(int i=0;i<TV::m;i++){if(distance(i)<0)distance(i)*=-1;if(distance(i)>grid.numbers_of_cells(i))distance(i)=2*grid.numbers_of_cells(i)-distance(i);}
    distance-=1;distance/=2;
    int total_distance=distance.Sum();
    
    return total_distance;
}
//#####################################################################
// Function Index_Of_Uncompressed_Node
//#####################################################################
template<class TV> VECTOR<int,TV::m> CONSTRAINT_AGGREGATION_COLOR<TV>::
Index_Of_Uncompressed_Node(int node)
{   
    
    int node_temp=node; TV_INT padded_cells;for(int i=0;i<TV::m;i++)padded_cells(i)=2*cdi->padding;padded_cells+=grid.numbers_of_cells;
    TV_INT node_index;for(int i=TV::m-1;i>=0;i--){node_index(i)=node_temp % padded_cells(i);node_temp-=node_index(i);node_temp/=padded_cells(i);}
    //Remove the rest of the padding
    for(int i=0;i<TV::m;i++)node_index(i)--;
    return node_index;
}
//#####################################################################
// Function Aggregate_Constraints
//#####################################################################
template<class TV> void CONSTRAINT_AGGREGATION_COLOR<TV>::
Aggregate_Constraints(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix_uu,const ARRAY<int,TV_INT>& phi_color)
{

    int all_dofs(cm->uncompressed.m);
    int total_cells(grid.numbers_of_cells.Product());
    weights.Resize(all_dofs);
    indices.Resize(all_dofs);for(int i=0;i<all_dofs;i++)indices(i)=i;
    is_idof_combined.Resize(all_dofs,true,true,-1);
    is_idof.Resize(cdi->colors);
    is_non_idof.Resize(cdi->colors);
    int this_index=0;
    ARRAY<int> nodes_per_color(cdi->colors);
    ARRAY<int> carried_nodes(cdi->colors);
    ARRAY<int> neighboring_nodes(1<<(2*TV::m)),adjacent_nodes(1<<TV::m);
    ARRAY<VECTOR<int,2> > indices_represent; indices_represent.Resize(all_dofs);
    
    Set_Neighbor_Node_Offsets();
    
    ARRAY<ARRAY<VECTOR<int,2> > >& accp=biu.all_constraint_color_pairs;
    //ARRAY<ARRAY<int> >& aco=biu.all_constraint_offsets; //Commented out until we use it
    
    for(int c=0;c<cdi->colors;c++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=full_constraint_matrix(c);
        is_idof(c).Resize(cm->dofs(c),true,true,-1);//Set all to not IDOF
        is_non_idof(c).Resize(cm->dofs(c),true,true,-1);//Set all to IDOF
        //Go through all DOFs and fill in a vector of weights
        
        for(int j=0;j<M.m;j++)for(int k=0;k<M.offsets(j+1)-M.offsets(j);k++){
            int column_index_jk=M.A(M.offsets(j)+k).j;
            weights(this_index+column_index_jk)+=fabs(M.A(M.offsets(j)+k).a);
        }
        for(int j=0;j<cm->dofs(c);j++){indices_represent(this_index+j)=VECTOR<int,2>(c,j);}
        carried_nodes(c)=this_index;this_index+=cm->dofs(c);nodes_per_color(c)=cm->dofs(c);
    }
    
    //for(int i=0;i<all_dofs;i++)std::cout<<i<<" "<<weights(i)<<" "<<std::endl;
    //for(int i=0;i<all_dofs;i++)std::cout<<indices(i)<<" "<<cm->uncompressed(i)(1)<<" "<<weights(i)<<std::endl;
    //Set weights to be zero for non-virtual nodes
//    for(CELL_ITERATOR<TV> it(grid);it.Valid();it.Next()){
//        int c=phi_color(it.index);
//        if(c>=0){
//            std::cout<<it.index<<" "<<it.Location()<<" "<<c<<std::endl;
//            int k=cm->Get_Index(it.index,c);
//            assert(k>=0);
//            weights(carried_nodes(c)+k)=(T)-1;}}
    
    // Initialize two arrays of ints or bools: one saying whether the DOFs are adjacent to some IDOF, another saying whether each bdy cell is within a 4x4(x3) block of grid cells centered on an IDOF in I (these are both from Algorithm 1 of Jeff's paper)
    // Here's what we do with cell_checker. Initially cell_checker will have a value of 1 for a non-cut cell, or 2 (or 3 or 4, etc., but probably 2) for a cut cell. Once the algorithm starts, all cells within the block of the new IDOF will be reduced to zero and hit_cells will be decremented.
    ARRAY<int> dof_adj_to_an_idof(all_dofs),cell_checker; cell_checker.Resize(total_cells,true,true,0);
    int hit_cells(0);
        for(int i=0;i<total_cells;i++)if(accp(i).m)
            {cell_checker(i)=accp(i).m;hit_cells++;
                if(cell_checker(i)>1)//If we are at a triple junction we can't allow any of the incident nodes to be IDOF
                {}//Do this later
            }

    //Order the dofs by their weight
    //indices.Sort([&](int a,int b){return weights(a)<weights(b);};);
    //for(int i=0;i<all_dofs;i++)std::cout<<weights(i)<<" "<<indices(i)<<" "<<cm->uncompressed(indices(i))(1)<<" "<<Index_Of_Uncompressed_Node(cm->uncompressed(indices(i))(1))<<std::endl;
    
    //Ugly temporary sort, to be replaced by the thing above ---^ at some point...this works just fine for now
    ARRAY<VECTOR<T,2> >weights_and_indices; weights_and_indices.Resize(all_dofs);
    for(int i=0;i<all_dofs;i++){weights_and_indices(i)(0)=weights(i);weights_and_indices(i)(1)=(T)i+1e-8;}
    weights_and_indices.Sort(LEXICOGRAPHIC_COMPARE());
    for(int i=0;i<all_dofs;i++){ indices(i)=(int)(weights_and_indices(i)(1));weights(i)=weights_and_indices(i)(0);}
    
    // Loop through the virtual nodes
    for(int i=all_dofs-1;i>=0;i--){
        if(weights(i)<=1e-8) break; //We should no longer get to the point of dealing with zero weight nodes
        if(!dof_adj_to_an_idof(indices(i))){
            idofs.Append(indices(i)); //Add the next independent degree of freedom
            is_idof(indices_represent(indices(i))(0))(indices_represent(indices(i))(1))=idofs.m-1;
            is_idof_combined(indices(i))=idofs.m-1;
            //Fill in the dofs adjacent to this one (for all colors!)
            int this_uncompressed_index=cm->uncompressed(indices(i))(1);
            for(int c=0;c<cdi->colors;c++)
                for(int j=0;j<neighbor_node_offsets.m;j++){
                    int neighbor_node=cm->compressed(c)(this_uncompressed_index+neighbor_node_offsets(j));
                    if(neighbor_node>=0){
                        dof_adj_to_an_idof(carried_nodes(c)+neighbor_node)=1;//carried_nodes(c)+cm->compressed(this_uncompressed_index+neighbor_node_offsets(j)))=1;
                    }
                }
            
            //Fill in the cells within a 4x4x4 block of this dof
            Get_Neighboring_Cells_From_Padded_Node(this_uncompressed_index,2,neighboring_nodes);
            //Check to see if all the cells are in the 4x4x4 block yet
            for(int j=0;j<neighboring_nodes.m;j++){if(cell_checker(neighboring_nodes(j))>0){cell_checker(neighboring_nodes(j))*=-1;hit_cells--;}}
            if(hit_cells<0)PHYSBAM_FATAL_ERROR(); if(hit_cells==0)break;
        }
        
    }    
    //Now that we have our independent degrees of freedom, we decide which individual constraints we will associate with
    //each IDOF
    
    //***TODO: In the future we need to assign all the idofs to a color pair based on weighting
    idof_color_pairs.Resize(idofs.m);

    corresponding_singlewide_constraints.Resize(idofs.m);
    
    
    //Fill this array (of arrays, one for each color pair) with the singlewide constraints indexes corresponding to the
    //color pair and location. -1 means not used.
    ARRAY<ARRAY<int> > singlewide_constraint_to_be_picked; singlewide_constraint_to_be_picked.Append(ARRAY<int>());
    singlewide_constraint_to_be_picked(0).Resize(total_cells,true,true,-1);//***Update this for when we worry abt triple junctions
    
    int next_index=0;
    for(int i=0;i<accp.m;i++)for(int j=0;j<accp(i).m;j++){
        if(j>0)PHYSBAM_FATAL_ERROR();//We shouldn't have more than one constraint per cell yet ***TODO: Fix this
        int this_color_pair=0;
        singlewide_constraint_to_be_picked(this_color_pair)(i)=(next_index++);
    }
    
    
//    //Go through each IDOF and get the cells adjacent and claim them
//    for(int i=0;i<idofs.m;i++){
//        Get_Neighboring_Cells_From_Padded_Node(cm->uncompressed(idofs(i))(1),1,adjacent_nodes);
//        for(int j=0;j<adjacent_nodes.m;j++){
//            int k=singlewide_constraint_to_be_picked(0)(adjacent_nodes(j));if(k>=0){corresponding_singlewide_constraints(i).Append(k);singlewide_constraint_to_be_picked(0)(adjacent_nodes(j))=-1;}
//        }
//            
//    }
    
    //For all the unclaimed cells, attach them to the closest IDOF(***TODO: Of that same color)
    for(int i=0;i<singlewide_constraint_to_be_picked.m;i++)for(int j=0;j<singlewide_constraint_to_be_picked(i).m;j++)
    {//i is the color pair
        int k = singlewide_constraint_to_be_picked(i)(j); if(k>=0){
            int closest_dof=-1;int closest_distance=grid.numbers_of_cells.Sum();//Should be enough right? lol
            for(int l=0;l<idofs.m;l++)if(true){//if the IDOF is the right color pair
                int d0=Distance_Between_Node_And_Cell(cm->uncompressed(idofs(l))(1),j);
                if(d0<closest_distance){closest_distance=d0;closest_dof=l;}
                }
        //corres...reset k to 0. we don't need to do that right?
        corresponding_singlewide_constraints(closest_dof).Append(k);
        }
    }
    
//        for(int i=0;i<corresponding_singlewide_constraints.m;i++){std::cout<<idofs(i)<<"("<<cm->uncompressed(idofs(i))(0)<<","<<Index_Of_Uncompressed_Node(cm->uncompressed(idofs(i))(1))<<"): ";
//            for(int j=0;j<corresponding_singlewide_constraints(i).m;j++){
//                std::cout<<corresponding_singlewide_constraints(i)(j)<<" ";
//            }std::cout<<std::endl;
//        }
    
    //We need to get the addresses of all the non-IDOFS
    int this_non_idof=0;    
    for(int c=0;c<cdi->colors;c++)for(int i=0;i<is_non_idof(c).m;i++)if(is_idof(c)(i)<0)is_non_idof(c)(i)=this_non_idof++;

}
//#####################################################################
// Function Build_Condensed_Constraint_Matrix
//#####################################################################
template<class TV> void CONSTRAINT_AGGREGATION_COLOR<TV>::
Build_Condensed_Constraint_Matrix(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& matrix)
{
    matrix.Resize(cdi->colors);
    
    int m=idofs.m;
    agg_constraint_rhs.Resize(m);
    
    //This thing doesn't handle multiple color pairs at all just yet...needs to be written in
    
    for(int c=0;c<cdi->colors;c++){
        SPARSE_MATRIX_FLAT_MXN<T>& M_old=full_constraint_matrix(c);
        SPARSE_MATRIX_FLAT_MXN<T>& M_new=matrix(c);
        int n=cm->dofs(c);
        M_new.Reset(n);
        M_new.offsets.Resize(m+1);
        M_new.m=m;
        
        for(int i=0;i<m;i++){
            ARRAY<SPARSE_MATRIX_ENTRY<T> > condensed_row_entries;
            ARRAY<T> condensed_row_values;condensed_row_values.Resize(n);
            for(int j=0;j<corresponding_singlewide_constraints(i).m;j++){
                int this_full_row_number=corresponding_singlewide_constraints(i)(j);
                if(c==0)agg_constraint_rhs(i)+=full_constraint_rhs(this_full_row_number);//Let's not take this out of the loop. Lazy
                ARRAY<SPARSE_MATRIX_ENTRY<T> > these_full_row_entries;
                M_old.Get_Row(these_full_row_entries,this_full_row_number);
                for(int k=0;k<these_full_row_entries.m;k++){
                    condensed_row_values(these_full_row_entries(k).j)+=these_full_row_entries(k).a;
                }
            }
            for(int j=0;j<n;j++){
                T value=condensed_row_values(j);
                if(value){M_new.offsets(i+1)++;condensed_row_entries.Append(SPARSE_MATRIX_ENTRY<T>(j,value));}
            }
            M_new.A.Append_Elements(condensed_row_entries);
        }
        for(int i=0;i<M_new.offsets.m-1;i++){ M_new.offsets(i+1)+=M_new.offsets(i);}
    }
    condensed_constraint_matrix=matrix;
}
//#####################################################################
// Function Build_Nullspace_And_Specific_Constraint_Solution
//#####################################################################
template<class TV> void CONSTRAINT_AGGREGATION_COLOR<TV>::
Build_Nullspace_And_Specific_Constraint_Solution(ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& a_matrix,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& z,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& z_transpose,ARRAY<SPARSE_MATRIX_FLAT_MXN<T> >& z_t_a,SPARSE_MATRIX_FLAT_MXN<T>& ztaz,VECTOR_T& specific_solution)
{

    int n_minus_m=cm->uncompressed.m-idofs.m;
    
    specific_solution.Zero_Out();

    z.Resize(cdi->colors);
    int this_non_idof=0;
    for(int c=0;c<cdi->colors;c++){
        SPARSE_MATRIX_FLAT_MXN<T>& M=z(c);M.Reset(n_minus_m);
        int offsets=0;
        int mc = cm->dofs(c);
        M.offsets.Resize(mc+1);M.offsets(0)=0;
        M.m=mc;
        for(int i=0;i<mc;i++){
            if(is_idof(c)(i)>-1)
            {
                T diagonal_term=condensed_constraint_matrix(c)(is_idof(c)(i),i);
                specific_solution.u(c)(i)=agg_constraint_rhs(is_idof(c)(i))/diagonal_term;
                //Copy the constraint row corresponding to this IDOF. We need to combine several different constraint matrices
                for(int c0=0;c0<cdi->colors;c0++){
                    SPARSE_MATRIX_FLAT_MXN<T>& M_constraints_c0=condensed_constraint_matrix(c0);
                    ARRAY<SPARSE_MATRIX_ENTRY<T> >this_row;
                    M_constraints_c0.Get_Row(this_row,is_idof(c)(i));
                    for(int j=0;j<this_row.m;j++){
                        if(is_non_idof(c0)(this_row(j).j)<0){this_row.Remove_Index(j);j--;}
                        else{this_row(j).j=is_non_idof(c0)(this_row(j).j);
                        this_row(j).a/=-diagonal_term;//Left multiply by -Bm(inverse), equation 20 of Jeff's paper
                        }
                    }
                    M.A.Append_Elements(this_row);
                    offsets+=this_row.m;
                }
                M.offsets(i+1)=offsets;
            }
            else{
                //Put the 1 where it belongs
                M.A.Append(SPARSE_MATRIX_ENTRY<T>(this_non_idof++,(T)1));
                offsets++;M.offsets(i+1)=offsets;
            }
        }
    }
    z_transpose.Resize(cdi->colors);
#pragma omp parallel for
    for(int c=0;c<cdi->colors;c++)
        z(c).Transpose(z_transpose(c));
    
    //Now we build AZ..actually with the functions available to use let's build Z'A
    z_t_a.Resize(cdi->colors);
    for(int c=0;c<cdi->colors;c++){
        z_t_a(c)=z_transpose(c)*(a_matrix(c));//Because A better be symmetric...
    }
    
    //Now we build Z'AZ
    ztaz.Reset(n_minus_m);
    ztaz.m=n_minus_m;
    ztaz.offsets.Resize(ztaz.m+1);
    
    for(int c=0;c<cdi->colors;c++){
        ARRAY<T> ones_diag;ones_diag.Resize(z(c).m,true,true,(T)1); //Because I can't find a regular matrix multiply
        ztaz=ztaz+z_t_a(c)*(z(c));
    }
//    OCTAVE_OUTPUT<T>("z0.txt").Write("z0",z(0));//,*vectors(0),*vectors(1));
//    OCTAVE_OUTPUT<T>("z1.txt").Write("z1",z(1));//,*vectors(0),*vectors(1));
//    OCTAVE_OUTPUT<T>("b1.txt").Write("b1",condensed_constraint_matrix(0));//,*vectors(0),*vectors(1));
//    OCTAVE_OUTPUT<T>("b2.txt").Write("b2",condensed_constraint_matrix(1));
//    OCTAVE_OUTPUT<T>("b10.txt").Write("b10",full_constraint_matrix(0));//,*vectors(0),*vectors(1));
//    OCTAVE_OUTPUT<T>("b20.txt").Write("b20",full_constraint_matrix(1));
//    OCTAVE_OUTPUT<T>("a1.txt").Write("a1",a_matrix(0));
//    OCTAVE_OUTPUT<T>("a2.txt").Write("a2",a_matrix(1));
//    OCTAVE_OUTPUT<T>("ztaz.txt").Write("ztaz",ztaz);
//    OCTAVE_OUTPUT<T>("c.txt").Write("c",specific_solution);
//    OCTAVE_OUTPUT<T>("agg.txt").Write("agg",agg_constraint_rhs);
}
namespace PhysBAM{
template class CONSTRAINT_AGGREGATION_COLOR<VECTOR<float,2> >;
template class CONSTRAINT_AGGREGATION_COLOR<VECTOR<float,3> >;
template class CONSTRAINT_AGGREGATION_COLOR<VECTOR<double,2> >;
template class CONSTRAINT_AGGREGATION_COLOR<VECTOR<double,3> >;
}
