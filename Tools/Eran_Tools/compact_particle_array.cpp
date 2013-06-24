#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <iostream>

using namespace PhysBAM;

template<class DEFORMABLE_PARTICLES,class RW>
void Compact(const std::string &filename)
{
    ARRAYS<VECTOR<DEFORMABLE_PARTICLES*,3> > particles_per_cell;
    FILE_UTILITIES::Read_From_File<RW>(filename,particles_per_cell);
    int total_cells=particles_per_cell.m*particles_per_cell.n*particles_per_cell.mn;
    int non_null_cells=0,zero_particle_cells=0;
    for(int i=particles_per_cell.m_start;i<=particles_per_cell.m_end;i++) for(int j=particles_per_cell.n_start;j<=particles_per_cell.n_end;j++) for(int k=particles_per_cell.mn_start;k<=particles_per_cell.mn_end;k++){
        if(particles_per_cell(i,j,k)){
            non_null_cells++;
            if(!particles_per_cell(i,j,k)->number){
                zero_particle_cells++;delete particles_per_cell(i,j,k);particles_per_cell(i,j,k)=0;
            }
        }
    }

    std::cout << "Cells = " << total_cells << " total, " << non_null_cells << " non-null, " << zero_particle_cells << " zero particle" << std::endl;
    std::cout << "Writing converted file..." << std::endl;
    FILE_UTILITIES::Write_To_File<RW>(filename,particles_per_cell);
}


int main(int argc,char *argv[])
{
    bool removed=false;
    PARSE_ARGS parse_args(argc,argv);
    std::string filename;

    parse_args.Add("-levelset",&opt_levelset, "levelset particle");
    parse_args.Extra(&filename,"filename","filename");
    parse_args.Parse();

    if(removed){
        std::cout << "Removed particle" << std::endl;
        Compact<PARTICLE_LEVELSET_REMOVED_PARTICLES<float,VECTOR<float,3> >,float>(filename);
    }
    else{
        std::cout << "Levelset particle" << std::endl;
        Compact<PARTICLE_LEVELSET_PARTICLES<float,VECTOR<float,3> >,float>(filename);
    }
}
