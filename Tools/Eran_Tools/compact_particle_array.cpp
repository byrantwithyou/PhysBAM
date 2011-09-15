#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <iostream>

using namespace PhysBAM;

template<class PARTICLES,class RW>
void Compact(const std::string &filename)
{
    ARRAYS<VECTOR<PARTICLES*,3> > particles_per_cell;
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
    PARSE_ARGS parse_args;

    parse_args.Add_Option_Argument("-levelset", "levelset particle");
    parse_args.Add_Option_Argument("-removed", "removed particle");
    parse_args.Set_Extra_Arguments(1, "<filename>");
    parse_args.Parse(argc, argv);

    std::string filename;

    if (parse_args.Num_Extra_Args() < 1) return 1;
    else filename=parse_args.Extra_Arg(1);

    if(parse_args.Get_Option_Value("-removed")){
        std::cout << "Removed particle" << std::endl;
        Compact<PARTICLE_LEVELSET_REMOVED_PARTICLES<float,VECTOR<float,3> >,float>(filename);
    }
    else{
        std::cout << "Levelset particle" << std::endl;
        Compact<PARTICLE_LEVELSET_PARTICLES<float,VECTOR<float,3> >,float>(filename);
    }
}
