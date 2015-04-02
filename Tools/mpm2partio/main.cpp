//#####################################################################
// Copyright 2015, Andre Pradhana and Gergely Klar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
#include <Partio.h>

using namespace PhysBAM;

template<class T>
void writePartio(const std::string& output_filename,const MPM_PARTICLES<VECTOR<T,3> >& particles,bool dump_valid)
{
    Partio::ParticlesDataMutable* parts=Partio::create();

    Partio::ParticleAttribute posH,vH,FxH,FyH,FzH,mH,volH,validH;
    posH=parts->addAttribute("position",Partio::VECTOR,3);
    vH=parts->addAttribute("v",Partio::VECTOR,3);
    FxH=parts->addAttribute("Fx",Partio::VECTOR,3);
    FyH=parts->addAttribute("Fy",Partio::VECTOR,3);
    FzH=parts->addAttribute("Fz",Partio::VECTOR,3);
    volH=parts->addAttribute("vol",Partio::VECTOR,1);
    mH=parts->addAttribute("m",Partio::VECTOR,1);
    if(dump_valid) validH=parts->addAttribute("valid",Partio::INT,1);

    for(int i=0;i<particles.number;++i){
        if(!dump_valid && !particles.valid(i)) continue;
        int idx=parts->addParticle();
        float *p=parts->dataWrite<float>(posH,idx);
        float *v=parts->dataWrite<float>(vH,idx);
        float *Fx=parts->dataWrite<float>(FxH,idx);
        float *Fy=parts->dataWrite<float>(FyH,idx);
        float *Fz=parts->dataWrite<float>(FzH,idx);
        float *vol=parts->dataWrite<float>(volH,idx);
        float *m=parts->dataWrite<float>(mH,idx);

        MATRIX<T,3> F=particles.F(i);
        for(int k=0;k<3;++k) Fx[k]=F(k,0);
        for(int k=0;k<3;++k) Fy[k]=F(k,1);
        for(int k=0;k<3;++k) Fz[k]=F(k,2);
        vol[0]=particles.volume(i);

        for(int k=0;k<3;++k) p[k]=particles.X(i)(k);
        for(int k=0;k<3;++k) v[k]=particles.V(i)(k);

        m[0]=particles.mass(i);

        if(dump_valid) *parts->dataWrite<int>(validH,idx)=particles.valid(i);}

    Partio::write(output_filename.c_str(),*parts);
    parts->release();

    LOG::printf("Exported %d particles to %s\n",particles.number,output_filename);
}

template<class T>
void Convert(const std::string& input_directory,const std::string& output_filename_pattern,bool dump_valid)
{
    MPM_PARTICLES<VECTOR<T,3> > particles;
    int last_frame;
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);
    for(int i=0;i<=last_frame;++i){
        FILE_UTILITIES::Read_From_File<T>(LOG::sprintf("%s/%d/mpm_particles.gz",input_directory,i),particles);
        writePartio<T>(LOG::sprintf(output_filename_pattern.c_str(),i),particles,dump_valid);}
}

int main(int argc,char *argv[])
{
    bool type_double=true;
    bool dump_valid=false;
    std::string input_directory,output_filename_pattern;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-dump_valid_flag",&dump_valid,"Dump all particles along with a valid flag; don't prune invalid particles");
    parse_args.Extra(&input_directory,"sim-dir","simulation output directory");
    parse_args.Extra(&output_filename_pattern,"output-pattern","output partio file name pattern, use %d");
    parse_args.Parse();

    if(type_double){
        Convert<double>(input_directory,output_filename_pattern,dump_valid);}
    else{
        Convert<float>(input_directory,output_filename_pattern,dump_valid);}

    return 0;
}
