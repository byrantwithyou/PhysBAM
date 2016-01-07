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
#include <regex>

using namespace PhysBAM;

template<class T,int N>
void writePartio(const std::string& output_filename,const MPM_PARTICLES<VECTOR<T,N> >& particles,bool dump_valid)
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

        VECTOR<T,3> pFx(particles.F(i).Column(0));
        VECTOR<T,3> pFy(particles.F(i).Column(1));
        VECTOR<T,3> pFz(N==3?VECTOR<T,3>(particles.F(i).Column(N-1)):VECTOR<T,3>(0,0,1));
        VECTOR<T,3> pX(particles.X(i));
        VECTOR<T,3> pV(particles.V(i));
        for(int k=0;k<3;++k){
            Fx[k]=pFx(k);
            Fy[k]=pFy(k);
            Fz[k]=pFz(k);

            p[k]=pX(k);
            v[k]=pV(k);}

        vol[0]=particles.volume(i);
        m[0]=particles.mass(i);

        if(dump_valid) *parts->dataWrite<int>(validH,idx)=particles.valid(i);}

    Partio::write(output_filename.c_str(),*parts);
    parts->release();

    LOG::printf("Exported %d particles to %s\n",particles.number,output_filename);
}

template<class T,int N>
void Convert(const std::string& input_directory,const std::string& output_filename_pattern,int start_at,bool dump_valid,bool attentive)
{
    MPM_PARTICLES<VECTOR<T,N> > particles;
    int last_frame;
    FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);
    PHYSBAM_ASSERT(start_at>=0 && start_at<=last_frame);
    for(int i=start_at;i<=last_frame;++i){
        FILE_UTILITIES::Read_From_File<T>(LOG::sprintf("%s/%d/mpm_particles.gz",input_directory,i),particles);
        writePartio<T,N>(LOG::sprintf(output_filename_pattern.c_str(),i),particles,dump_valid);
        if(attentive)
            FILE_UTILITIES::Read_From_Text_File(input_directory+"/common/last_frame",last_frame);}
}

int main(int argc,char *argv[])
{
    bool type_double=true;
    bool dump_valid=false;
    bool use_3d=false;
    bool attentive=true;
    int start_at=0;
    std::string input_directory,output_filename_pattern;

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Add("-dump_valid_flag",&dump_valid,"Dump all particles along with a valid flag; don't prune invalid particles");
    parse_args.Add("-3d",&use_3d,"Convert 3D examples");
    parse_args.Add("-start",&start_at,"start_at","Start conversion from this frame");
    parse_args.Add_Not("-no_attentive",&attentive,"Don't look for new frames during conversion");
    parse_args.Extra(&input_directory,"sim-dir","simulation output directory");
    parse_args.Extra(&output_filename_pattern,"output-pattern","output partio file name pattern, use %d");
    parse_args.Parse();

    if(!std::regex_match(output_filename_pattern,std::regex(".*%[0-9]*d.*"))){
        LOG::printf("Format string missing from filename! Eg. %%04d\n");
        return -1;}

    if(type_double)
        if(use_3d) Convert<double,3>(input_directory,output_filename_pattern,start_at,dump_valid,attentive);
        else Convert<double,2>(input_directory,output_filename_pattern,start_at,dump_valid,attentive);
    else
        if(use_3d) Convert<float,3>(input_directory,output_filename_pattern,start_at,dump_valid,attentive);
        else Convert<float,2>(input_directory,output_filename_pattern,start_at,dump_valid,attentive);

    return 0;
}
