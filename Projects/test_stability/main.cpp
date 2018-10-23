//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Arrays/PROJECTED_ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>

#define lapack_complex_double std::complex<double>
#include <lapacke.h>

using namespace PhysBAM;

void Compute_Eigenvalues(MATRIX_MXN<std::complex<double> > M,ARRAY<std::complex<double> >& eig)
{
    typedef double T;
    PHYSBAM_ASSERT(M.m==M.n);
    char jobvl = 'N';
    char jobvr = 'N';
    int n = M.m, size_work=10*n;
    ARRAY<std::complex<T> > work(size_work);
    ARRAY<T> rwork(2*n);
    int info = 0;

    LAPACK_zgeev(&jobvl, &jobvr, &n, &M(0,0), &n, &eig(0), 0, &n, 0, &n,
        &work(0), &size_work, &rwork(0), &info);

    if(info) printf("ZGEEV FAILED: %i\n",info);
    PHYSBAM_ASSERT(info==0);
}
void Compute_Eigenvalues(MATRIX_MXN<std::complex<float> > M,ARRAY<std::complex<float> >& eig)
{
    PHYSBAM_FATAL_ERROR();
}


template<class TV,class T>
void Step(SOLID_BODY_COLLECTION<TV>& sbc,T time,T dt,T alpha,const std::string& out_dir,int frame,int step)
{
    sbc.Update_Simulated_Particles();
    sbc.Update_Position_Based_State(time,false,false);
    DEFORMABLE_PARTICLES<TV>& particles=sbc.deformable_body_collection.particles;
    MATRIX_MXN<T> A(particles.number*TV::m,particles.number*TV::m),M(particles.number*TV::m*2,particles.number*TV::m*2);
    ARRAY<TV> V(particles.number),col(particles.number);
    auto idx=[](int p,int dim){return p*TV::m+dim;};
    int base=particles.number*TV::m;
    for(int i=0;i<particles.number;i++){
        for(int j=0;j<TV::m;j++){
            V*=0;
            V(i)(j)=1;
            col*=0;
            sbc.deformable_body_collection.Add_Implicit_Velocity_Independent_Forces(V,col,time);
            for(int k=0;k<particles.number;k++){
                for(int l=0;l<TV::m;l++){
                    T a=col(k)(l)/particles.mass(k);
                    M(idx(k,l),idx(i,j))=(idx(k,l)==idx(i,j)?1:0)+dt*dt*alpha*a;
                    M(base+idx(k,l),idx(i,j))=dt*a;
                    if(idx(k,l)==idx(i,j)){
                        M(idx(k,l),base+idx(i,j))=dt;
                        M(base+idx(k,l),base+idx(i,j))=1;}

                    A(idx(k,l),idx(i,j))=col(k)(l);}}}}

    ARRAY<T> U(particles.number*TV::m*2),U1(particles.number*TV::m*2);
    for(int i=0;i<particles.number;i++){
        for(int j=0;j<TV::m;j++){
            U(idx(i,j))=particles.X(i)(j);
            U(base+idx(i,j))=particles.V(i)(j);}}
    U1=M*U;

    OCTAVE_OUTPUT<T>(LOG::sprintf("%s/A%d.txt",out_dir,step).c_str()).Write("A",A);
    OCTAVE_OUTPUT<T>(LOG::sprintf("%s/M%d.txt",out_dir,step).c_str()).Write("M",M);
    ARRAY<std::complex<T> > eig(M.m);
    MATRIX_MXN<std::complex<T> > cM(M.m);
    for(int i=0;i<M.x.m;i++) cM.x(i)=(std::complex<T>)M.x(i);
    Compute_Eigenvalues(cM,eig);
    LOG::printf("eigenvalues: %P\n",eig);

    col*=0;
    sbc.Add_Velocity_Independent_Forces(col,ARRAY_VIEW<TWIST<TV> >(),time);
    for(int i=0;i<particles.number;i++){
        particles.X(i)+=dt*dt*alpha/particles.mass(i)*col(i)+dt*particles.V(i);
        particles.V(i)+=dt/particles.mass(i)*col(i);}

    T x_error=-1,v_error=-1;
    for(int i=0;i<particles.number;i++){
        for(int j=0;j<TV::m;j++){
            x_error=max(x_error,abs(U1(idx(i,j))-particles.X(i)(j)));
            v_error=max(v_error,abs(U1(base+idx(i,j))-particles.V(i)(j)));
        }
    }
    LOG::printf("L-inf error x: %P v: %P\n",x_error,v_error);
}

template<class TV>
void Run(PARSE_ARGS& parse_args,STREAM_TYPE stream_type,const std::string& output_dir)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    T frame_dt=1.0/24,alpha=0,dt=0.01;
    int last_frame=1;
    parse_args.Add("-frame_dt",&frame_dt,"rate","Number of frames per second");
    parse_args.Add("-dt",&dt,"time","dt");
    parse_args.Add("-last_frame",&last_frame,"frame","number of frames to simulate");
    parse_args.Add("-a",&alpha,"alpha","alpha");
    parse_args.Parse();

    std::string data_directory="../../Public_Data";
    VIEWER_OUTPUT<TV> vo(stream_type,GRID<TV>(),output_dir);

    T density=1;
    RIGID_BODY_STATE<TV> initial_state;
    SOLID_BODY_COLLECTION<TV> sbc;
    SOLIDS_STANDARD_TESTS<TV> tests(stream_type,data_directory,sbc);
    TETRAHEDRALIZED_VOLUME<T>& st=tests.Create_Tetrahedralized_Volume(
        data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
        initial_state,true,true,density,1);

    // stiffness, poissons_ratio, damping
    FINITE_VOLUME<TV,TV::m>* force=Create_Finite_Volume(st,new COROTATED_FIXED<T,TV::m>(1e5,.45,0));
    sbc.Add_Force(force);

    sbc.deformable_body_collection.particles.X.template Project<T,&TV::y>()*=0.9;
    Flush_Frame<TV>("init");
    sbc.Write(stream_type,output_dir,vo.frame-1,0,true,true,true,true,true);

    T time=0;
    int step=0;
    for(int frame=0;frame<last_frame;frame++){
        T frame_end=(frame+1)*frame_dt;
        while(time<frame_end){
            T step_dt=dt;
            T end_time=time+dt;
            if(end_time>=frame_end){
                step_dt=frame_end-time;
                end_time=frame_end;}
            Step(sbc,time,step_dt,alpha,output_dir,frame,step);
            time=end_time;
            step++;}
        Flush_Frame<TV>("Frame");
        sbc.Write(stream_type,output_dir,vo.frame-1,0,true,true,true,true,true);}
}

int main(int argc, char* argv[])
{
    bool type_double=true;
    std::string output_dir="output";

    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add("-o",&output_dir,"dir","Output directory");
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Parse(true);

    Create_Directory(output_dir+"/common");
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::printf("%s\n",parse_args.Print_Arguments());
    if(type_double){
        STREAM_TYPE stream_type((double()));
        Run<VECTOR<double,3> >(parse_args,stream_type,output_dir);}
    else{
        STREAM_TYPE stream_type((float()));
        Run<VECTOR<float,3> >(parse_args,stream_type,output_dir);}
    LOG::Finish_Logging();
    return 0;
}

