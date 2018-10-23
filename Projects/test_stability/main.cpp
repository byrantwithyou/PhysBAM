//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Arrays/PROJECTED_ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
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

template<class TV>
struct EXPLICIT_EXAMPLE
{
    typedef typename TV::SCALAR T;
    std::string out_dir;
    STREAM_TYPE stream_type;
    int last_frame,frame,step;
    T time,default_dt,frame_dt,dt;
    T alpha;
    bool test_matrix,dump_matrix;

    EXPLICIT_EXAMPLE():stream_type((T())),last_frame(100),default_dt(1e-3),frame_dt(1.0/24),alpha(1),
        test_matrix(false),dump_matrix(false)
    {}
};

template<class TV>
void Step(SOLID_BODY_COLLECTION<TV>& sbc,EXPLICIT_EXAMPLE<TV>& example)
{
    typedef typename TV::SCALAR T;
    sbc.Update_Position_Based_State(example.time,false,false);
    DEFORMABLE_PARTICLES<TV>& particles=sbc.deformable_body_collection.particles;
    ARRAY<TV> col(particles.number);
    auto idx=[](int p,int dim){return p*TV::m+dim;};
    int base=particles.number*TV::m;
    MATRIX_MXN<T> M;
    if(example.dump_matrix || example.test_matrix){
        MATRIX_MXN<T> A(particles.number*TV::m,particles.number*TV::m);
        M.Resize(particles.number*TV::m*2,particles.number*TV::m*2);
        ARRAY<TV> V(particles.number);
        for(int i=0;i<particles.number;i++){
            for(int j=0;j<TV::m;j++){
                V*=0;
                V(i)(j)=1;
                col*=0;
                sbc.deformable_body_collection.Add_Implicit_Velocity_Independent_Forces(V,col,example.time);
                for(int k=0;k<particles.number;k++){
                    for(int l=0;l<TV::m;l++){
                        T a=col(k)(l)/particles.mass(k);
                        M(idx(k,l),idx(i,j))=(idx(k,l)==idx(i,j)?1:0)+example.dt*example.dt*example.alpha*a;
                        M(base+idx(k,l),idx(i,j))=example.dt*a;
                        if(idx(k,l)==idx(i,j)){
                            M(idx(k,l),base+idx(i,j))=example.dt;
                            M(base+idx(k,l),base+idx(i,j))=1;}
                        A(idx(k,l),idx(i,j))=col(k)(l);}}}}
        OCTAVE_OUTPUT<T>(LOG::sprintf("%s/A%d.txt",example.out_dir,example.step).c_str()).Write("A",A);
        OCTAVE_OUTPUT<T>(LOG::sprintf("%s/M%d.txt",example.out_dir,example.step).c_str()).Write("M",M);
        ARRAY<std::complex<T> > eig(M.m);
        MATRIX_MXN<std::complex<T> > cM(M.m);
        for(int i=0;i<M.x.m;i++) cM.x(i)=(std::complex<T>)M.x(i);
        Compute_Eigenvalues(cM,eig);
        LOG::printf("eigenvalues: %P\n",eig);}

    ARRAY<T> U,U1;
    if(example.test_matrix){
        U.Resize(particles.number*TV::m*2);
        U1.Resize(particles.number*TV::m*2);
        for(int i=0;i<particles.number;i++){
            for(int j=0;j<TV::m;j++){
                U(idx(i,j))=particles.X(i)(j);
                U(base+idx(i,j))=particles.V(i)(j);}}
        U1=M*U;}

    col*=0;
    sbc.Add_Velocity_Independent_Forces(col,ARRAY_VIEW<TWIST<TV> >(),example.time);
    for(int i=0;i<particles.number;i++){
        particles.X(i)+=example.dt*example.dt*example.alpha/particles.mass(i)*col(i)+example.dt*particles.V(i);
        particles.V(i)+=example.dt/particles.mass(i)*col(i);}

    if(example.test_matrix){
        T x_error=-1,v_error=-1;
        for(int i=0;i<particles.number;i++){
            for(int j=0;j<TV::m;j++){
                x_error=max(x_error,abs(U1(idx(i,j))-particles.X(i)(j)));
                v_error=max(v_error,abs(U1(base+idx(i,j))-particles.V(i)(j)));
            }
        }
        LOG::printf("L-inf error x: %P v: %P\n",x_error,v_error);}
}

template<class TV>
void Setup(SOLID_BODY_COLLECTION<TV>& sbc,SOLIDS_STANDARD_TESTS<TV>& tests,int test_case)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    T density=1;
    RIGID_BODY_STATE<TV> initial_state;
    std::string data_directory="../../Public_Data";
    switch(test_case){
    case 0:{
            TETRAHEDRALIZED_VOLUME<T>& st=tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
                initial_state,true,true,density,1);
            // stiffness, poissons_ratio, damping
            FINITE_VOLUME<TV,TV::m>* force=Create_Finite_Volume(st,new COROTATED_FIXED<T,TV::m>(1e5,.45,0));
            sbc.Add_Force(force);
            sbc.deformable_body_collection.particles.X.template Project<T,&TV::y>()*=0.9;}
        break;
    case 1:{
            GRID<TV> mattress_grid(TV_INT()+2,RANGE<TV>(TV(),TV()+1),true);
            TETRAHEDRALIZED_VOLUME<T>& st=tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            // stiffness, poissons_ratio, damping
            FINITE_VOLUME<TV,TV::m>* force=Create_Finite_Volume(st,new COROTATED_FIXED<T,TV::m>(1e5,.45,0));
            sbc.Add_Force(force);
            sbc.deformable_body_collection.particles.X.template Project<T,&TV::y>()*=0.9;}
        break;
    default:PHYSBAM_FATAL_ERROR();}
}

template<class TV>
void Run(PARSE_ARGS& parse_args,STREAM_TYPE stream_type,const std::string& output_dir)
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    EXPLICIT_EXAMPLE<TV> example;
    example.out_dir=output_dir;
    int test_number;
    parse_args.Add("-test_matrix",&example.test_matrix,"Test matrix");
    parse_args.Add("-dump_matrix",&example.dump_matrix,"Dump transit matrix");
    parse_args.Add("-frame_dt",&example.frame_dt,"rate","Number of frames per second");
    parse_args.Add("-dt",&example.default_dt,"time","dt");
    parse_args.Add("-last_frame",&example.last_frame,"frame","number of frames to simulate");
    parse_args.Add("-a",&example.alpha,"alpha","alpha");
    parse_args.Extra(&test_number,"example number","example number to run");
    parse_args.Parse();

    std::string data_directory="../../Public_Data";
    VIEWER_OUTPUT<TV> vo(stream_type,GRID<TV>(),output_dir);
    SOLID_BODY_COLLECTION<TV> sbc;
    SOLIDS_STANDARD_TESTS<TV> tests(stream_type,data_directory,sbc);
    Setup(sbc,tests,test_number);
    Flush_Frame<TV>("init");
    sbc.Write(stream_type,output_dir,vo.frame-1,0,true,true,true,true,true);

    sbc.Update_Simulated_Particles();
    example.time=0;
    example.step=0;
    for(example.frame=0;example.frame<example.last_frame;example.frame++){
        T frame_end=(example.frame+1)*example.frame_dt;
        while(example.time<frame_end){
            example.dt=example.default_dt;
            T end_time=example.time+example.dt;
            if(end_time>=frame_end){
                example.dt=frame_end-example.time;
                end_time=frame_end;}
            Step(sbc,example);
            example.time=end_time;
            example.step++;}
        Flush_Frame<TV>("Frame");
        sbc.Write(stream_type,output_dir,vo.frame-1,0,true,true,true,true,true);}
}

int main(int argc, char* argv[])
{
    bool type_double=true;
    std::string output_dir="output";
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);

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

