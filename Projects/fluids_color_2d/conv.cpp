#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>

using namespace PhysBAM;

template<class T>
T Approx_Exact(const ARRAY<int>& n,const ARRAY<T>& x,int order)
{
    ARRAY<int> w(x.m);
    MATRIX_MXN<T> A(x.m,order+1);
    for(int i=0;i<x.m;i++){
        T w=(T)1/n(i), p=1;
        for(int j=0;j<order+1;j++){
            A(i,j)=p;
            p*=w;}}

    ARRAY<T> coeffs=(A.Transpose_Times(A)).Solve_Linear_System(A.Transpose_Times(x));
    return coeffs(0);
}

int main(int argc,char *argv[])
{
    typedef double T;
    typedef double RW;
    STREAM_TYPE stream_type((RW()));
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;

    TV X0(TV()+FLT_MAX);
    int frame=1;
    ARRAY<std::string> sim_dirs;
    PARSE_ARGS parse_args(argc,argv);
    LOG::Initialize_Logging(false,false,1<<30,true);
    LOG::cout<<parse_args.Print_Arguments()<<std::endl;
    parse_args.Add("-frame",&frame,"frame","Frame to test");
    parse_args.Add("-x",&X0.x,"frame","Frame to test");
    parse_args.Add("-y",&X0.y,"frame","Frame to test");
    parse_args.Extra_Optional(&sim_dirs,"sim dirs","simulation directories");
    parse_args.Parse();

    T time;
    ARRAY<GRID<TV> > grids(sim_dirs.m);
    ARRAY<VECTOR<GRID<TV>,TV::m> > face_grids(sim_dirs.m);
    ARRAY<LEVELSET<TV>*> levelsets(sim_dirs.m);
    ARRAY<ARRAY<T,TV_INT> > pressures(sim_dirs.m),phis(sim_dirs.m);
    ARRAY<ARRAY<int,FACE_INDEX<TV::dimension> > > face_color(sim_dirs.m);
    ARRAY<ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > > > face_velocities(sim_dirs.m);
    ARRAY<int> res;

    for(int i=0;i<sim_dirs.m;i++){
        ARRAY<int,FACE_INDEX<TV::dimension> > prev_face_color;
        ARRAY<ARRAY<T,FACE_INDEX<TV::dimension> > > prev_face_velocities;
        FILE_UTILITIES::Read_From_File(stream_type,sim_dirs(i)+"/common/grid",grids(i));
        levelsets(i)=new LEVELSET<TV>(grids(i),phis(i));
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/levelset_%d",sim_dirs(i).c_str(),frame,0),*levelsets(i));
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/pressure",sim_dirs(i).c_str(),frame),pressures(i));
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/restart_data",sim_dirs(i).c_str(),frame),
            time,face_color(i),prev_face_color,face_velocities(i),prev_face_velocities);
        res.Append(grids(i).counts.x);
        for(int d=0;d<TV::m;d++) face_grids(i)(d)=grids(i).Get_Face_Grid(d);}

    ARRAY<int> samples(sim_dirs.m);
    ARRAY<T> u_linf(sim_dirs.m);
    ARRAY<T> u_2(sim_dirs.m);
    int index_min_res=res.Arg_Min();
    CUBIC_MN_INTERPOLATION_UNIFORM<GRID<TV>,T> interp;
    for(CELL_ITERATOR<TV> it(grids(index_min_res),-1);it.Valid();it.Next()){
        int col=phis(index_min_res)(it.index)>0;
        for(int d=0;d<TV::m;d++){
            ARRAY<T> values;
            for(int i=0;i<sim_dirs.m;i++){
                values.Append(interp.Periodic(face_grids(i)(d),face_velocities(i)(col).Component(d),it.Location()));}
            T exact=Approx_Exact(res,values,2);
            for(int i=0;i<sim_dirs.m;i++){
                T v=abs(values(i)-exact);
                samples(i)++;
                u_linf(i)=std::max(u_linf(i),v);
                u_2(i)+=sqr(v);}}}

    T Sx=0,Sy=0,Sz=0,Sxx=0,Sxy=0,Sxz=0;
    for(int i=0;i<sim_dirs.m;i++){
        T u2=sqrt(u_2(i)/samples(i)),uinf=u_linf(i),x=log10((T)res(i)),y=log10(u2),z=log10(uinf);
        Sx+=x;
        Sy+=y;
        Sz+=z;
        Sxx+=x*x;
        Sxy+=x*y;
        Sxz+=x*z;
        printf("%g %g\n", u2, uinf);}

    printf("order %g %g\n", -(Sy*Sx-sim_dirs.m*Sxy)/(-Sxx*sim_dirs.m+sqr(Sx)), -(Sz*Sx-sim_dirs.m*Sxz)/(-Sxx*sim_dirs.m+sqr(Sx)));

    if(X0.Max()<FLT_MAX){
        LOG::cout<<"sample at "<<X0<<std::endl;
        int col=0;
        for(int d=0;d<TV::m;d++){
            ARRAY<T> values;
            for(int i=0;i<sim_dirs.m;i++) values.Append(interp.Periodic(face_grids(i)(d),face_velocities(i)(col).Component(d),X0));
            LOG::cout<<"u("<<d<<")  "<<values<<"    "<<Approx_Exact(res,values,2)<<std::endl;}

        ARRAY<T> values,values2;
        for(int i=0;i<sim_dirs.m;i++) values.Append(interp.Periodic(grids(i),levelsets(i)->phi,X0));
        LOG::cout<<"level set  "<<values<<"    "<<Approx_Exact(res,values,2)<<std::endl;

        for(int i=0;i<sim_dirs.m;i++) values2.Append(interp.Periodic(grids(i),pressures(i),X0));
        LOG::cout<<"pressure  "<<values2<<"    "<<Approx_Exact(res,values2,2)<<std::endl;}
    LOG::Finish_Logging();
    return 0;
}
