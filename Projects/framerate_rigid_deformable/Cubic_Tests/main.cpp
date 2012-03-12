#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <iomanip>

using namespace PhysBAM;

void Check_Intervals(CUBIC<double>& cubic,int num_intervals,VECTOR<INTERVAL<double>,3> intervals,int mask,int case_index,VECTOR<double,3> roots,double a,double t,double u) {
    if(!(((1<<num_intervals) & mask) || (num_intervals==0 && mask==0))){
        LOG::cout << "error: the number of intervals is wrong " << case_index << std::endl;
        LOG::cout << "cubic " << cubic.c3 << " " << cubic.c2 << " " << cubic.c1 << " " << cubic.c0 << std::endl;
        LOG::cout << "cubic " << cubic(roots(1)) << " "  << cubic(roots(2)) << " "  << cubic(roots(3)) << std::endl;
        LOG::cout << "a " << a << std::endl;
        LOG::cout << "cubic " << roots(1) << " "  << roots(2) << " " << roots(3) << std::endl;
        LOG::cout << "num_intervals " << num_intervals << std::endl;
        LOG::cout << "intervals " << intervals << std::endl;
        LOG::cout << "passed interval " << t << " " << u << std::endl;
        LOG::cout << "mask " << mask << std::endl;
        return;
    }
    for(int i=0;i<num_intervals;i++){
        double interval_endpoint1=cubic(intervals(i).min_corner);
        double interval_endpoint2=cubic(intervals(i).max_corner);
        if(interval_endpoint1*interval_endpoint2>0 && (abs(interval_endpoint1) > 1e-11 && abs(interval_endpoint2) > 1e-11)){
            LOG::cout << "case "<<case_index<<": didn't find root " << i << std::endl;
            LOG::cout << interval_endpoint1 << " " << interval_endpoint2 << std::endl;
            LOG::cout << "cubic " << cubic.c3 << " " << cubic.c2 << " " << cubic.c1 << " " << cubic.c0 << std::endl;
            LOG::cout << "cubic " << cubic(roots(1)) << " "  << cubic(roots(2)) << " "  << cubic(roots(3)) << std::endl;
            LOG::cout << "a " << a << std::endl;
            LOG::cout << "cubic " << roots(1) << " "  << roots(2) << " " << roots(3) << std::endl;
            LOG::cout << "num_intervals " << num_intervals << std::endl;
            LOG::cout << "intervals " << intervals << std::endl;
            LOG::cout << "passed interval " << t << " " << u << std::endl;
            LOG::cout << "mask " << mask << std::endl;
        }}
}

int main(int argc,char* argv[])
{
    LOG::cout << std::setprecision(16);
/*    CUBIC<double> cubic(atof(argv[1]),atof(argv[2]),atof(argv[3]),atof(argv[4]));
    double r1=atof(argv[5]),r2=atof(argv[6]),r3=atof(argv[7]);
    INTERVAL<double> blah(atof(argv[8]),atof(argv[9]));
    if(r1!=100) blah.Enlarge_To_Include_Point(r1);
    if(r2!=100) blah.Enlarge_To_Include_Point(r2);
    if(r3!=100) blah.Enlarge_To_Include_Point(r3);
    blah*=1.5;
    blah=blah.Thickened(.1);
    RANDOM_NUMBERS<double> rn;
    for(int i=0;i<1000000;i++){
        double t=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
        double u=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
        if(t>u) exchange(t,u);
        int intervals=0;
        INTERVAL<double> interval1,interval2,interval3;
        cubic.Compute_Intervals(t,u,intervals,interval1,interval2,interval3);
        
        int total=0;
        if(INTERVAL<double>(t,u).Lazy_Inside(r1)) total++;
        if(INTERVAL<double>(t,u).Lazy_Inside(r2)) total++;
        if(INTERVAL<double>(t,u).Lazy_Inside(r3)) total++;
        if(total!=intervals) LOG::cout << "error: " << t << " " << u << std::endl;
    }
*/
    RANDOM_NUMBERS<double> rn;
    for(int i=0;i<1000000000;i++){
        int index=rn.Get_Uniform_Integer(1,4);
        VECTOR<double,3> roots;
        roots(1)=rn.Get_Uniform_Number(-10,10);
        roots(2)=rn.Get_Uniform_Number(-10,10);
        roots(3)=rn.Get_Uniform_Number(-10,10);
        double a=rn.Get_Uniform_Number(-10,10);

        switch(index){
            case 1:
                {// a(x-r)(x-s)(x-t), roots r,s,t
                CUBIC<double> cubic(a,-a*(roots(1)+roots(2)+roots(3)),a*(roots(1)*roots(2)+roots(1)*roots(3)+roots(2)*roots(3)),-a*roots(1)*roots(2)*roots(3));
                /*LOG::cout << "cubic " << cubic.c3 << " " << cubic.c2 << " " << cubic.c1 << " " << cubic.c0 << std::endl;
                  LOG::cout << "cubic " << cubic(roots(1)) << " "  << cubic(roots(2)) << " "  << cubic(roots(3)) << std::endl;*/
                INTERVAL<double> blah=INTERVAL<double>::Empty_Box();
                blah.Enlarge_To_Include_Point(roots(1));
                blah.Enlarge_To_Include_Point(roots(2));
                blah.Enlarge_To_Include_Point(roots(3));
                blah=blah.Thickened(5);

                double t=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
                double u=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
                if(t>u) exchange(t,u);
                INTERVAL<double> new_interval(t,u);
                int mask=0;
                int normal_roots=0;
                for(int i=0;i<3;i++) if(new_interval.Lazy_Inside(roots(i))) normal_roots++;
                mask|=(1<<normal_roots);
                int shrunk_roots=0;
                INTERVAL<double> shrunk_interval=new_interval.Thickened(1e-11);
                for(int i=0;i<3;i++) if(shrunk_interval.Lazy_Inside(roots(i))) shrunk_roots++;
                mask|=(1<<shrunk_roots);
                int expanded_roots=0;
                INTERVAL<double> expanded_interval=new_interval.Thickened(-1e-11);
                for(int i=0;i<3;i++) if(expanded_interval.Lazy_Inside(roots(i))) expanded_roots++;
                mask|=(1<<expanded_roots);
                int num_intervals=0;
                VECTOR<INTERVAL<double>,3> intervals;
                cubic.Compute_Intervals(t,u,num_intervals,intervals(1),intervals(2),intervals(3));
                Check_Intervals(cubic,num_intervals,intervals,mask,1,roots,a,t,u);
/*                LOG::cout << "case 1: a " << a << " r("<<roots(1)<<"), s("<<roots(2)<<"), t("<<roots(3)<<") " << std::endl;
                LOG::cout << "num_intervals " << num_intervals << std::endl;
                LOG::cout << "intervals " << intervals << std::endl;*/
                break;}
            case 2:
                {// a(x-r)^2(x-s) roots, r, s
                CUBIC<double> cubic(a,-a*(2*roots(1)+roots(2)),a*(roots(1)*roots(1)+2*roots(1)*roots(2)),-a*roots(1)*roots(1)*roots(2));
/*                LOG::cout << "a " << a << std::endl;
                  LOG::cout << "cubic " << roots(1) << " "  << roots(2) << std::endl;*/
                INTERVAL<double> blah=INTERVAL<double>::Empty_Box();
                blah.Enlarge_To_Include_Point(roots(1));
                blah.Enlarge_To_Include_Point(roots(2));
                blah=blah.Thickened(5);

                double t=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
                double u=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
                if(t>u) exchange(t,u);
                INTERVAL<double> new_interval(t,u);
                int mask=0;
                if(new_interval.Lazy_Inside(roots(1))){
                    if(new_interval.Lazy_Inside(roots(2))) mask|=10;
                    else mask|=5;}
                else if(new_interval.Lazy_Inside(roots(2))) mask|=2;
                INTERVAL<double> shrunk_interval=new_interval.Thickened(1e-11);
                if(shrunk_interval.Lazy_Inside(roots(1))){
                    if(shrunk_interval.Lazy_Inside(roots(2))) mask|=10;
                    else mask|=5;}
                else if(shrunk_interval.Lazy_Inside(roots(2))) mask|=2;
                INTERVAL<double> expanded_interval=new_interval.Thickened(-1e-11);
                if(expanded_interval.Lazy_Inside(roots(1))){
                    if(expanded_interval.Lazy_Inside(roots(2))) mask|=10;
                    else mask|=5;}
                else if(expanded_interval.Lazy_Inside(roots(2))) mask|=2;

                int num_intervals=0;
                VECTOR<INTERVAL<double>,3> intervals;
                cubic.Compute_Intervals(t,u,num_intervals,intervals(1),intervals(2),intervals(3));
                Check_Intervals(cubic,num_intervals,intervals,mask,2,roots,a,t,u);
                break;}
            case 3:
                {// a(x-r)^3, roots r
//              LOG::cout << "cubic roots " << roots << std::endl;
                CUBIC<double> cubic(a,-3*a*roots(1),a*3*roots(1)*roots(1),-a*roots(1)*roots(1)*roots(1));
//              LOG::cout << "cubic at roots " << cubic(roots(1)) << " "  << cubic(roots(2)) << " "  << cubic(roots(3)) << std::endl;
                INTERVAL<double> blah=INTERVAL<double>::Empty_Box();
                blah.Enlarge_To_Include_Point(roots(1));
                blah=blah.Thickened(5);

                double t=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
                double u=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
                if(t>u) exchange(t,u);
                INTERVAL<double> new_interval(t,u);
                int mask=0;
                if(new_interval.Lazy_Inside(roots(1))) mask|=10;
                INTERVAL<double> shrunk_interval=new_interval.Thickened(1e-11);
                if(shrunk_interval.Lazy_Inside(roots(1))) mask|=10;
                INTERVAL<double> expanded_interval=new_interval.Thickened(-1e-11);
                if(expanded_interval.Lazy_Inside(roots(1))) mask|=10;

                int num_intervals=0;
                VECTOR<INTERVAL<double>,3> intervals;
                cubic.Compute_Intervals(t,u,num_intervals,intervals(1),intervals(2),intervals(3));
                Check_Intervals(cubic,num_intervals,intervals,mask,3,roots,a,t,u);
                break;}
            case 4:
                {// a(x-r)((x-s)^2+t^2) t!=0, roots r
                if(roots(3)==0) continue;
                CUBIC<double> cubic(a,-a*(2*roots(2)+roots(1)),a*(roots(2)*roots(2)+roots(3)*roots(3)+2*roots(2)*roots(1)),-a*(roots(1)*roots(2)*roots(2)+roots(1)*roots(3)*roots(3)));
                INTERVAL<double> blah=INTERVAL<double>::Empty_Box();
                blah.Enlarge_To_Include_Point(roots(1));
                blah=blah.Thickened(5);

                double t=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
                double u=rn.Get_Uniform_Number(blah.min_corner,blah.max_corner);
                if(t>u) exchange(t,u);
                INTERVAL<double> new_interval(t,u);
                int expected_number=0;
                if(new_interval.Lazy_Inside(roots(1))) expected_number++;

                int num_intervals=0;
                VECTOR<INTERVAL<double>,3> intervals;
                cubic.Compute_Intervals(t,u,num_intervals,intervals(1),intervals(2),intervals(3));
                if(num_intervals!=expected_number) LOG::cout << "case 4: num intervals isn't expected" << std::endl;
                if(num_intervals>0 && !intervals(1).Lazy_Inside(roots(1))){
                    LOG::cout << "case 4: didn't find root 1" << std::endl;
                    double interval_endpoint1=cubic(intervals(1).min_corner);
                    double interval_endpoint2=cubic(intervals(1).max_corner);
                    LOG::cout << interval_endpoint1 << " " << interval_endpoint2 << std::endl;
                    LOG::cout << "cubic " << cubic.c3 << " " << cubic.c2 << " " << cubic.c1 << " " << cubic.c0 << std::endl;
                    LOG::cout << "cubic " << cubic(roots(1)) << " "  << cubic(roots(2)) << " "  << cubic(roots(3)) << std::endl;
                    LOG::cout << "a " << a << std::endl;
                    LOG::cout << "cubic " << roots(1) << " "  << roots(2) << " " << roots(3) << std::endl;
                    LOG::cout << "num_intervals " << num_intervals << std::endl;
                    LOG::cout << "intervals " << intervals << std::endl;
                    LOG::cout << "passed interval " << t << " " << u << std::endl;}
                break;}
        }
    }
}
