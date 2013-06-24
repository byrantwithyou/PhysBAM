#include <Tools/Arrays/ARRAY.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/cube.h>

namespace PhysBAM{
template<class T,class T2>
class BEZIER_SPLINE{
public:
    int segments;
    ARRAY<T2> control_points;
    ARRAY<ARRAY<T> > arclength_table;
    int samples_per_segment;
    T one_over_samples_per_segment;
    T total_length;

    BEZIER_SPLINE()
    {
        samples_per_segment=1000;
        one_over_samples_per_segment=(T)1/(T)(samples_per_segment-1);
    }

    T2 f(int segment,T t)
    {int base=(segment-1)*4+1;T2 c1=control_points(base),c2=control_points(base+1),c3=control_points(base+2),c4=control_points(base+3);
    return cube(1-t)*c1+3*t*sqr(1-t)*c2+3*sqr(t)*(1-t)*c3+cube(t)*c4;}

    T2 f_prime(int segment,T t)
    {int base=(segment-1)*4+1;T2 c1=control_points(base),c2=control_points(base+1),c3=control_points(base+2),c4=control_points(base+3);
    return -(T)(3*sqr(1-t))*c1+(T)(3*sqr(1-t)-6*(1-t)*t)*c2+(T)(6*(1-t)*t-3*sqr(t))*c3+(T)(3*sqr(t))*c4;}

    T2 f_arc(T s)
    {int segment;T t;Find_Segment(s,segment,t);return f(segment,t);}

    void Find_Segment(const T s,int& segment,T& t)
    {assert(s>=0&&s<=total_length);
    for(int segment=0;segment<segments;segment++) if(arclength_table(segment)(samples_per_segment) >= s) 
        for(int i=0;i<samples_per_segment;i++) if(arclength_table(segment)(i) >= s){segment=segment;t=(i-1)*one_over_samples_per_segment;return;}
    assert(false);}

    void Compute_Arclength()
    {arclength_table.Resize(segments);
    T distance=0;
    for(int segment=0;segment<segments;segment++){
        arclength_table(segment).Resize(samples_per_segment);
        T2 old_value=f(segment,0);
        for(int i=0;i<samples_per_segment;i++){
            T t=(i-1)*one_over_samples_per_segment;T2 new_value=f(segment,t);
            distance+=(new_value-old_value).Magnitude();arclength_table(segment)(i)=distance;
            //std::cout<<"Distance "<<segment<<" i="<<i<<" t="<<t<<" "<<distance<<std::endl;
            old_value=new_value;}}
    total_length=distance;}
    
    int Add_Point(const T2& point)
    {control_points.Append(point);segments=control_points.m/4;return control_points.m;}
};
}
