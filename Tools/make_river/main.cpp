#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FLOOD_FILL_2D.h>
#include <PhysBAM_Tools/Images/BMP_FILE.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include "../../Public_Library/Level_Sets/EXTRAPOLATION_2D.h"
#include "BEZIER_SPLINE.h"

using namespace PhysBAM; 

template<class T,class T2> void Make_River_From_Splines(GRID_2D<T>& grid,LEVELSET_2D<T>& levelset,ARRAYS<VECTOR<T,2> >& distances,ARRAYS<VECTOR<T,2> >& height,BEZIER_SPLINE<T,T2>& spline,BEZIER_SPLINE<T,T2>& bed_spline)
{
    ARRAYS<VECTOR<int,2> > colors(grid,1);ARRAYS<VECTOR<int,2> >::copy(0,colors);ARRAYS<VECTOR<int,2> >::put_ghost(-1,colors,grid,1);
    ARRAYS<VECTOR<VECTOR_2D<T> ,2> > tangents(grid,0);
    // rasterize the tangents, and arclength and the spline
    for(int segment=0;segment<spline.segments;segment++){
        ARRAY<T>& spline_distances=spline.arclength_table(segment);
        for(int k=0;k<spline.samples_per_segment;k++){
            T t=(k-1)*spline.one_over_samples_per_segment;
            VECTOR_2D<T> location=spline.f(segment,t);VECTOR_2D<T> tangent=spline.f_prime(segment,t);
            int i,j;grid.Closest_Node(location,i,j);
            tangents(i,j)=tangent;
            colors(i,j)=-1;distances(i,j)=spline_distances(k);}}
    ARRAYS<VECTOR<bool,2> > edge_is_blocked(grid,3);ARRAYS<VECTOR<bool,2> >::copy(false,edge_is_blocked);
    // flood fill to choose signs for each side of the river
    FLOOD_FILL_2D().Flood_Fill(colors,edge_is_blocked,edge_is_blocked);
    std::cout<<"colors.m="<<colors.m<<std::endl;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) if(colors(i,j)==1) levelset.phi(i,j)=grid.dx;else levelset.phi(i,j)=-grid.dx;

    // extrapolate tangents and arclength
    levelset.Fast_Marching_Method();
    EXTRAPOLATION_2D<T> extrapolation(grid,levelset.phi,distances);extrapolation.Set_Band_Width_To_Entire_Domain();
    EXTRAPOLATION_2D<T,VECTOR_2D<T> > extrapolation_tangent(grid,levelset.phi,tangents);extrapolation_tangent.Set_Band_Width_To_Entire_Domain();
    extrapolation.Extrapolate();extrapolation_tangent.Extrapolate();
    levelset.phi*=-1;extrapolation.Extrapolate();extrapolation_tangent.Extrapolate();levelset.phi*=-1;

    // Rasterize the profile spline for easy interpolation
    GRID<VECTOR<T,1> > profile_grid(100,0,1);ARRAYS<VECTOR<T,1> > profile(profile_grid,0);
    for(int segment=0;segment<bed_spline.segments;segment++){
        for(int k=0;k<bed_spline.samples_per_segment;k++){
            T t=(k-1)*bed_spline.one_over_samples_per_segment;
            VECTOR_2D<T> location=bed_spline.f(segment,t);
            int i;profile_grid.Clamped_Index(VECTOR_1D<T>(location.x/4),i);profile(i)=location.y;}}
    
    // choose height for spline.
    LINEAR_INTERPOLATION<T,T> profile_interpolation;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) height(i,j)=(T).4*(spline.total_length-distances(i,j))/spline.total_length+profile_interpolation.Clamped_To_Array(profile_grid,profile,VECTOR_1D<T>(fabs(levelset.phi(i,j))));

    // adjust based on arclength (make river slope downward)
    ARRAYS<VECTOR<T,2> > river_bed_height(grid,0);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) river_bed_height(i,j)=(T).4*(spline.total_length-distances(i,j))/spline.total_length+profile_interpolation.Clamped_To_Array(profile_grid,profile,VECTOR_1D<T>(0));
    
    FILE_UTILITIES::Write_To_File<T>("river_bed_height.gz",river_bed_height);
    FILE_UTILITIES::Write_To_File<T>("river_bed_grid.gz",grid);
    FILE_UTILITIES::Write_To_File<T>("river_bed_tangents.gz",tangents);
    FILE_UTILITIES::Write_To_File<T>("velocities.0.gz",tangents);
    FILE_UTILITIES::Write_To_File<T>("levelset.0.gz",levelset);
    FILE_UTILITIES::Write_To_File<T>("colors.0.gz",colors);
    FILE_UTILITIES::Write_To_File<T>("density.0.gz",distances);
    FILE_UTILITIES::Write_To_File<T>("temperature.0.gz",height);
    // write heightfield for photoshopping
    ARRAYS<VECTOR<VECTOR_3D<T> ,2> > image_out(grid,0);for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) image_out(i,j)=(T)height(i,j)*VECTOR_3D<T>(1,1,1);
    BMP_FILE().Write("image.bmp",image_out);
}

template<class T> void Make_River(){
    BEZIER_SPLINE<T,VECTOR_2D<T> > spline;
    spline.Add_Point(VECTOR_2D<T>(-0.008848,0.787145));
    spline.Add_Point(VECTOR_2D<T>(0.135392,0.784474));
    spline.Add_Point(VECTOR_2D<T>(0.322371,0.829883));
    spline.Add_Point(VECTOR_2D<T>(0.485309,0.883305));
    spline.Add_Point(VECTOR_2D<T>(0.485309,0.883305));
    spline.Add_Point(VECTOR_2D<T>(0.648247,0.936728));
    spline.Add_Point(VECTOR_2D<T>(0.720367,0.931386));
    spline.Add_Point(VECTOR_2D<T>(0.789816,0.851252));
    spline.Add_Point(VECTOR_2D<T>(0.789816,0.851252));
    spline.Add_Point(VECTOR_2D<T>(0.859265,0.771118));
    spline.Add_Point(VECTOR_2D<T>(0.811185,0.634891));
    spline.Add_Point(VECTOR_2D<T>(0.715025,0.594825));
    spline.Add_Point(VECTOR_2D<T>(0.715025,0.594825));
    spline.Add_Point(VECTOR_2D<T>(0.618865,0.554758));
    spline.Add_Point(VECTOR_2D<T>(0.330384,0.578798));
    spline.Add_Point(VECTOR_2D<T>(0.268948,0.479967));
    spline.Add_Point(VECTOR_2D<T>(0.268948,0.479967));
    spline.Add_Point(VECTOR_2D<T>(0.207512,0.381135));
    spline.Add_Point(VECTOR_2D<T>(0.250250,0.247579));
    spline.Add_Point(VECTOR_2D<T>(0.391820,0.210184));
    spline.Add_Point(VECTOR_2D<T>(0.391820,0.210184));
    spline.Add_Point(VECTOR_2D<T>(0.533389,0.172788));
    spline.Add_Point(VECTOR_2D<T>(0.682972,0.287646));
    spline.Add_Point(VECTOR_2D<T>(0.781803,0.290317));
    spline.Add_Point(VECTOR_2D<T>(0.781803,0.290317));
    spline.Add_Point(VECTOR_2D<T>(0.880634,0.292988));
    spline.Add_Point(VECTOR_2D<T>(1.030217,0.287646));
    spline.Add_Point(VECTOR_2D<T>(1.030217,0.307646));
    spline.Add_Point(VECTOR_2D<T>(1.030217,0.287646));
    spline.Add_Point(VECTOR_2D<T>(1.030217,0.287646));

    BEZIER_SPLINE<T,VECTOR_2D<T> > bed_spline;
    bed_spline.Add_Point(VECTOR_2D<T>(-0.019533,0.023205));
    bed_spline.Add_Point(VECTOR_2D<T>(0.060601,0.020534));
    bed_spline.Add_Point(VECTOR_2D<T>(0.097997,0.017863));
    bed_spline.Add_Point(VECTOR_2D<T>(0.151419,0.028548));
    bed_spline.Add_Point(VECTOR_2D<T>(0.151419,0.028548));
    bed_spline.Add_Point(VECTOR_2D<T>(0.204841,0.039232));
    bed_spline.Add_Point(VECTOR_2D<T>(0.218197,0.041903));
    bed_spline.Add_Point(VECTOR_2D<T>(0.250250,0.097997));
    bed_spline.Add_Point(VECTOR_2D<T>(0.250250,0.097997));
    bed_spline.Add_Point(VECTOR_2D<T>(0.282304,0.154090));
    bed_spline.Add_Point(VECTOR_2D<T>(0.295659,0.258264));
    bed_spline.Add_Point(VECTOR_2D<T>(0.309015,0.335726));
    bed_spline.Add_Point(VECTOR_2D<T>(0.309015,0.335726));
    bed_spline.Add_Point(VECTOR_2D<T>(0.322371,0.413189));
    bed_spline.Add_Point(VECTOR_2D<T>(0.389149,0.479967));
    bed_spline.Add_Point(VECTOR_2D<T>(0.469282,0.495993));
    bed_spline.Add_Point(VECTOR_2D<T>(0.469282,0.495993));
    bed_spline.Add_Point(VECTOR_2D<T>(0.549416,0.512020));
    bed_spline.Add_Point(VECTOR_2D<T>(0.592154,0.525375));
    bed_spline.Add_Point(VECTOR_2D<T>(0.656260,0.557429));
    bed_spline.Add_Point(VECTOR_2D<T>(0.656260,0.557429));
    bed_spline.Add_Point(VECTOR_2D<T>(0.720367,0.589482));
    bed_spline.Add_Point(VECTOR_2D<T>(0.781803,0.544073));
    bed_spline.Add_Point(VECTOR_2D<T>(0.829883,0.560100));
    bed_spline.Add_Point(VECTOR_2D<T>(0.829883,0.560100));
    bed_spline.Add_Point(VECTOR_2D<T>(0.877963,0.576127));
    bed_spline.Add_Point(VECTOR_2D<T>(1.046244,0.573456));
    bed_spline.Add_Point(VECTOR_2D<T>(1.110350,0.586811));
    bed_spline.Add_Point(VECTOR_2D<T>(1.110350,0.586811));
    bed_spline.Add_Point(VECTOR_2D<T>(1.174457,0.600167));

    GRID_2D<T> grid(201,201,0,1,0,1);ARRAYS<VECTOR<T,2> > phi(grid,0);ARRAYS<VECTOR<T,2> > distances(grid,0);ARRAYS<VECTOR<T,2> > height(grid,0);LEVELSET_2D<T> levelset(grid,phi);
    spline.Compute_Arclength();
    Make_River_From_Splines<T,VECTOR_2D<T> >(grid,levelset,distances,height,spline,bed_spline);
}

int main(int argc,char* argv[])
{
    Make_River<float>();
}

