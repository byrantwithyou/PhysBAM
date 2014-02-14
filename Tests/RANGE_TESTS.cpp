//#####################################################################
// Copyright 2008, Justin Solomon.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANGE_TESTS
//#####################################################################
// 1.  Tests statics
// 2.  Tests box arithmetic/boolean operators, basic operations
// 3.  Tests box "advanced" operators
// 4.  Tests functions in box .cpp file
//#####################################################################

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cyclic_shift.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Utilities/TEST_BASE.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry_Intersections/BOX_SPHERE_INTERSECTION.h>
#include <Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <limits>
namespace PhysBAM{

template<class T_input>
class RANGE_TESTS:public TEST_BASE
{
    typedef T_input T;
    mutable RANDOM_NUMBERS<T> rand;
    const static int num_iterations_per_test=1000;
    const static int num_iterations_per_subtest=100;
    T epsilon;

    bool Equal_To_Tolerance(T a,T b)
    {if(!(abs(a-b)<=epsilon*max(abs(a),abs(b),(T)1))){
        LOG::cout<<"failed on (a): "<<a<<std::endl;
        LOG::cout<<"failed on (b): "<<b<<std::endl;
        LOG::cout<<"failed on (c): "<<abs(a-b)<<std::endl;
        LOG::cout<<"failed on (d): "<<epsilon*max(abs(a),abs(b),(T)1)<<std::endl;
        LOG::cout<<"failed on (e): "<<(abs(a-b)<=epsilon*max(abs(a),abs(b),(T)1))<<std::endl;}
    return abs(a-b)<=epsilon*max(abs(a),abs(b),(T)1);}

    template<class TV>
    bool Equal_To_Tolerance(TV a,TV b)
    {return Equal_To_Tolerance((a-b).Magnitude(),(T)0);}

    template<class TV>
    RANGE<TV> Make_Random_Valid_Box()
    {TV random_edge1,random_edge2;
    rand.Fill_Uniform(random_edge1,-(T)1,(T)1); // bounds allow box to grow without overflow problems
    rand.Fill_Uniform(random_edge2,-(T)1,(T)1);
    return Construct_Box<TV>(random_edge1,random_edge2);}

    template<class TV>
    RANGE<TV> Construct_Box(const TV& vector1,const TV& vector2)
    {return RANGE<TV>(TV::Componentwise_Min(vector1,vector2),TV::Componentwise_Max(vector1,vector2));}

    bool Test(const bool test,const std::string& message,bool& ok) const
    {if(test) return true;LOG::cout<<"FAILED: "<<message<<std::endl;ok=false;return false;}

public:
    RANGE_TESTS()
        :TEST_BASE("box"),epsilon((T)20*std::numeric_limits<T>::epsilon())
    {rand.Set_Seed(17);}

    virtual ~RANGE_TESTS(){}

    TEST_RESULT Run_Test(int n)
    {switch(n){
            case 0:
            if(!Test_Static_Members<VECTOR<T,1> >() || !Test_Static_Members<VECTOR<T,2> >() || !Test_Static_Members<VECTOR<T,3> >()) return failure;
            return success;
        case 1:
            if(!Test_Box_Booleans<VECTOR<T,1> >() || !Test_Box_Booleans<VECTOR<T,2> >() || !Test_Box_Booleans<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Scaling<VECTOR<T,1> >() || !Test_Box_Scaling<VECTOR<T,2> >() || !Test_Box_Scaling<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Shift<VECTOR<T,1> >() || !Test_Box_Shift<VECTOR<T,2> >() || !Test_Box_Shift<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Volume<VECTOR<T,1> >() || !Test_Box_Volume<VECTOR<T,2> >() || !Test_Box_Volume<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Bounds_Adjustment<VECTOR<T,1> >() || !Test_Box_Bounds_Adjustment<VECTOR<T,2> >() || !Test_Box_Bounds_Adjustment<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Expansion<VECTOR<T,1> >() || !Test_Box_Expansion<VECTOR<T,2> >() || !Test_Box_Expansion<VECTOR<T,3> >()) return failure;
            return success;
        case 2:
            if(!Test_Box_Containment<VECTOR<T,1> >() || !Test_Box_Containment<VECTOR<T,2> >() || !Test_Box_Containment<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Projection<VECTOR<T,1> >() || !Test_Box_Projection<VECTOR<T,2> >() || !Test_Box_Projection<VECTOR<T,3> >()) return failure;
            if(!Test_Bounding_Box<VECTOR<T,1> >() || !Test_Bounding_Box<VECTOR<T,2> >() || !Test_Bounding_Box<VECTOR<T,3> >()) return failure;
            if(!Test_Remove_Dimension<VECTOR<T,2> >() || !Test_Remove_Dimension<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Intersection<VECTOR<T,1> >() || !Test_Box_Intersection<VECTOR<T,2> >() || !Test_Box_Intersection<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Project_Onto_Line<VECTOR<T,1> >() || !Test_Box_Project_Onto_Line<VECTOR<T,2> >() || !Test_Box_Project_Onto_Line<VECTOR<T,3> >()) return failure;
            return success;
        case 3:
            if(!Test_Box_Sphere_Intersection<VECTOR<T,2> >() || !Test_Box_Sphere_Intersection<VECTOR<T,3> >()) return failure;
            if(!Test_Box_Signed_Distance<VECTOR<T,1> >() || !Test_Box_Signed_Distance<VECTOR<T,2> >() || !Test_Box_Signed_Distance<VECTOR<T,3> >()) return failure;
            if(!Test_Ray_Box_Intersection<VECTOR<T,1> >() || !Test_Ray_Box_Intersection<VECTOR<T,2> >() || !Test_Ray_Box_Intersection<VECTOR<T,3> >()) return failure;
            return success;}
    return failure;}

    template<class TV>
    bool Test_Empty()
    {bool ok=true;
    Test(!RANGE<TV>::Unit_Box().Empty(),"Unit box is not empty.",ok);
    Test(!RANGE<TV>::Zero_Box().Empty(),"Zero box is not empty.",ok);
    Test(!RANGE<TV>::Full_Box().Empty(),"Full box is not empty.",ok);
    Test(RANGE<TV>::Empty_Box().Empty(),"Empty box is empty.",ok);
    return ok;}

    template<class TV>
    bool Test_Static_Members()
    {TV random_vector;bool ok=true;
    Test(RANGE<TV>::Unit_Box().Size()==1,"Unit box has exactly unit volume.",ok);
    Test(RANGE<TV>::Zero_Box().Lazy_Inside(TV()),"Zero box contains origin.",ok);
    for(int i=0;i<num_iterations_per_test && ok;i++){
        rand.Fill_Uniform(random_vector,(T)0,(T)1);
        Test(RANGE<TV>::Unit_Box().Lazy_Inside(random_vector),"Unit box contains random vector with components in [0,1).",ok);
        rand.Fill_Uniform(random_vector,-(T)1e10,(T)1e10);
        Test(!RANGE<TV>::Empty_Box().Lazy_Inside(random_vector),"Empty box does not contain a random vector.",ok);
        Test(RANGE<TV>::Full_Box().Lazy_Inside(random_vector),"Full box contains a random vector.",ok);
        Test(random_vector.Magnitude()<=epsilon || !RANGE<TV>::Zero_Box().Lazy_Inside(random_vector),"Zero box does not contain vectors outside the origin.",ok);}
    return ok;}

    template<class TV>
    bool Test_Box_Booleans()
    {VECTOR<TV,4> random_vectors;bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        for(int j=0;j<4;j++) rand.Fill_Uniform(random_vectors(j),-(T)1e10,(T)1e10);
        if(random_vectors(0)!=random_vectors(2) || random_vectors(1)!=random_vectors(3)){ // _highly_ unlikely this if would fail...
            RANGE<TV> box1(random_vectors(0),random_vectors(1)),box2(random_vectors(2),random_vectors(3));
            Test(box1!=box2,"Two disjoint boxes are not equal.",ok);}
        RANGE<TV> box(random_vectors(0),random_vectors(1)),box_copy(random_vectors(0),random_vectors(1));
        Test(box==box_copy,"Box equals its copy.",ok);}
    return ok;}

    template<class TV>
    bool Test_Box_Scaling()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>();
        box.max_corner+=(T).1; // Avoid sliver boxes, as these cannot computed accurately enough.
        T random_value=rand.Get_Uniform_Number((T).1,(T)10);
        Test(Equal_To_Tolerance((box*random_value).Size(),box.Size()*pow<TV::dimension,T>(random_value)),"Multiplication affects scaling of volume correctly.",ok);
        Test(Equal_To_Tolerance((box/random_value).Size(),box.Size()/pow<TV::dimension,T>(random_value)),"Division affects scaling of volume correctly.",ok);}
    return ok;}

    template<class TV>
    bool Test_Box_Shift()
    {TV test_vector,test_shift;bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>();
        T bound=max(box.Minimum_Corner().Magnitude(),box.Maximum_Corner().Magnitude())*2;
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            rand.Fill_Uniform(test_vector,-bound,bound);
            rand.Fill_Uniform(test_shift,-bound,bound);
            Test(box.Lazy_Inside(test_vector-test_shift)==(box+test_shift).Lazy_Inside(test_vector),"Shifting can be moved from box to vector in opposite directions.",ok);}}
    return ok;}

    // TODO: test Corners()/Surface_Area() functions -- dependent on dimension

    template<class TV>
    bool Test_Box_Volume() // also tested against scaling in Test_Box_Scaling() method and elsewhere
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> base_box=Make_Random_Valid_Box<TV>();
        TV shifted_min=base_box.Minimum_Corner(),shifted_max=base_box.Maximum_Corner();
        for(int j=0;j<4 && ok;j++){
            Test(Equal_To_Tolerance(RANGE<TV>(shifted_min,shifted_max).Size(),base_box.Size()),"Volume invariant to cyclic permutation.",ok);
            cyclic_shift(shifted_min);
            cyclic_shift(shifted_max);}}
    return ok;}

    template<class TV>
    bool Test_Box_Bounds_Adjustment()
    {TV random_vector,box_vector;bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>(),copy=box;
        rand.Fill_Uniform(random_vector,(T)-10000,(T)10000);
        box.Reset_Bounds(random_vector);
        Test(box.Lazy_Inside(random_vector),"Reset_bounds moves bounds correctly.",ok);
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            box_vector=rand.Get_Uniform_Vector(copy);
            Test(box_vector==random_vector || !box.Lazy_Inside(box_vector),"Reset_bounds does not include any points but center.",ok);}}

    ARRAY<TV> points(num_iterations_per_test);
    for(int i=0;i<num_iterations_per_subtest;i++)
        rand.Fill_Uniform(points(i),(T)-1000,(T)1000);
    for(int i=0;i<num_iterations_per_subtest && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>(),copy=box;
        // first test the array method
        box.Enlarge_Nonempty_Box_To_Include_Points(points);
        for(int j=0;j<num_iterations_per_test && ok;j++)
            Test(box.Lazy_Inside(points(j)),"Completely enlarged box contains all points for which it was enlarged.",ok);
        Test(box.Contains(copy),"Enlarged box contains original box.",ok);
        RANGE<TV> new_box=copy;
        for(int j=0;j<num_iterations_per_subtest && ok;j++){ // subtest here to save time...
            new_box.Enlarge_Nonempty_Box_To_Include_Point(points(j));
            Test(new_box.Lazy_Inside(points(j)),"Incrementally enlarged box contains all previous points.",ok);}
        Test(new_box.Contains(copy),"Incrementally enlarged box contains orignal box.",ok);
        Test(box==new_box,"Box enlargement approaches agree.",ok);}

    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box1=Make_Random_Valid_Box<TV>(),box2=Make_Random_Valid_Box<TV>(),copy0=box1,copy1=box2,combined=RANGE<TV>::Combine(box1,box2);
        copy0.Enlarge_To_Include_Box(box2);
        copy1.Enlarge_To_Include_Box(box1);
        Test(copy0.Contains(box1) && copy0.Contains(box2) && copy1.Contains(box1) && copy1.Contains(box2) && combined.Contains(box1) && combined.Contains(box2),
            "Simple box containment tests passed.",ok);}
    return ok;}

    template<class TV>
    bool Test_Box_Expansion()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){ // TODO: test Enlarge_By_Sign
        RANGE<TV> box=Make_Random_Valid_Box<TV>()/2,copy=box; // make box smaller so it's more likely that random points end on the border
        T random_delta=rand.Get_Uniform_Number((T)1e-3,(T)10);
        copy.Change_Size(random_delta);
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            TV random_inside_vector=rand.Get_Uniform_Vector(box),shift;
            rand.Fill_Uniform(shift,(T)1e-5,random_delta);
            Test(copy.Lazy_Inside(random_inside_vector+shift),"Scalar expansion allows for sufficiently-close vectors outside boundary.",ok);}
        copy=box;
        TV random_expansion;
        rand.Fill_Uniform(random_expansion,(T)1e-3,(T)4);
        copy.Change_Size(random_expansion);
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            TV random_inside_vector=rand.Get_Uniform_Vector(box),shift=rand.Get_Uniform_Vector(-random_expansion,random_expansion);
            Test(copy.Lazy_Inside(random_inside_vector+shift),"Vector expansion allows for sufficiently-close vectors outside boundary.",ok);}
        box.max_corner+=(T).1; // Make sure volumes can be computed accurately to pass tolerance.
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            copy=box;
            T random_delta_value=rand.Get_Uniform_Number((T).01,(T)10);
            copy.Change_Size(random_delta_value);
            Test(Equal_To_Tolerance(copy.Size(),box.Thickened(random_delta_value).Size()),"Change_Size and Thicken produce same volume.",ok);
            copy=box;
            copy.Scale_About_Center(random_delta_value);
            Test(Equal_To_Tolerance(copy.Size(),box.Size()*pow<TV::dimension,T>(random_delta_value)),"Scaling about the center affects volume correctly.",ok);
            TV random_delta_vector;
            rand.Fill_Uniform(random_delta_vector,(T)1e-3,(T)10);
            copy=box;
            copy.Scale_About_Center(random_delta_vector);
            Test(Equal_To_Tolerance(copy.Size(),box.Size()*random_delta_vector.Product()),"Scaling about the center by a vector affects volume correctly.",ok);}}
    return ok;}

    template<class TV>
    bool Test_Box_Containment()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){ // TODO: test boundary better
        RANGE<TV> box=Make_Random_Valid_Box<TV>();
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            TV random_vector=rand.Get_Uniform_Vector(box*(T)2);
            T random_value=rand.Get_Uniform_Number((T)0.01,(T)1);
            bool lazy_inside=box.Lazy_Inside(random_vector),lazy_outside=box.Lazy_Outside(random_vector),inside=box.Inside(random_vector,random_value),
                outside=box.Outside(random_vector,random_value),boundary=box.Boundary(random_vector,random_value);
            Test(!(lazy_inside==lazy_outside || (!lazy_inside && inside) || (!lazy_outside && outside) || (inside && outside) || (inside && !boundary && !lazy_inside)
                || (outside && !boundary && !lazy_outside) || (boundary && inside) || (boundary && outside)),"Compatibility tests for outside/inside booleans correct.",ok);
            random_vector=rand.Get_Uniform_Vector(box); // guaranteed to be inside this time
            Test(box.Lazy_Inside(random_vector),"Random vector generated to be inside box is indeed inside.",ok);}
        T random_delta=rand.Get_Uniform_Number((T)0.001,(T)1);
        TV min=box.Minimum_Corner(),max=box.Maximum_Corner();
        Test(box.Boundary(min,random_delta) && box.Boundary(max,random_delta) && box.Lazy_Inside(min) && !box.Lazy_Outside(min) && box.Lazy_Inside(max) && !box.Lazy_Outside(max),
            "Minimum/maximum corners pass inside/outside tests.",ok);
        Test(box.Lazy_Inside_Half_Open(min) && !box.Lazy_Inside_Half_Open(max),"Inside_Half_Open behaves correctly on corners.",ok);}
    return true;}

    template<class TV>
    bool Test_Box_Projection()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>();
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            TV random_vector=rand.Get_Uniform_Vector(box*(T)2);
            Test(!((box.Inside(random_vector,epsilon) && box.Clamp(random_vector)!=random_vector) || !box.Lazy_Inside(box.Clamp(random_vector))),
                "Clamp only does something if the vector is outside the box.",ok);
            rand.Fill_Uniform(random_vector,(T)-2,(T)2); // -2 to 2 so that sometimes it is in the unit box, sometimes not
            TV weighted_vector=box.Point_From_Normalized_Coordinates(random_vector);
            bool contained=box.Lazy_Inside(weighted_vector),in_unit_box=RANGE<TV>::Unit_Box().Inside(random_vector,epsilon),out_of_unit_box=RANGE<TV>::Unit_Box().Outside(random_vector,epsilon);
            Test(!((contained && out_of_unit_box) || (!contained && in_unit_box)),"Vector construction from weights places points in right locations.",ok);}}
    return ok;}

    template<class TV>
    bool Test_Box_Project_Onto_Line()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        TV center1=rand.template Get_Vector_In_Unit_Sphere<TV>(),center2=rand.template Get_Vector_In_Unit_Sphere<TV>();
        if((center1-center2).Magnitude()<(T)0.1){i--;continue;} // too close -- try again
        T distance=(center1-center2).Magnitude(),radius1=rand.Get_Uniform_Number((T)0.01,distance-(T)0.01),radius2=rand.Get_Uniform_Number((T)0.01,distance-radius1-(T)0.01);
        radius1/=(T)2;radius2/=(T)2; // because just because two points of a box are in the sphere does not mean the whole box is
        VECTOR<TV,4> random_vectors;
        for(int j=0;j<4;j++) random_vectors(j)=rand.template Get_Vector_In_Unit_Sphere<TV>();
        random_vectors(0)*=radius1;random_vectors(1)*=radius1,random_vectors(2)*=radius2;random_vectors(3)*=radius2;
        random_vectors(0)+=center1;random_vectors(1)+=center1;random_vectors(2)+=center2;random_vectors(3)+=center2;
        RANGE<TV> box1=Construct_Box<TV>(random_vectors(0),random_vectors(1)),box2=Construct_Box<TV>(random_vectors(2),random_vectors(3));
        VECTOR<T,2> region1,region2;
        box1.Project_Points_Onto_Line(center1-center2,region1(0),region1(1));
        box2.Project_Points_Onto_Line(center1-center2,region2(0),region2(1));
        Test(region1(0)<=region1(1) && region2(0)<=region2(1),"Regions from Test_Box_Project_Onto_Line have appropriate endpoints.",ok);
        Test(region1(1)<region2(0) || region2(1)<region1(0),"Test_Box_Project_Onto_Line regions from disjoint boxes do not overlap.",ok);}
    return ok;}

    template<class TV>
    bool Test_Box_Intersection()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box1=Make_Random_Valid_Box<TV>(),box2=Make_Random_Valid_Box<TV>(),combination=RANGE<TV>::Combine(box1,box2),intersection=RANGE<TV>::Intersect(box1,box2);
        T volume1=box1.Size(),volume2=box2.Size(),combination_volume=combination.Size(),intersection_volume=intersection.Size();
        Test(!(combination_volume<volume1 || combination_volume<volume2 || (!intersection.Empty() && (intersection_volume>volume1 || intersection_volume>volume2))),
            "\"Boolean\" volume identities hold",ok);
        Test(!(intersection.Empty() && (combination_volume<volume1+volume2 || intersection_volume>combination_volume)),"Volume sum/max properties hold under intersection and combination.",ok);
        TV min1=box1.Minimum_Corner(),min2=box2.Minimum_Corner(),max1=box1.Maximum_Corner(),max2=box2.Maximum_Corner();
        for(int j=0;j<TV::dimension && ok;j++){
            cyclic_shift(min1);cyclic_shift(min2);cyclic_shift(max1);cyclic_shift(max2);
            RANGE<TV> new_box1=RANGE<TV>(min1,max1),new_box2=RANGE<TV>(min2,max2);
            Test(Equal_To_Tolerance(intersection_volume,RANGE<TV>::Intersect(new_box1,new_box2).Size()) && Equal_To_Tolerance(combination_volume,RANGE<TV>::Combine(new_box1,new_box2).Size()),
                "Intersection volume preserved under cyclic shift.",ok);}

        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            TV random_vector=rand.Get_Uniform_Vector(RANGE<TV>::Unit_Box());
            bool in1=box1.Lazy_Inside(random_vector),in2=box2.Lazy_Inside(random_vector),in_intersection=intersection.Lazy_Inside(random_vector),
                in_combination=combination.Lazy_Inside(random_vector);
            Test(!in1 || !in2 || in_intersection,"If a random vector is in two boxes, it is in their intersection.",ok);
            Test((!in1 && !in2) || in_combination,"If a random vector is in a box, it is in that box combined with another.",ok);}
        bool intersection_empty=intersection.Empty(),intersection_1_2=!box1.Intersection(box2),intersection_2_1=!box2.Intersection(box1);
        Test(intersection_empty==intersection_1_2 && intersection_empty==intersection_2_1,"Intsersection size gives same empty/non-empty result as the boolean Intersection function.",ok);
        Test(intersection.Empty() || (Equal_To_Tolerance(box1.Intersection_Area(box2),intersection_volume) && Equal_To_Tolerance(box2.Intersection_Area(box1),intersection_volume)),
            "Intersection areas match.",ok);}
    return ok;}

    template<class TV>
    bool Test_Bounding_Box()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        ARRAY<TV> points(20);
        for(int j=0;j<20;j++) rand.Fill_Uniform(points(j),(T)-1000,(T)1000);
        RANGE<TV> box_two_points=RANGE<TV>::Bounding_Box(points(0),points(1)),box_three_points=RANGE<TV>::Bounding_Box(points(2),points(3),points(4)),
            box_four_points=RANGE<TV>::Bounding_Box(points(5),points(6),points(7),points(8)),box_full=RANGE<TV>::Bounding_Box(points);
        Test(box_two_points.Lazy_Inside(points(0)) && box_two_points.Lazy_Inside(points(1)),"Bounding box of two points contains those points.",ok);
        Test(box_three_points.Lazy_Inside(points(2)) && box_three_points.Lazy_Inside(points(3)) && box_three_points.Lazy_Inside(points(4)),"Bounding box of two points contains those points.",ok);
        Test(box_four_points.Lazy_Inside(points(5)) && box_four_points.Lazy_Inside(points(6)) && box_four_points.Lazy_Inside(points(7)) && box_four_points.Lazy_Inside(points(9)),
            "Bounding box of two points contains those points.",ok);
        for(int j=0;j<20 && ok;j++) Test(box_full.Lazy_Inside(points(j)),"Full bounding box contains all points from which it was constructed.",ok);}
    return ok;}

    template<class TV>
    bool Test_Remove_Dimension()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>();
        TV bounds=box.Edge_Lengths();
        T volume=box.Size();
        for(int j=0;j<TV::dimension && ok;j++) Test(Equal_To_Tolerance(volume/bounds(j),box.Remove_Dimension(j).Size()),"Volume changes according to removed dimension.",ok);
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            TV inside_vector=rand.Get_Uniform_Vector(box);
            Test(box.Lazy_Inside(inside_vector),"Box contains random vector inside of it.",ok);
            for(int k=0;k<TV::dimension && ok;k++) Test(box.Remove_Dimension(k).Lazy_Inside(inside_vector.Remove_Index(k)),"Box with removed dimension contains vector with removed dimension.",ok);}}
    return ok;}

    template<class TV>
    void Test_Ray_Box_Intersection_3D(const RANGE<TV>& box,const TV& direction,const RAY<TV>& ray0,const RAY<TV>& ray1,const RAY<TV>& ray2,RAY<TV>& ray3,bool& ok)
    {STATIC_ASSERT(TV::m!=3);}

    template<class T>
    void Test_Ray_Box_Intersection_3D(const RANGE<VECTOR<T,3> >& box,const VECTOR<T,3>& direction,RAY<VECTOR<T,3> >& ray0,RAY<VECTOR<T,3> >& ray1,
        RAY<VECTOR<T,3> >& ray2,RAY<VECTOR<T,3> >& ray3,bool& ok)
    {Test(INTERSECTION::Lazy_Intersects(ray0,box,-epsilon) && INTERSECTION::Lazy_Intersects(ray1,box,-epsilon),"Box intersects rays with one endpoint inside; uses \"Lazy Inside\" routine.",ok);
    if(direction.Magnitude()>epsilon) Test(!INTERSECTION::Lazy_Intersects(ray2,box,epsilon),"Ray outside radius of box does not intersect box; uses \"Lazy Inside\" routine.",ok);
    Test(!INTERSECTION::Lazy_Intersects(ray3,box,epsilon),"Radial ray does not intersect box; uses \"Lazy Inside\" routine.",ok);
    Test(!INTERSECTION::Lazy_Outside(ray0,box) && !INTERSECTION::Lazy_Outside(ray1,box),"Lazy_Outside gives complement of inside value.",ok);} // TODO: Test this function more

    template<class TV>
    bool Test_Ray_Box_Intersection()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>();
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            TV random_inside1=rand.Get_Uniform_Vector(box),random_inside2=rand.Get_Uniform_Vector(box),random_vector;
            rand.Fill_Uniform(random_vector,(T)-2,(T)2);
            RAY<TV> ray0(random_vector,random_inside1-random_vector),ray1(random_inside2,random_inside1-random_inside2),ray1_copy=ray0,ray2_copy=ray1; // both of these intersect the box
            Test(INTERSECTION::Intersects(ray0,box,epsilon) && INTERSECTION::Intersects(ray1,box,epsilon),"Box intersects rays with one endpoint inside.",ok);
            T radius=(box.Maximum_Corner()-box.Minimum_Corner()).Magnitude()/(T)2;
            TV endpoint=box.Center()+rand.template Get_Direction<TV>()*rand.Get_Uniform_Number(radius+(T).1,radius+(T)3),direction=rand.template Get_Direction<TV>();
            direction-=direction.Projected(endpoint-box.Center());
            RAY<TV> ray2(endpoint,direction),ray3(endpoint,endpoint-box.Center()),ray3_copy=ray2,ray4_copy=ray3;
            if(direction.Magnitude()>epsilon) Test(!INTERSECTION::Intersects(ray2,box,epsilon),"Ray outside radius of box does not intersect box.",ok);
            Test(!INTERSECTION::Intersects(ray3,box,epsilon),"Radial ray does not intersect box.",ok);
            Test_Ray_Box_Intersection_3D(box,direction,ray1_copy,ray2_copy,ray3_copy,ray4_copy,ok);}}
    return ok;}

    template<class TV>
    bool Test_Box_Sphere_Intersection()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>();
        T box_radius=(box.Maximum_Corner()-box.Minimum_Corner()).Magnitude();
        TV box_center=box.Center();
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            T distance_from_box=rand.Get_Uniform_Number(box_radius+(T)100*epsilon,(T)10);
            TV random_center=box_center+rand.template Get_Direction<TV>()*distance_from_box;
            T random_radius=rand.Get_Uniform_Number(epsilon,distance_from_box-box_radius-(T)10*epsilon);
            SPHERE<TV> sphere(random_center,random_radius);
            Test(!INTERSECTION::Intersects(box,sphere,epsilon),"Box and far-away sphere do not intersect.",ok);
            TV intersecting_center;
            rand.Fill_Uniform(intersecting_center,(T)-2,(T)2);
            T intersecting_radius=(intersecting_center-rand.Get_Uniform_Vector(box)).Magnitude();
            SPHERE<TV> intersecting_sphere(intersecting_center,intersecting_radius);
            Test(INTERSECTION::Intersects(box,intersecting_sphere,epsilon),"Box intersects with sphere containing at least one box point.",ok);
            TV random_outside_center;
            rand.Fill_Uniform(random_outside_center,(T)-2,(T)2);
            while(box.Inside(random_outside_center,(T)2*epsilon)) rand.Fill_Uniform(random_outside_center,(T)-2,(T)2);
            SPHERE<TV> outside_sphere(random_outside_center,(random_outside_center-box.Clamp(random_outside_center)).Magnitude()-(T)2*epsilon);
            Test(!INTERSECTION::Intersects(box,outside_sphere,epsilon),"Box does not intersect with sphere constructed to be on the outside.",ok);
            TV random_inside_center=rand.Get_Uniform_Vector(box); // also tests "Surface" function
            SPHERE<TV> inside_sphere(random_inside_center,(random_inside_center-box.Surface(random_inside_center)).Magnitude()-(T)2*epsilon);
            Test(INTERSECTION::Intersects(box,inside_sphere,epsilon),"Box intersects sphere completely on the inside.",ok);}}
    return ok;}

    template<class TV>
    bool Test_Box_Signed_Distance()
    {bool ok=true;
    for(int i=0;i<num_iterations_per_test && ok;i++){
        RANGE<TV> box=Make_Random_Valid_Box<TV>();
        for(int j=0;j<num_iterations_per_subtest && ok;j++){
            TV random_inside_vector=rand.Get_Uniform_Vector(box);
            TV random_outside_vector;
            rand.Fill_Uniform(random_outside_vector,(T)-2,(T)2);
            while(box.Inside(random_outside_vector,(T)2*epsilon)) rand.Fill_Uniform(random_outside_vector,(T)-2,(T)2);
            Test(box.Signed_Distance(random_inside_vector)<=(T)0,"Interior point has negative signed distance from box.",ok);
            Test(box.Signed_Distance(random_outside_vector)>=(T)0,"Exterior point has positive signed distance from box.",ok);
            Test(box.Signed_Distance(random_inside_vector)+box.Signed_Distance(random_outside_vector)<(random_inside_vector-random_outside_vector).Magnitude(),
                "Signed distances should violate the triangle inequality for inside/outside pair.",ok);
            TV inside_projection=box.Surface(random_inside_vector),outside_projection=box.Surface(random_outside_vector);
            Test(Equal_To_Tolerance(box.Signed_Distance(random_outside_vector),(outside_projection-random_outside_vector).Magnitude()),
                "Signed distance for outside point is same as distance to projection onto surface of box.",ok);
            Test(Equal_To_Tolerance(box.Signed_Distance(random_inside_vector),-(inside_projection-random_inside_vector).Magnitude()),
                "Signed distance for inside point is same as distance to projection onto surface of box times -1.",ok);}}
    return ok;}

    int Number_Of_Tests() const
    {return 4;}
//#####################################################################
};
static RANGE_TESTS<double> box_tests;
}
