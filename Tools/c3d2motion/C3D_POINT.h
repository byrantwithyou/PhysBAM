//##########################################################################################
// Class  C3D_POINT
// Records a single point
//##########################################################################################
#ifndef _C3D_POINT_h
#define _C3D_POINT_h

class C3D_POINT{
public:
    float x,y,z;
    C3D_BYTE camera_visible;
    float residual;
    void Print() const{printf("point (%g %g %g);\tcamera_visible %d;\tresidual %g\n",x,y,z,camera_visible,residual);}
};

#endif
