#include "utils.h"

bool Compute_Barycentric(const VECTOR p1, const VECTOR p2, const VECTOR p3, 
                         const VECTOR p, double &u1, double &u2)
{
    //compute the barycentric of the point p in the triangle p1p2p3. 
    double base_det = (p1.x-p3.x)*(p2.y-p3.y) - (p1.y-p3.y)*(p2.x-p3.x);
    u1= ( (p.x-p3.x)*(p2.y-p3.y)-(p.y-p3.y)*(p2.x-p3.x) ) /base_det;
    u2= -( (p.x-p3.x)*(p1.y-p3.y)-(p.y-p3.y)*(p1.x-p3.x) ) /base_det;
    if( u1 < 0 && u1 > -EPSILON) u1=0;
    if( u2 < 0 && u2 > -EPSILON) u2=0;
    if( u1 > 1 && u1 < 1+EPSILON) u1=1;
    if( u2 > 1 && u2 < 1+EPSILON) u2=1;
    if( u1+u2 > 1 && u1+u2< 1+EPSILON) u1=1-u2;
    if(  u1>=0 && u1<=1 && u2>=0 && u2<=1 && (u1+u2)<=1)    return true;
    return false;
}

int SameTriples(int ti, int tj, int tk, int tti, int ttj, int ttk, int &notused)
{
    int count=0;
    notused=tti+ttj+ttk;
    if( tti==ti || tti==tj || tti==tk)
    {
        count++;
        notused-=tti;
    }
    if( ttj==ti || ttj==tj || ttj==tk)
    {
        count++;
        notused-=ttj;
    }
    if( ttk==ti || ttk==tj || ttk==tk)
    {
        count++;
        notused-=ttk;
    }
    return count;
}

