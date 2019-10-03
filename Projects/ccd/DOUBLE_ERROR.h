#ifndef __DOUBLE_ERROR__
#define __DOUBLE_ERROR__

#include <Core/Log/LOG.h>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include "quadmath.h"
#include <type_traits>

namespace PhysBAM {
typedef __float128 quad;

std::string qdtostr(__float128 q)
{
    int width=46;
    char buf[128];
    quadmath_snprintf (buf, sizeof buf, "%0.20Qe", width, q);
    return std::string(buf);
}

struct DOUBLE_ERROR
{
    std::string val;    // String representation of all operations leading to the derived values of dbl and qd
    double dbl;         // Double precision (64 bits) representation of DOUBLE_ERROR
    quad  qd;           // Quad precision (128 bits) representation of DOUBLE_ERROR

    static const constexpr double EPSILON = DBL_EPSILON;
    DOUBLE_ERROR():val("0"),dbl(0),qd(0){}

    template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    DOUBLE_ERROR(T t):val(""),dbl(t),qd(t)
    {
        std::ostringstream out;
        if(t>=0)out<<" ";
        out<<std::setprecision(20)<<std::fixed<<t;
        val=out.str();
    }
    DOUBLE_ERROR(std::string a,double b,quad c):val(a),dbl(b),qd(c){}
    DOUBLE_ERROR(const DOUBLE_ERROR& d):val(d.val),dbl(d.dbl),qd(d.qd){}

    DOUBLE_ERROR operator+(const DOUBLE_ERROR& d) const{return DOUBLE_ERROR("ADD("+val+","+d.val+")",dbl+d.dbl,qd+d.qd);}
    DOUBLE_ERROR operator-(const DOUBLE_ERROR& d) const{return DOUBLE_ERROR("SUB("+val+","+d.val+")",dbl-d.dbl,qd-d.qd);}
    DOUBLE_ERROR operator*(const DOUBLE_ERROR& d) const{return DOUBLE_ERROR("MUL("+val+","+d.val+")",dbl*d.dbl,qd*d.qd);}
    DOUBLE_ERROR operator/(const DOUBLE_ERROR& d) const{return DOUBLE_ERROR("DIV("+val+","+d.val+")",dbl/d.dbl,qd/d.qd);}

    DOUBLE_ERROR& operator+=(const DOUBLE_ERROR& d){return *this=*this+d;}
    DOUBLE_ERROR& operator-=(const DOUBLE_ERROR& d){return *this=*this-d;}
    DOUBLE_ERROR& operator*=(const DOUBLE_ERROR& d){return *this=*this*d;}
    DOUBLE_ERROR& operator/=(const DOUBLE_ERROR& d){return *this=*this/d;}

    DOUBLE_ERROR operator-() const {return DOUBLE_ERROR("NEG("+val+")",-1*dbl,-1*qd);}

    bool operator<(const DOUBLE_ERROR& d) const{if((dbl<d.dbl)!=(qd<d.qd)) LOG::printf("Conditional Disagreement:< \n\t%P\n\t%P\n",val,d.val);return qd<d.qd;}
    bool operator>(const DOUBLE_ERROR& d) const{if((dbl>d.dbl)!=(qd>d.qd)) LOG::printf("Conditional Disagreement:> \n\t%P\n\t%P\n",val,d.val);return qd>d.qd;}
    bool operator<=(const DOUBLE_ERROR& d) const{if((dbl<=d.dbl)!=(qd<=d.qd)) LOG::printf("Conditional Disagreement:<=\n\t%P\n\t%P\n",val,d.val);return qd<=d.qd;}
    bool operator>=(const DOUBLE_ERROR& d) const{if((dbl>=d.dbl)!=(qd>=d.qd)) LOG::printf("Conditional Disagreement:>=\n\t%P\n\t%P\n",val,d.val);return qd>=d.qd;}
    bool operator==(const DOUBLE_ERROR& d) const{if((dbl==d.dbl)!=(qd==d.qd)) LOG::printf("Conditional Disagreement:==\n\t%P\n\t%P\n",val,d.val);return qd==d.qd;}
    bool operator!=(const DOUBLE_ERROR& d) const{if((dbl!=d.dbl)!=(qd!=d.qd)) LOG::printf("Conditional Disagreement:!=\n\t%P\n\t%P\n",val,d.val);return qd!=d.qd;}

    quad error() const {return fabsq(qd-(quad)dbl);};

    operator bool(){return dbl!=0;}  // Needed for VECTOR3D.h:238...? //TODO: Why...
};

std::ostream& operator<<(std::ostream& os,const DOUBLE_ERROR& d)
{
    os<<std::scientific<<std::setprecision(21)<<d.dbl;
    return os;
}

// Binary Operator:std::is_arithmetic<T> and DOUBLE_ERROR (+,-,*,/,<,>,<=,>=,==,!=)
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
DOUBLE_ERROR operator+(const DOUBLE_ERROR& d,const T& t){return DOUBLE_ERROR("ADD("+d.val+","+std::to_string(t)+")",d.dbl+t,d.qd+t);}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
DOUBLE_ERROR operator+(const T& t,const DOUBLE_ERROR& d){return d+t;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
DOUBLE_ERROR operator-(const DOUBLE_ERROR& d,const T& t){return DOUBLE_ERROR("SUB("+d.val+","+std::to_string(t)+")",d.dbl-t,d.qd-t);}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
DOUBLE_ERROR operator-(const T& t,const DOUBLE_ERROR& d){return DOUBLE_ERROR("SUB("+std::to_string(t)+","+d.val+")",t-d.dbl,t-d.qd);}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
DOUBLE_ERROR operator*(const DOUBLE_ERROR& d,const T& t){return DOUBLE_ERROR("MUL("+d.val+","+std::to_string(t)+")",d.dbl*t,d.qd*t);}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
DOUBLE_ERROR operator*(const T& t,const DOUBLE_ERROR& d){return d*t;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
DOUBLE_ERROR operator/(const DOUBLE_ERROR& d,const T& t){return DOUBLE_ERROR("DIV("+d.val+","+std::to_string(t)+")",d.dbl/t,d.qd/t);}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
DOUBLE_ERROR operator/(const T& t,const DOUBLE_ERROR& d){return DOUBLE_ERROR("DIV("+std::to_string(t)+","+d.val+")",t/d.dbl,t/d.qd);}

template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator>(const T& t,const DOUBLE_ERROR& d){if((t>d.dbl)!=(t>d.qd)) LOG::printf("Conditional Disagreement:>\n\t%P\n\t%P\n",t,d.val);return t>d.qd;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator>(const DOUBLE_ERROR& d,const T& t){if((d.dbl>t)!=(d.qd>t)) LOG::printf("Conditional Disagreement:>\n\t%P\n\t%P\n",d.val,t);return d.qd>t;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator<(const T& t,const DOUBLE_ERROR& d){if((t<d.dbl)!=(t<d.qd)) LOG::printf("Conditional Disagreement:>\n\t%P\n\t%P\n",t,d.val);return t<d.qd;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator<(const DOUBLE_ERROR& d,const T& t){if((d.dbl<t)!=(d.qd<t)) LOG::printf("Conditional Disagreement:>\n\t%P\n\t%P\n",d.val,t);return d.qd<t;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator<=(const T& t,const DOUBLE_ERROR& d){if((t<= d.dbl)!=(t<=d.qd)) LOG::printf("Conditional Disagreement:<=\n\t%P\n\t%P\n",t,d.val);return t<=d.qd;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator<=(const DOUBLE_ERROR& d,const T& t){if((d.dbl<= t)!=(d.qd<=t)) LOG::printf("Conditional Disagreement:<=\n\t%P\n\t%P\n",d.val,t);return d.qd<=t;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator>=(const T& t,const DOUBLE_ERROR& d){if((t >= d.dbl)!=(t>=d.qd)) LOG::printf("Conditional Disagreement:>=\n\t%P\n\t%P\n",t,d.val);return t>=d.qd;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator>=(const DOUBLE_ERROR& d,const T& t){if((d.dbl >= t)!=(d.qd>=t)) LOG::printf("Conditional Disagreement:>=\n\t%P\n\t%P\n",d.val,t);return d.qd>=t;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator==(const T& t,const DOUBLE_ERROR& d){if((t == d.dbl)!=(t==d.qd)) LOG::printf("Conditional Disagreement:==\n\t%P\n\t%P\n",t,d.val);return t==d.qd;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator==(const DOUBLE_ERROR& d,const T& t){return t == d;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator!=(const T& t,const DOUBLE_ERROR& d){if((t != d.dbl)!=(t!=d.qd)) LOG::printf("Conditional Disagreement:!=\n\t%P\n\t%P\n",t,d.val);return t!=d.qd;}
template<class T,class = typename std::enable_if<std::is_arithmetic<T>::value>::type>
bool operator!=(const DOUBLE_ERROR& d,const T& t){return t != d;}
 
DOUBLE_ERROR pow(const DOUBLE_ERROR& d1,const DOUBLE_ERROR& d2){return DOUBLE_ERROR("POW("+d1.val+","+d2.val + ")",std::pow(d1.dbl,d2.dbl),powq(d1.qd,d2.qd));}
DOUBLE_ERROR abs(const DOUBLE_ERROR& d){ return DOUBLE_ERROR("std::abs("+d.val+")",std::abs(d.dbl),fabsq(d.qd)); } 
#define UNARY_MATH(X) \
DOUBLE_ERROR X(const DOUBLE_ERROR& d){ return DOUBLE_ERROR("std::X("+d.val+")",std::X(d.dbl),X##q(d.qd)); } 
UNARY_MATH(fabs);
UNARY_MATH(sqrt);
UNARY_MATH(log);
UNARY_MATH(exp);
UNARY_MATH(sin);
UNARY_MATH(cos);
UNARY_MATH(tan);
UNARY_MATH(asin);
UNARY_MATH(acos);
UNARY_MATH(atan);
}
#endif
