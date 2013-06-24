#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <cmath>
#include <iomanip>
#include "ACCURACY_INFO.h"
#include "BOUNDARY_CONDITIONS.h"
#include "OBJECTS_COMMON.h"

template<class TV> OBJECTS_COMMON<TV>::
OBJECTS_COMMON()
    :bc(0)
{
}

template<class TV> OBJECTS_COMMON<TV>::
~OBJECTS_COMMON()
{
    delete bc;
}

template struct OBJECTS_COMMON<VECTOR<double,1> >;
template struct OBJECTS_COMMON<VECTOR<double,2> >;
