#include "EXECUTE_HELPER.h"
#include "CACHED_ELIMINATION_MATRIX.h"
//#####################################################################
// Function Execute_Helper
//#####################################################################
namespace PhysBAM
{
void Execute_Helper(CACHED_ELIMINATION_MATRIX<double>* cem,int num_threads)
{
    cem->Execute_Jobs(num_threads);
}
}