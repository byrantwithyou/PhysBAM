#ifndef __EXECUTE_HELPER_H__
#define __EXECUTE_HELPER_H__
namespace PhysBAM
{
template<class T> class CACHED_ELIMINATION_MATRIX;
//void Execute_Helper(CACHED_ELIMINATION_MATRIX<float>* cem,int num_threads);
void Execute_Helper(CACHED_ELIMINATION_MATRIX<double>* cem,int num_threads);
}
#endif