#ifndef _MATH_FUNCTION
#define _MATH_FUNCTION

#include <tplmp/base.h>

namespace math
{
/**
 * 函数
 */
template<typename _RetType, typename ..._ArgTypes>
using function = _RetType(*)(_ArgTypes ...);

}

#endif//_MATH_FUNCTION
