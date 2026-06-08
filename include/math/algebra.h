#ifndef _MATH_ALGEBRA
#define _MATH_ALGEBRA

#include <math/elementary_function.h>

/**
 * 算子及基本函数族
 */
namespace math
{
/**
 * 恒等元
 */
template<typename _T>
struct __identity_impl
{
	static _T value()
	{
		return _T(1);
	}
};

template<typename _T>
inline _T identity()
{
	return __identity_impl<_T>::value();
}

/**
 * 基
 */
template<typename _T>
struct __basis_impl
{
	template<typename ..._Indices>
	static _T value(_Indices ... indices)
	{
		static_assert(false, "basis not defined");
	}
};

template<typename _T, typename ..._Indices>
inline _T basis(_Indices ... indices)
{
	return __basis_impl<_T>::value(indices...);
}

/**
 * 可带容差地判断是否是零元
 */
template<typename _T, typename ..._Args>
inline bool is_zero(_T value, _Args ... args)
{
	return value == 0;
}

template<typename _T>
inline bool is_zero(_T value, _T eps)
{
	return math::abs(value) < eps;
}

}

#endif//_MATH_ALGEBRA
