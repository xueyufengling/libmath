#ifndef _MATH_ALGEBRA
#define _MATH_ALGEBRA

#include <math/elementary_function.h>

/**
 * 算子及基本函数族
 */
namespace math
{
struct add;
struct mul;

/**
 * 恒等元
 */
template<typename _T, typename _Op>
struct __identity_impl
{
	inline static _T value()
	{
		return _T(1);
	}
};

template<typename _T, typename _Op>
inline _T identity()
{
	return __identity_impl<_T, _Op>::value();
}

/**
 * 零元
 */
template<typename _T, typename _Op>
struct __zero_impl
{
	inline static _T value()
	{
		return _T(0);
	}
};

template<typename _T, typename _Op>
inline _T zero()
{
	return __zero_impl<_T, _Op>::value();
}

template<typename T>
struct precision
{
	inline static T& default_val()
	{
		static T value = T(0);
		return value;
	}
};

template<>
struct precision<float>
{
	inline static float& default_val()
	{
		static float value = 1e-6f;
		return value;
	}
};

template<>
struct precision<double>
{
	inline static double& default_val()
	{
		static double value = 1e-14;
		return value;
	}
};

template<>
struct precision<__float128>
{
	inline static __float128& default_val()
	{
		static __float128 value = 1e-28;
		return value;
	}
};

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

}

#endif//_MATH_ALGEBRA
