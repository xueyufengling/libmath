#ifndef _MATH_OP
#define _MATH_OP

#include <math.h>
#include <quadmath.h>

/**
 * 算子及基本函数族
 */
namespace math
{
/**
 * 平方根
 */
template<typename _T>
inline _T sqrt(_T value)
{
	return (_T)std::sqrt(value);
}

template<>
inline __float128 sqrt(__float128 value)
{
	return sqrtq(value);
}

/**
 * 绝对值
 */
template<typename _T>
inline _T abs(_T value)
{
	return (_T)std::abs(value);
}

template<>
inline __float128 abs(__float128 value)
{
	return fabsq(value);
}

/**
 * 最大、最小值
 */
template<typename _T>
inline _T max(_T value1, _T value2)
{
	return (_T)std::max(value1, value2);
}

template<>
inline __float128 max(__float128 value1, __float128 value2)
{
	return fmaxq(value1, value2);
}

template<typename _T>
inline _T min(_T value1, _T value2)
{
	return (_T)std::min(value1, value2);
}

template<>
inline __float128 min(__float128 value1, __float128 value2)
{
	return fminq(value1, value2);
}

/**
 * 三角函数族
 */
//正弦
template<typename _T>
inline _T sin(_T value)
{
	return (_T)std::sin(value);
}

template<>
inline __float128 sin(__float128 value)
{
	return sinq(value);
}

template<typename _T>
inline _T arcsin(_T value)
{
	return (_T)std::asin(value);
}

template<>
inline __float128 arcsin(__float128 value)
{
	return asinq(value);
}

//余弦
template<typename _T>
inline _T cos(_T value)
{
	return (_T)std::cos(value);
}

template<>
inline __float128 cos(__float128 value)
{
	return cosq(value);
}

template<typename _T>
inline _T arccos(_T value)
{
	return (_T)std::acos(value);
}

template<>
inline __float128 arccos(__float128 value)
{
	return acosq(value);
}

//正切
template<typename _T>
inline _T tan(_T value)
{
	return (_T)std::tan(value);
}

template<>
inline __float128 tan(__float128 value)
{
	return tanq(value);
}

template<typename _T>
inline _T arctan(_T value)
{
	return (_T)std::atan(value);
}

template<>
inline __float128 arctan(__float128 value)
{
	return atanq(value);
}

//余切
template<typename _T>
inline _T cot(_T value)
{
	return _T(1) / math::tan(value);
}

template<typename _T>
inline _T arccot(_T value)
{
	return (_T)std::atan2(value);
}

template<>
inline __float128 arccot(__float128 value)
{
	return atan2q(1, value);
}

//正割
template<typename _T>
inline _T sec(_T value)
{
	return _T(1) / math::cos(value);
}

template<typename _T>
inline _T arcsec(_T value)
{
	return math::arccos(_T(1) / value);
}

//余割
template<typename _T>
inline _T csc(_T value)
{
	return _T(1) / math::sin(value);
}

template<typename _T>
inline _T arccsc(_T value)
{
	return math::arcsin(_T(1) / value);
}

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

#endif//_MATH_OP
