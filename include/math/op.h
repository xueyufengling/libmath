#ifndef _MATH_OP
#define _MATH_OP

#include <math.h>
#include <quadmath.h>

/**
 * 算子及函数族
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
	static _T value(size_t i);
};

template<typename _T>
inline _T basis(size_t i)
{
	return __basis_impl<_T>::value(i);
}

}

#endif//_MATH_OP
