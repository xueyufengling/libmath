#ifndef _MATH_ELEMENTARYFUNCTION
#define _MATH_ELEMENTARYFUNCTION

#include <math.h>
#include <quadmath.h>

/**
 * 初等函数族及常用函数
 */
namespace math
{
/**
 * 平方根
 */
template<typename _T>
inline _T sqrt(_T x)
{
	return (_T)std::sqrt(x);
}

template<>
inline __float128 sqrt(__float128 x)
{
	return sqrtq(x);
}

/**
 * 绝对值
 */
template<typename _T>
inline _T abs(_T x)
{
	return (_T)std::abs(x);
}

template<>
inline __float128 abs(__float128 x)
{
	return fabsq(x);
}

/**
 * 自然指数
 */
template<typename _T>
inline _T exp(_T x)
{
	return (_T)std::exp(x);
}

template<>
inline __float128 exp(__float128 x)
{
	return expq(x);
}

/**
 * 自然对数
 */
template<typename _T>
inline _T ln(_T x)
{
	return (_T)std::log(x);
}

template<>
inline __float128 ln(__float128 x)
{
	return logq(x);
}

/**
 * 常用对数
 */
template<typename _T>
inline _T lg(_T x)
{
	return (_T)std::log10(x);
}

template<>
inline __float128 lg(__float128 x)
{
	return log10q(x);
}

/**
 * 任意底对数
 */
template<typename _T>
inline _T log(_T a, _T x)
{
	return (_T)(math::ln(x) / math::ln(a));
}

/**
 * 幂函数
 */
template<typename _T>
inline _T pow(_T a, _T x)
{
	return (_T)std::pow(a, x);
}

template<>
inline __float128 pow(__float128 a, __float128 x)
{
	return powq(a, x);
}

/**
 * 最大、最小值
 */
template<typename _T>
inline _T max(_T x1, _T x2)
{
	return (_T)std::max(x1, x2);
}

template<>
inline __float128 max(__float128 x1, __float128 x2)
{
	return fmaxq(x1, x2);
}

template<typename _T>
inline _T min(_T x1, _T x2)
{
	return (_T)std::min(x1, x2);
}

template<>
inline __float128 min(__float128 x1, __float128 x2)
{
	return fminq(x1, x2);
}

/**
 * 三角函数族
 */

//正弦
template<typename _T>
inline _T sin(_T x)
{
	return (_T)std::sin(x);
}

template<>
inline __float128 sin(__float128 x)
{
	return sinq(x);
}

template<typename _T>
inline _T arcsin(_T x)
{
	return (_T)std::asin(x);
}

template<>
inline __float128 arcsin(__float128 x)
{
	return asinq(x);
}

//余弦
template<typename _T>
inline _T cos(_T x)
{
	return (_T)std::cos(x);
}

template<>
inline __float128 cos(__float128 x)
{
	return cosq(x);
}

template<typename _T>
inline _T arccos(_T x)
{
	return (_T)std::acos(x);
}

template<>
inline __float128 arccos(__float128 x)
{
	return acosq(x);
}

//正切
template<typename _T>
inline _T tan(_T x)
{
	return (_T)std::tan(x);
}

template<>
inline __float128 tan(__float128 x)
{
	return tanq(x);
}

template<typename _T>
inline _T arctan(_T x)
{
	return (_T)std::atan(x);
}

template<>
inline __float128 arctan(__float128 x)
{
	return atanq(x);
}

//余切
template<typename _T>
inline _T cot(_T x)
{
	return _T(1) / math::tan(x);
}

template<typename _T>
inline _T arccot(_T x)
{
	return (_T)std::atan2(x);
}

template<>
inline __float128 arccot(__float128 x)
{
	return atan2q(1, x);
}

//正割
template<typename _T>
inline _T sec(_T x)
{
	return _T(1) / math::cos(x);
}

template<typename _T>
inline _T arcsec(_T x)
{
	return math::arccos(_T(1) / x);
}

//余割
template<typename _T>
inline _T csc(_T x)
{
	return _T(1) / math::sin(x);
}

template<typename _T>
inline _T arccsc(_T x)
{
	return math::arcsin(_T(1) / x);
}

/**
 * 双曲三角函数族
 */

//双曲正弦
template<typename _T>
inline _T sinh(_T x)
{
	return (_T)std::sinh(x);
}

template<>
inline __float128 sinh(__float128 x)
{
	return sinhq(x);
}

template<typename _T>
inline _T arcsinh(_T x)
{
	return (_T)std::asinh(x);
}

template<>
inline __float128 arcsinh(__float128 x)
{
	return asinhq(x);
}

//双曲余弦
template<typename _T>
inline _T cosh(_T x)
{
	return (_T)std::cosh(x);
}

template<>
inline __float128 cosh(__float128 x)
{
	return coshq(x);
}

template<typename _T>
inline _T arccosh(_T x)
{
	return (_T)std::acosh(x);
}

template<>
inline __float128 arccosh(__float128 x)
{
	return acoshq(x);
}

//双曲正切
template<typename _T>
inline _T tanh(_T x)
{
	return (_T)std::tanh(x);
}

template<>
inline __float128 tanh(__float128 x)
{
	return tanhq(x);
}

template<typename _T>
inline _T arctanh(_T x)
{
	return (_T)std::atanh(x);
}

template<>
inline __float128 arctanh(__float128 x)
{
	return atanhq(x);
}

//双曲余切
template<typename _T>
inline _T coth(_T x)
{
	return _T(1) / math::tanh(x);
}

template<typename _T>
inline _T arccoth(_T x)
{
	return math::arctanh(_T(1) / x);
}

//双曲正割
template<typename _T>
inline _T sech(_T x)
{
	return _T(1) / math::cosh(x);
}

template<typename _T>
inline _T arcsech(_T x)
{
	return math::arccosh(_T(1) / x);
}

//双曲余割
template<typename _T>
inline _T csch(_T x)
{
	return _T(1) / math::sinh(x);
}

template<typename _T>
inline _T arccsch(_T x)
{
	return math::arcsinh(_T(1) / x);
}

}

#endif//_MATH_ELEMENTARYFUNCTION
