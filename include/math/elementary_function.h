#ifndef _MATH_ELEMENTARYFUNCTION
#define _MATH_ELEMENTARYFUNCTION

#include <math.h>
#include <quadmath.h>

#include <math/complex.h>

/**
 * 实数域及复数域的初等函数族及常用函数
 */

namespace math
{
//实数域

/**
 *平方根
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
 *绝对值
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
 *自然指数
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
 *自然对数
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
 *常用对数
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
 *任意底对数
 */
template<typename _T>
inline _T log(_T a, _T x)
{
	return (_T)(math::ln(x) / math::ln(a));
}

/**
 *幂函数
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
 *最大、最小值
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
 *三角函数族
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
 *双曲三角函数族
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

//复数域

/**
 *平方根
 *sqrt(z) = sqrt(|z|)*e^(i arg(z)/2)
 *		   = sqrt(|z|)*[cos(arg(z)/2)+i sin(arg(z)/2)]
 */
template<typename _T>
inline complex<_T> sqrt(const complex<_T>& z)
{
	_T sqrt_r = sqrt(sqrt(z.re * z.re + z.im * z.im));
	_T arg_2 = arg(z) / 2;
	return
	{	sqrt_r*cos(arg_2), sqrt_r*sin(arg_2)};
}

/**
 *模
 *|z| = sqrt(Re(z)²+Im(z)²)
 */
template<typename _T>
inline _T abs(const complex<_T>& z)
{
	return sqrt(z.re * z.re + z.im * z.im);
}

/**
 *自然指数
 *e^z = e^(Re(z))*e^(i Im(z))
 *	 	 e^(Re(z))*(cos(Im(z))+i sin(Im(z)))
 */
template<typename _T>
inline complex<_T> exp(const complex<_T>& z)
{
	return
	{	exp(z.re)*cos(z.im), exp(z.re)*sin(z.im)};
}

/**
 *自然对数主值
 *ln(z) = ln|z|+i arg(z)
 */
template<typename _T>
inline complex<_T> ln(const complex<_T>& z)
{
	return
	{	ln(abs(z)), arg(z)};
}

/**
 *常用对数主值
 *lg(z) = ln(z)/ln(10)
 */
template<typename _T>
inline complex<_T> lg(const complex<_T>& z)
{
	return ln(z) / constant<_T>::ln10;
}

/**
 *任意底对数主值
 *log_a(z) = ln(z)/ln(a)
 */
template<typename _T>
inline complex<_T> log(const complex<_T>& a, const complex<_T>& z)
{
	return ln(z) / ln(a);
}

template<typename _T>
inline complex<_T> log(const complex<_T>& a, _T x)
{
	return ln(x) / ln(a);
}

template<typename _T>
inline complex<_T> log(_T a, const complex<_T>& z)
{
	return ln(z) / ln(a);
}

/**
 *幂函数
 *z^a = e^(a*ln(z))
 */
template<typename _T>
inline complex<_T> pow(const complex<_T>& z, const complex<_T>& a)
{
	return exp(a * ln(z));
}

template<typename _T>
inline complex<_T> pow(const complex<_T>& z, _T n)
{
	return exp(n * ln(z));
}

template<typename _T>
inline complex<_T> pow(_T a, const complex<_T>& x)
{
	return exp(x * ln(a));
}

/**
 *三角函数族
 */

/**
 *正弦
 *sin(z) = (e^(iz)-e^(-iz))/(2i)
 */
template<typename _T>
inline complex<_T> sin(const complex<_T>& z)
{
	complex<_T> iz = math::iz(z);
	return (exp(iz) - exp(-iz)) / (2_i);
}

/**
 *反正弦
 *arcsin(z) = -i ln(iz+sqrt(1-z²))
 */
template<typename _T>
inline complex<_T> arcsin(const complex<_T>& z)
{
	return -math::iz(ln(math::iz(z) + sqrt(1 - z * z)));
}

/**
 *余弦
 *cos(z) = (e^(iz)+e^(-iz))/2
 */
template<typename _T>
inline complex<_T> cos(const complex<_T>& z)
{
	complex<_T> iz = math::iz(z);
	return (exp(iz) + exp(-iz)) / _T(2);
}

/**
 *反余弦
 *arccos(z) = -i ln(z+i sqrt(1-z²))
 */
template<typename _T>
inline complex<_T> arccos(const complex<_T>& z)
{
	return -math::iz(ln(z + math::iz(sqrt(1 - z * z))));
}

/**
 *正切
 *tan(z) = sin(z)/cos(z)
 */
template<typename _T>
inline complex<_T> tan(const complex<_T>& z)
{
	return sin(z) / cos(z);
}

/**
 *反正切
 *arctan(z) = (i/2)*(ln(1-i*z)-ln(1+i*z))
 */
template<typename _T>
inline complex<_T> arctan(const complex<_T>& z)
{
	complex<_T> iz = math::iz(z);
	return 0.5_i * (ln(1 - iz) - ln(1 + iz));
}

/**
 *余切
 *cot(z) = cos(z)/sin(z)
 */
template<typename _T>
inline complex<_T> cot(const complex<_T>& z)
{
	return cos(z) / sin(z);
}

/**
 *反余切
 *arccot(z) = arctan(1/z)
 */
template<typename _T>
inline complex<_T> arccot(const complex<_T>& z)
{
	return arctan(1 / z);
}

/**
 *正割
 *sec(z) = 1/cos(z)
 */
template<typename _T>
inline complex<_T> sec(const complex<_T>& z)
{
	return 1 / cos(z);
}

/**
 *反正割
 *arcsec(z) = arccos(1/z)
 */
template<typename _T>
inline complex<_T> arcsec(const complex<_T>& z)
{
	return arccos(1 / z);
}

/**
 *余割
 *csc(z) = 1/sin(z)
 */
template<typename _T>
inline complex<_T> csc(const complex<_T>& z)
{
	return 1 / sin(z);
}

/**
 *反余割
 *arccsc(z) = arcsin(1/z)
 */
template<typename _T>
inline complex<_T> arccsc(const complex<_T>& z)
{
	return arcsin(1 / z);
}

/**
 *双曲三角函数族
 */

/**
 *双曲正弦
 *sinh(z) = (e^z-e^(-z))/2
 */
template<typename _T>
inline complex<_T> sinh(const complex<_T>& z)
{
	return (exp(z) - exp(-z)) / _T(2);
}

/**
 *反双曲正弦
 *arcsinh(z) = ln(z+sqrt(z²+1))
 */
template<typename _T>
inline complex<_T> arcsinh(const complex<_T>& z)
{
	return ln(z + sqrt(z * z + 1));
}

/**
 *双曲余弦
 *cosh(z) = (e^z+e^(-z))/2
 */
template<typename _T>
inline complex<_T> cosh(const complex<_T>& z)
{
	return (exp(z) + exp(-z)) / _T(2);
}

/**
 *反双曲余弦
 *arccosh(z) = ln(z+sqrt(z²-1))
 */
template<typename _T>
inline complex<_T> arccosh(const complex<_T>& z)
{
	return ln(z + sqrt(z * z - 1));
}

/**
 *双曲正切
 *tanh(z) = sinh(z)/cosh(z)
 */
template<typename _T>
inline complex<_T> tanh(const complex<_T>& z)
{
	return sinh(z) / cosh(z);
}

/**
 *反双曲正切
 *arctanh(z) = (1/2)*(ln(1+z)-ln(1-z))
 */
template<typename _T>
inline complex<_T> arctanh(const complex<_T>& z)
{
	return (ln(1 + z) - ln(1 - z)) / _T(2);
}

/**
 *双曲余切
 *coth(z) = cosh(z)/sinh(z)
 */
template<typename _T>
inline complex<_T> coth(const complex<_T>& z)
{
	return cosh(z) / sinh(z);
}

/**
 *反双曲余切
 *arccoth(z) = (ln(1+1/z)-ln(1-1/z))/2
 */
template<typename _T>
inline complex<_T> arccoth(const complex<_T>& z)
{
	return (ln(1 + 1 / z) - ln(1 - 1 / z)) / _T(2);
}

/**
 *双曲正割
 *sech(z) = 1/cosh(z)
 */
template<typename _T>
inline complex<_T> sech(const complex<_T>& z)
{
	return 1 / cosh(z);
}

/**
 *反双曲正割
 *arcsech(z) = ln((1+sqrt(1-z²))/z)
 */
template<typename _T>
inline complex<_T> arcsech(const complex<_T>& z)
{
	return ln((1 + sqrt(1 - z * z)) / z);
}

/**
 *双曲余割
 *csch(z) = 1/sinh(z)
 */
template<typename _T>
inline complex<_T> csch(const complex<_T>& z)
{
	return 1 / sinh(z);
}

/**
 *反双曲余割
 *arccsch(z) = ln((1+sqrt(1+z²))/z)
 */
template<typename _T>
inline complex<_T> arccsch(const complex<_T>& z)
{
	return ln((1 + sqrt(1 + z * z)) / z);
}

/**
 * 其他常用函数
 */

/**
 * 阶乘，仅自然数有定义
 */
template<typename _T>
inline constexpr _T factorial(_T k)
{
	return k == 0 ? 1 : factorial<_T>(k - 1) * k;
}

template<typename _T, size_t _K>
inline constexpr _T factorial(void)
{
	return _K == 0 ? 1 : factorial<_T, _K - 1> * _K;
}

}

#endif//_MATH_ELEMENTARYFUNCTION
