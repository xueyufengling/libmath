#ifndef _MATH_COMPLEX
#define _MATH_COMPLEX

#include <tplmp/base.h>

#include <math.h>
#include <quadmath.h>
#include <sstream>

#include <math/constant.h>

/**
 * 复数
 */
namespace math
{

template<typename _T>
class complex
{
public:
	/**
	 * @brief 实数部分
	 */
	_T re;

	/**
	 * @brief 虚数部分
	 */
	_T im;

	template<typename _T2>
	inline operator complex<_T2>()
	{
		return
		{	(_T2)re, (_T2)im};
	}

	inline complex<_T>& operator=(const complex<_T>& other) = default;

	template<typename _T2>
	inline complex<_T>& operator=(const complex<_T2>& other)
	{
		re = (_T)other.re;
		im = (_T)other.im;
		return *this;
	}

	template<typename _T2>
	inline complex<_T>& operator=(_T2 r)
	{
		re = (_T)r;
		im = _T(0);
		return *this;
	}

	template<typename _T2>
	inline complex<_T>& operator+=(const complex<_T2>& other)
	{
		re += other.re;
		im += other.im;
		return *this;
	}

	template<typename _T2>
	inline complex<_T>& operator+=(_T2 r)
	{
		re += r;
		return *this;
	}

	inline complex<_T>& operator-=(const complex& other)
	{
		re -= other.re;
		im -= other.im;
		return *this;
	}

	template<typename _T2>
	inline complex<_T>& operator-=(_T2 r)
	{
		re -= r;
		return *this;
	}

	template<typename _T2>
	inline complex<_T>& operator*=(const complex<_T2>& other)
	{
		re = re * other.re - im * other.im;
		im = re * other.im + im * other.re;
		return *this;
	}

	template<typename _T2>
	inline complex<_T>& operator*=(_T2 r)
	{
		re *= r;
		im *= r;
		return *this;
	}

	template<typename _T2>
	inline complex<_T>& operator/=(const complex<_T2>& other)
	{
		auto sqr = other.re * other.re + other.im * other.im;
		re = (_T)((re * other.re + im * other.im) / sqr);
		im = (_T)((im * other.re - re * other.im) / sqr);
		return *this;
	}

	inline complex<_T>& operator/=(_T r)
	{
		re /= r;
		im /= r;
		return *this;
	}

	inline complex<_T> operator+() const
	{
		return *this;
	}

	inline complex<_T> operator-() const
	{
		return
		{	-re, -im};
	}

	template<typename _T2>
	inline bool operator==(const complex<_T2>& other) const
	{
		return re == other.re && im == other.im;
	}

	template<typename _T2>
	inline bool operator!=(const complex<_T2>& other) const
	{
		return re != other.re || im != other.im;
	}

	/**
	 * @brief 共轭
	 */
	inline complex<_T> conj() const
	{
		return
		{	re, -im};
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << re;
		if(im > 0)
			oss << '+';
		if(im != 0)
			oss << im << 'i';
		return oss.str();
	}
};

template<typename _T>
class complex<complex<_T> >
{
	static_assert(false, "invalid complex type");
};

template<typename _T>
inline std::ostream& operator<<(std::ostream& os, const complex<_T>& z)
{
	os << ((std::string)z);
	return os;
}

//实数 op 复数
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() + tplmp::decl<_T2>::val())>
inline complex<_Result> operator+(_T1 lhs, const complex<_T2>& rhs)
{
	return
	{	lhs + rhs.re, rhs.im};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() - tplmp::decl<_T2>::val())>
inline complex<_Result> operator-(_T1 lhs, const complex<_T2>& rhs)
{
	return
	{	lhs - rhs.re, -rhs.im};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline complex<_Result> operator*(_T1 lhs, const complex<_T2>& rhs)
{
	return
	{	lhs * rhs.re, lhs * rhs.im};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline complex<_Result> operator/(_T1 lhs, const complex<_T2>& rhs)
{
	auto sqr = rhs.re * rhs.re + rhs.im * rhs.im;
	return
	{
		_Result(lhs * rhs.re / sqr),
		_Result(-lhs * rhs.im / sqr)
	};
}

//复数 op 实数
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() + tplmp::decl<_T2>::val())>
inline complex<_Result> operator+(const complex<_T1>& lhs, _T2 rhs)
{
	return
	{	lhs.re + rhs, lhs.im};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() - tplmp::decl<_T2>::val())>
inline complex<_Result> operator-(const complex<_T1>& lhs, _T2 rhs)
{
	return
	{	lhs.re - rhs, lhs.im};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline complex<_Result> operator*(const complex<_T1>& lhs, _T2 rhs)
{
	return
	{	lhs.re * rhs, lhs.im * rhs};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline complex<_Result> operator/(const complex<_T1>& lhs, _T2 rhs)
{
	return
	{	lhs.re / rhs, lhs.im / rhs};
}

//复数 op 复数
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() + tplmp::decl<_T2>::val())>
inline complex<_Result> operator+(const complex<_T1>& lhs, const complex<_T2>& rhs)
{
	return
	{	lhs.re + rhs.re, lhs.im + rhs.im};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() - tplmp::decl<_T2>::val())>
inline complex<_Result> operator-(const complex<_T1>& lhs, const complex<_T2>& rhs)
{
	return
	{	lhs.re - rhs.re, lhs.im - rhs.im};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline complex<_Result> operator*(const complex<_T1>& lhs, const complex<_T2>& rhs)
{
	return
	{
		lhs.re * rhs.re - lhs.im * rhs.im,
		lhs.re * rhs.im + lhs.im * rhs.re
	};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline complex<_Result> operator/(const complex<_T1>& lhs, const complex<_T2>& rhs)
{
	auto sqr = rhs.re * rhs.re + rhs.im * rhs.im;
	return
	{
		_Result( (lhs.re * rhs.re + lhs.im * rhs.im) / sqr),
		_Result( (lhs.im * rhs.re - lhs.re * rhs.im) / sqr)
	};
}

//值域(-π, π]，从x轴正向逆时针为正方向
template<typename _T>
inline _T arctan(_T y, _T x)
{
	return (_T)std::atan2(y, x);
}

template<>
inline __float128 arctan(__float128 y, __float128 x)
{
	return atan2q(y, x);
}

template<typename _T>
struct __complex_op_impl
{
	typedef _T type;
	typedef complex<_T> complex_type;

	// 实数
	inline static constexpr _T Re(_T r)
	{
		return r;
	}

	inline static constexpr _T Im(_T r)
	{
		return 0;
	}

	inline static constexpr _T arg(_T r)
	{
		return r < 0 ? constant<_T>::pi : 0;
	}

	/**
	 * @brief 纯虚数
	 */
	inline static constexpr complex<_T> iz(_T r)
	{
		return
		{	0, r};
	}
};

template<typename _T>
struct __complex_op_impl<complex<_T> >
{
	typedef complex<_T> type;
	typedef complex<_T> complex_type;

	// 复数特化
	inline static constexpr _T Re(const complex<_T>& z)
	{
		return z.re;
	}

	inline static constexpr _T Im(const complex<_T>& z)
	{
		return z.im;
	}

	inline static constexpr _T arg(const complex<_T>& z)
	{
		return arctan(z.im, z.re);
	}

	inline static constexpr complex<_T> iz(const complex<_T>& z)
	{
		return
		{	-z.im, z.re};
	}
};

/**
 * @brief 虚数单位字面量。
 * 		  字面量重载不能是模板
 */
#define __def_imaginary_unit__(type)\
inline constexpr math::complex<type> operator "" _i(type im)\
{\
	return {0, im};\
}

__def_imaginary_unit__(unsigned long long)
__def_imaginary_unit__(long double)

#undef __def_imaginary_unit__

/**
 * @brief 对复数取实部
 */
template<typename _T>
inline constexpr typename __complex_op_impl<_T>::type Re(const _T& z)
{
	return __complex_op_impl<_T>::Re(z);
}

/**
 * @brief 对复数取虚部
 */
template<typename _T>
inline constexpr typename __complex_op_impl<_T>::type Im(const _T& z)
{
	return __complex_op_impl<_T>::Im(z);
}

/**
 * @brief 复数的幅角，值域(-π, π]
 */
template<typename _T>
inline constexpr typename __complex_op_impl<_T>::type arg(const _T& z)
{
	return __complex_op_impl<_T>::arg(z);
}

/**
 * @brief 高效计算i*z的结果
 */
template<typename _T>
inline constexpr
typename __complex_op_impl<_T>::complex_type iz(const _T& z)
{
	return __complex_op_impl<_T>::iz(z);
}

}

#endif//_MATH_COMPLEX
