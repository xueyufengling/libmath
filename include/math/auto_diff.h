#ifndef _MATH_AUTODIFF
#define _MATH_AUTODIFF

#include <tplmp/base.h>
#include <math/elementary_function.h>
#include <math/vector.h>

/**
 * @brief 自动微分定义及初等函数族
 */

namespace math
{
namespace auto_diff
{
/**
 * @brief 计算中间值，包括函数值和偏导数值
 * _Dim为线性无关的基底数量，本类缓存由基函数支成的点的值
 */
template<size_t _Dim, typename _T>
struct calc_val
{
	_T value; //函数值
	math::vector<_Dim, _T> grad; //梯度，即偏导数向量

	inline calc_val(_T v, math::vector<_Dim, _T> g) :
			value(v), grad(g)
	{
	}

	inline operator _T()
	{
		return value;
	}

	/**
	 * @brief 提取偏导数
	 */
	inline calc_val<1, _T> partial(size_t i) const
	{
		return calc_val<1, _T>(value, {grad[i]});
	}

	template<typename _T2>
	inline operator calc_val<_Dim, _T2>() const
	{
		return calc_val<_Dim, _T2>((_T2)value, (math::vector<_Dim, _T2>)grad);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline calc_val<_Dim, _Result> operator+(const calc_val<_Dim, _T2>& other) const
	{
		//d(x+y)=dx+dy
		return calc_val(value + other.value, grad + other.grad);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline calc_val<_Dim, _Result> operator-(const calc_val<_Dim, _T2>& other) const
	{
		//d(x-y)=dx-dy
		return calc_val(value - other.value, grad - other.grad);
	}

	inline calc_val<_Dim, _T> operator-() const
	{
		return calc_val<_Dim, _T>(-value, -grad);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline calc_val<_Dim, _Result> operator*(const calc_val<_Dim, _T2>& other) const
	{
		//d(xy)=xdy+ydx
		return calc_val(value * other.value, grad * other.value + value * other.grad);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline calc_val<_Dim, _Result> operator/(const calc_val<_Dim, _T2>& other) const
	{
		//d(x/y)=(ydx-xdy)/y²
		return calc_val(value / other.value, (grad * other.value - value * other.grad) / (other.value * other.value));
	}
};

//自变量函数
template<size_t _Dim, typename _T>
inline calc_val<_Dim, _T> x(_T val)
{
	return calc_val<_Dim, _T>(val, math::vector<_Dim, _T>::one());
}

template<typename _T>
inline calc_val<1, _T> exp(const calc_val<1, _T>& x)
{
	_T exp_val = math::exp(x.value);
	return calc_val<1, _T>(exp_val, exp_val * x.grad);
}

template<typename _T>
inline calc_val<1, _T> ln(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::ln(x.value), x.grad / x.value);
}

template<typename _T>
inline calc_val<1, _T> lg(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::lg(x.value), x.grad / (x.value * math::ln(_T(10))));
}

template<typename _T>
inline calc_val<1, _T> sqrt(const calc_val<1, _T>& x)
{
	double sqrt_val = math::sqrt(x.value);
	return calc_val<1, _T>(sqrt_val, x.grad / (_T(2) * sqrt_val));
}

template<typename _T>
inline calc_val<2, _T> pow(const calc_val<1, _T>& a, const calc_val<1, _T>& x)
{
	_T pow_val = math::pow(a.value, x.value);
	return calc_val<2, _T>(pow_val, {
			x.value * math::pow(a.value, x.value - _T(1)), // ∂f/∂a
			math::ln(a.value) * pow_val   // ∂f/∂x
	});
}

template<typename _T>
inline calc_val<1, _T> abs(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::abs(x.value), (x.value >= _T(0) ? _T(1) : _T(-1)) * x.grad);
}

template<typename _T>
inline calc_val<1, _T> sin(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::sin(x.value), math::cos(x.value) * x.grad);
}

template<typename _T>
inline calc_val<1, _T> cos(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::cos(x.value), -math::sin(x.value) * x.grad);
}

template<typename _T>
inline calc_val<1, _T> tan(const calc_val<1, _T>& x)
{
	//d(tanx)/dx = sec²x
	_T cos_val = math::cos(x.value);
	return calc_val<1, _T>(math::tan(x.value), x.grad / (cos_val * cos_val));
}

template<typename _T>
inline calc_val<1, _T> arcsin(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::arcsin(x.value), x.grad / math::sqrt(_T(1) - x.value * x.value));
}

template<typename _T>
inline calc_val<1, _T> arccos(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::arccos(x.value), -x.grad / math::sqrt(_T(1) - x.value * x.value));
}

template<typename _T>
inline calc_val<1, _T> arctan(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::arctan(x.value), x.grad / (_T(1) + x.value * x.value));
}

template<typename _T>
inline calc_val<1, _T> sinh(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::sinh(x.value), math::cosh(x.value) * x.grad);
}

template<typename _T>
inline calc_val<1, _T> cosh(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::cosh(x.value), math::sinh(x.value) * x.grad);
}

template<typename _T>
inline calc_val<1, _T> tanh(const calc_val<1, _T>& x)
{
	double tanh_val = math::tanh(x.value);
	return calc_val<1, _T>(tanh_val, (_T(1) - tanh_val * tanh_val) * x.grad);
}

template<typename _T>
inline calc_val<1, _T> arcsinh(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::arcsinh(x.value), x.grad / math::sqrt(x.value * x.value + _T(1)));
}

template<typename _T>
inline calc_val<1, _T> arccosh(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::arccosh(x.value), x.grad / math::sqrt(x.value * x.value - _T(1)));
}

template<typename _T>
inline calc_val<1, _T> arctanh(const calc_val<1, _T>& x)
{
	return calc_val<1, _T>(math::arctanh(x.value), x.grad / (_T(1) - x.value * x.value));
}

}
}

#endif//_MATH_AUTODIFF
