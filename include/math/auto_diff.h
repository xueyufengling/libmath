#ifndef _MATH_AUTODIFF
#define _MATH_AUTODIFF

#include <tplmp/base.h>
#include <tplmp/array.h>

#include <math/elementary_function.h>
#include <math/vector.h>
#include <math/matrix.h>

/**
 * @brief 自动微分定义及初等函数族
 */

namespace math
{
namespace auto_diff
{
/**
 * @brief 计算中间值，包括函数值和偏导数值。
 * 表示_InDim元_OutDim维向量函数。当_OutDim==1时退化为标量函数。
 * 每行都是各自独立的，只与对应下标的分量有关。
 */
template<size_t _OutDim, size_t _InDim, typename _T>
struct calc_point
{
	math::vector<_OutDim, _T> value; //多元函数值，即坐标，列向量
	math::matrix<_OutDim, _InDim, _T> derivative; //多元向量函数导数，即Jacobi矩阵，每个行向量都是对应坐标的梯度

	inline calc_point() = default;

	inline calc_point(math::vector<_OutDim, _T> c) :
			value(c)
	{
	}

	inline calc_point(math::matrix<_OutDim, _InDim, _T> d) :
			derivative(d)
	{
	}

	inline calc_point(math::vector<_OutDim, _T> c, math::matrix<_OutDim, _InDim, _T> d) :
			value(c), derivative(d)
	{
	}

	/**
	 * @brief 隐式转换到坐标
	 */
	inline operator math::vector<_OutDim, _T>()
	{
		return value;
	}

	/**
	 * @brief 提取指定坐标的梯度
	 */
	inline math::vector<_InDim, _T> grad(size_t i) const
	{
		return derivative[i];
	}

	/**
	 * @brief 提取偏导数的值
	 * @param i 输出分量下标
	 * @param j 输入自变量下标
	 */
	inline _T partial(size_t i, size_t j) const
	{
		return derivative[i][j];
	}

	inline math::vector<_OutDim, _T> partial(size_t j) const
	{
		return derivative.column(j);
	}

	/**
	 * Jacobi矩阵为指定列赋值
	 * 即对所有坐标的某个自变量偏导数赋值
	 */
	inline void partial_assign(size_t j, const math::vector<_OutDim, _T>& partial)
	{
		for(size_t i = 0; i < _OutDim; ++i)
			derivative.elem[i][j] = partial[i];
	}

	template<typename _T2>
	inline operator calc_point<_OutDim, _InDim, _T2>() const
	{
		return calc_point<_OutDim, _InDim, _T2>((math::vector<_OutDim, _T2>)value, (math::matrix<_OutDim, _InDim, _T2>)derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline calc_point<_OutDim, _InDim, _Result> operator+(const calc_point<_OutDim, _InDim, _T2>& other) const
	{
		// ∇(u+v)=∇u+∇v
		return calc_point(value + other.value, derivative + other.derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline calc_point<_OutDim, _InDim, _Result> operator-(const calc_point<_OutDim, _InDim, _T2>& other) const
	{
		// ∇(u-v)=∇u-∇v
		return calc_point(value - other.value, derivative - other.derivative);
	}

	inline calc_point<_OutDim, _InDim, _T> operator-() const
	{
		// ∇(-u)=-∇u
		return calc_point<_OutDim, _InDim, _T>(-value, -derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline calc_point<_OutDim, _InDim, _Result> operator*(const calc_point<_OutDim, _InDim, _T2>& other) const
	{
		// ∇(uv)=u∇v+v∇u
		calc_point cp(hadamard(value, other.value));
		for(size_t i = 0; i < _OutDim; ++i)
			cp.derivative[i] = derivative[i] * other.value[i] + other.derivative[i] * value[i];
		return cp;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline calc_point<_OutDim, _InDim, _Result> operator/(const calc_point<_OutDim, _InDim, _T2>& other) const
	{
		// ∇(u/v)=(v∇u-u∇v)/v²
		calc_point cp(value / other.value);
		for(size_t i = 0; i < _OutDim; ++i)
			cp.derivative[i] = (derivative[i] * other.value[i] - other.derivative[i] * value[i]) / (other.value[i] * other.value[i]);
		return cp;
	}
};

/**
 * @brief 自变量函数。
 * @param basis 该自变量的值作为系数对应的基底，不同基底之间各自独立，故非basis所在行的value和derivative均为0.
 * @param subscript 该自变量的下标
 */
template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> x(size_t basis, size_t subscript, _T val)
{
	//仅对应下标的值为val，其偏导数全为1
	calc_point<_OutDim, _InDim, _T> cp(math::vector<_OutDim, _T>::zero(), math::matrix<_OutDim, _InDim, _T>::zero());
	cp.value[basis] = val;
	cp.derivative[basis][subscript] = _T(1); //仅basis坐标对应的梯度中关于subscript的偏导数为1，其余偏导数均为0
	return cp;
}

//函数都只映射指定下标的分量
template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> exp(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	_T exp_val = math::exp(x.value[subscript]);
	calc_point<_OutDim, _InDim, _T> result = x;
	x.value[subscript] = exp_val;
	x.derivative[subscript] = exp_val * x.derivative[subscript];
	return result;
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> ln(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::ln(x.value), x.derivative / x.value);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> lg(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::lg(x.value), x.derivative / (x.value * math::ln(_T(10))));
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> sqrt(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	double sqrt_val = math::sqrt(x.value);
	return calc_point<_OutDim, _InDim, _T>(sqrt_val, x.derivative / (_T(2) * sqrt_val));
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> pow(size_t subscript, const calc_point<_OutDim, _InDim, _T>& a, const calc_point<_OutDim, _InDim, _T>& x)
{
	_T pow_val = math::pow(a.value, x.value);
	return calc_point<_OutDim, _InDim, _T>(pow_val, {
			x.value * math::pow(a.value, x.value - _T(1)), // ∂f/∂a
			math::ln(a.value) * pow_val   // ∂f/∂x
	});
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> abs(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::abs(x.value), (x.value >= _T(0) ? _T(1) : _T(-1)) * x.derivative);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> sin(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::sin(x.value), math::cos(x.value) * x.derivative);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> cos(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::cos(x.value), -math::sin(x.value) * x.derivative);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> tan(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	//d(tanx)/dx = sec²x
	_T cos_val = math::cos(x.value);
	return calc_point<_OutDim, _InDim, _T>(math::tan(x.value), x.derivative / (cos_val * cos_val));
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> arcsin(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::arcsin(x.value), x.derivative / math::sqrt(_T(1) - x.value * x.value));
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> arccos(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::arccos(x.value), -x.derivative / math::sqrt(_T(1) - x.value * x.value));
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> arctan(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::arctan(x.value), x.derivative / (_T(1) + x.value * x.value));
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> sinh(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::sinh(x.value), math::cosh(x.value) * x.derivative);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> cosh(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::cosh(x.value), math::sinh(x.value) * x.derivative);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> tanh(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	double tanh_val = math::tanh(x.value);
	return calc_point<_OutDim, _InDim, _T>(tanh_val, (_T(1) - tanh_val * tanh_val) * x.derivative);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> arcsinh(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::arcsinh(x.value), x.derivative / math::sqrt(x.value * x.value + _T(1)));
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> arccosh(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::arccosh(x.value), x.derivative / math::sqrt(x.value * x.value - _T(1)));
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline calc_point<_OutDim, _InDim, _T> arctanh(size_t subscript, const calc_point<_OutDim, _InDim, _T>& x)
{
	return calc_point<_OutDim, _InDim, _T>(math::arctanh(x.value), x.derivative / (_T(1) - x.value * x.value));
}

}
}

#endif//_MATH_AUTODIFF
