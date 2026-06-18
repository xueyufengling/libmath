#ifndef _MATH_AUTODIFF
#define _MATH_AUTODIFF

#include <tplmp/base.h>
#include <tplmp/array.h>

#include <math/elementary_function.h>
#include <math/tensor.h>

/**
 * @brief 自动微分定义及初等函数族
 */

namespace math
{
/**
 * @brief 计算中间值，包括函数值和偏导数值。
 * 表示_InDim元_OutDim维向量函数。当_OutDim==1时退化为标量函数。
 * 每行都是各自独立的，只与对应下标的分量有关。
 */
template<typename _T, size_t _OutDim, size_t _InDim>
struct ad_point
{
	typedef vector<_T, _OutDim> value_type;
	typedef matrix<_T, _OutDim, _InDim> derivative_type;

	value_type value; //多元函数值，即坐标，列向量
	derivative_type derivative; //多元向量函数导数，即Jacobi矩阵，每个行向量都是对应坐标的梯度

	inline ad_point() = default;

	inline ad_point(const value_type& c) :
			value(c)
	{
	}

	inline ad_point(const derivative_type& d) :
			derivative(d)
	{
	}

	inline ad_point(const value_type& c, const derivative_type& d) :
			value(c), derivative(d)
	{
	}

	/**
	 * @brief 隐式转换到坐标
	 */
	inline operator value_type&()
	{
		return value;
	}

	inline operator const value_type&() const
	{
		return value;
	}

	/**
	 * @brief 强制转换赋值
	 */
	template<typename _T2>
	inline ad_point<_T, _OutDim, _InDim>& operator=(const ad_point<_T2, _OutDim, _InDim>& rhs)
	{
		value = (value_type)rhs.value;
		derivative = (derivative_type)rhs.derivative;
		return *this;
	}

	/**
	 * @brief 提取分量
	 */
	inline ad_point<_T, 1, _InDim> operator[](size_t i) const
	{
		return ad_point<_T, 1, _InDim>({value[i]}, (matrix<_T, 1, _InDim> )derivative[i]);
	}

	/**
	 * @brief 提取指定坐标的梯度
	 */
	inline vector<_T, _InDim>& grad(size_t i)
	{
		return derivative[i];
	}

	inline const vector<_T, _InDim>& grad(size_t i) const
	{
		return derivative[i];
	}

	/**
	 * @brief 雅可比矩阵值
	 */
	inline matrix<_T, _OutDim, _InDim>& jacobi(size_t i)
	{
		return derivative;
	}

	inline const matrix<_T, _OutDim, _InDim>& jacobi(size_t i) const
	{
		return derivative;
	}

	/**
	 * @brief 范数
	 */
	inline ad_point<_T, 1, _InDim> norm() const;

	/**
	 * @brief 提取偏导数的值
	 * @param i 输出分量下标
	 * @param j 输入自变量下标
	 */
	inline _T& partial(size_t i, size_t j)
	{
		return derivative[i][j];
	}

	inline const _T& partial(size_t i, size_t j) const
	{
		return derivative[i][j];
	}

	inline vector<_T, _OutDim> partial(size_t j) const
	{
		return derivative.column(j);
	}

	/**
	 * Jacobi矩阵为指定列赋值
	 * 即对所有坐标的某个自变量偏导数赋值
	 */
	inline void partial_assign(size_t j, const vector<_T, _OutDim>& partial)
	{
		for(size_t i = 0; i < _OutDim; ++i)
			derivative.elem[i][j] = partial[i];
	}

	template<typename _T2>
	inline operator ad_point<_T2, _OutDim, _InDim>() const
	{
		return ad_point<_T2, _OutDim, _InDim>((vector<_T2, _OutDim> )value, (matrix<_T2, _OutDim, _InDim> )derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline ad_point<_Result, _OutDim, _InDim> operator+(const ad_point<_T2, _OutDim, _InDim>& rhs) const
	{
		// ∇(u+v)=∇u+∇v
		return ad_point<_Result, _OutDim, _InDim>(value + rhs.value, derivative + rhs.derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline ad_point<_Result, _OutDim, _InDim> operator-(const ad_point<_T2, _OutDim, _InDim>& rhs) const
	{
		// ∇(u-v)=∇u-∇v
		return ad_point<_Result, _OutDim, _InDim>(value - rhs.value, derivative - rhs.derivative);
	}

	inline ad_point<_T, _OutDim, _InDim> operator-() const
	{
		// ∇(-u)=-∇u
		return ad_point<_T, _OutDim, _InDim>(-value, -derivative);
	}

	/**
	 * @brief 对值求内积，对偏导数求按位乘法再相加
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline ad_point<_Result, 1, _InDim> operator*(const ad_point<_T2, _OutDim, _InDim>& rhs) const
	{
		ad_point<_Result, 1, _InDim> F(vector<_Result, 1> {value * rhs.value});
		for(size_t j = 0; j < _InDim; ++j)
		{
			F.derivative[0][j] = 0;
			for(size_t i = 0; i < _OutDim; ++i)
			{
				// d(uv)=udv+vdu
				F.derivative[0][j] += derivative[i][j] * rhs.value[i] + rhs.derivative[i][j] * value[i];
			}
		}
		return F;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline ad_point<_Result, _OutDim, _InDim> operator/(const ad_point<_T2, _OutDim, _InDim>& rhs) const
	{
		// ∇(u/v)=(v∇u-u∇v)/v²
		ad_point<_Result, _OutDim, _InDim> F(value / rhs.value);
		for(size_t i = 0; i < _OutDim; ++i)
			F.derivative[i] = (derivative[i] * rhs.value[i] - rhs.derivative[i] * value[i]) / (rhs.value[i] * rhs.value[i]);
		return F;
	}

	/**
	 * @brief 所有坐标相加，Jacobi矩阵也是每列相加
	 */
	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline ad_point<_Result, 1, _InDim> sum() const
	{
		ad_point<_Result, 1, _InDim> F(vector<_Result, 1> {value.sum()});
		for(size_t j = 0; j < _InDim; ++j)
		{
			F.derivative[0][j] = derivative.sum_column(j);
		}
		return F;
	}

	/**
	 * @brief 向量值切片，区间[begin, begin + _Length)
	 * 		  对应的Jacobi矩阵将从(begin, begin)->(begin + _Length, begin + _Length)分块，得到新的_Length x _Length方阵
	 */
	template<size_t _Length>
	inline ad_point<_T, _Length, _Length> slice(size_t begin) const
	{
		return ad_point<_T, _Length, _Length>(value.template slice < _Length > (begin), partition<_Length, _Length>(derivative, begin, begin));
	}
};

template<typename _T1, typename _T2, size_t _OutDim, size_t _InDim, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline ad_point<_Result, _OutDim, _InDim> operator*(const _T1& lhs, const ad_point<_T2, _OutDim, _InDim>& rhs)
{
	// ∇(ku)=k∇u
	ad_point<_Result, _OutDim, _InDim> F;
	F.value = lhs * rhs.value;
	F.derivative = lhs * rhs.derivative;
	return F;
}

template<typename _T1, typename _T2, size_t _OutDim, size_t _InDim, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline ad_point<_Result, _OutDim, _InDim> operator*(const ad_point<_T1, _OutDim, _InDim>& lhs, const _T2& rhs)
{
	// ∇(uk)=k∇u
	ad_point<_Result, _OutDim, _InDim> F;
	F.value = lhs.value * rhs;
	F.derivative = lhs.derivative * rhs;
	return F;
}

template<typename _T1, typename _T2, size_t _OutDim, size_t _InDim, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline ad_point<_Result, _OutDim, _InDim> operator/(const _T1& lhs, const ad_point<_T2, _OutDim, _InDim>& rhs)
{
	// ∇(k/u)=k∇u⁻¹=-ku⁻²∇u
	ad_point<_Result, _OutDim, _InDim> F;
	F.value = lhs / rhs.value;
	for(size_t i = 0; i < _OutDim; ++i)
		F.derivative[i] = (-lhs / (rhs.value[i] * rhs.value[i])) * rhs.derivative[i];
	return F;
}

/**
 * @brief 按位相乘，即哈达玛积
 */
template<typename _T1, typename _T2, size_t _OutDim, size_t _InDim, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline ad_point<_Result, _OutDim, _InDim> hadamard(const ad_point<_T1, _OutDim, _InDim>& lhs, const ad_point<_T2, _OutDim, _InDim>& rhs)
{
	// ∇(uv)=u∇v+v∇u
	ad_point<_Result, _OutDim, _InDim> F(hadamard(lhs.value, rhs.value));
	for(size_t i = 0; i < _OutDim; ++i)
		F.derivative[i] = lhs.derivative[i] * rhs.value[i] + rhs.derivative[i] * lhs.value[i];
	return F;
}

/**
 * @brief 拼接向量和矩阵。
 * 		  向量直接拼接，矩阵则拼接为
 * 		  [S1  0]
 * 		  [0  S2]
 * 		  其中0为0矩阵，S1、S2分别为lhs、rhs的Jacobi矩阵
 */
template<typename _T, size_t _OutDim1, size_t _InDim1, size_t _OutDim2, size_t _InDim2>
inline ad_point<_T, _OutDim1 + _OutDim2, _InDim1 + _InDim2> cat(const ad_point<_T, _OutDim1, _InDim1>& lhs, const ad_point<_T, _OutDim2, _InDim2>& rhs)
{
	ad_point<_T, _OutDim1 + _OutDim2, _InDim1 + _InDim2> F;
	F.value = cat(lhs.value, rhs.value);
	F.derivative = cat(lhs.derivative, matrix<_T, _OutDim1, _InDim2>::zero(), matrix<_T, _OutDim2, _InDim1>::zero(), rhs.derivative);
	return F;
}

template<typename _T, size_t _OutDim, size_t _InDim1, size_t _InDim2>
inline ad_point<_T, _OutDim, _InDim1 + _InDim2> cat_derivative(const ad_point<_T, _OutDim, _InDim1>& F, const matrix<_T, _OutDim, _InDim2>& d)
{
	ad_point<_T, _OutDim, _InDim1 + _InDim2> cat_F(F.value);
	cat_F.derivative = cat_column(F.derivative, d);
	return cat_F;
}

template<typename _T, size_t _OutDim, size_t _InDim1, size_t _InDim2>
inline ad_point<_T, _OutDim, _InDim1 + _InDim2> cat_derivative(const matrix<_T, _OutDim, _InDim2>& d, const ad_point<_T, _OutDim, _InDim1>& F)
{
	ad_point<_T, _OutDim, _InDim1 + _InDim2> cat_F(F.value);
	cat_F.derivative = cat_column(d, F.derivative);
	return cat_F;
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline std::ostream& operator<<(std::ostream& os, const ad_point<_T, _OutDim, _InDim>& F)
{
	os << "value:\n" << F.value << "\nderivative:\n" << F.derivative;
	return os;
}

/**
 * @brief 常量。常量关于任何自变量的偏导数永远是0。在非O0调式级别下，会触发NRVO优化
 * @param basis 该常量的值作为系数对应的基底.
 */
template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> C(size_t basis, _T val)
{
	//仅对应下标的值为val，其偏导数全为1
	ad_point<_T, _OutDim, _InDim> cp(vector<_T, _OutDim>::zero(), matrix<_T, _OutDim, _InDim>::zero());
	cp.value[basis] = val;
	return cp;
}

/**
 * @brief 自变量函数。
 * @param i 该自变量的值作为系数对应的基底下标，不同基底之间各自独立，故非basis所在行的value和derivative均为0.
 * @param j 该自变量的下标
 */
template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> x(size_t i, size_t j, _T val)
{
	//仅对应下标的值为val，其偏导数全为1
	ad_point<_T, _OutDim, _InDim> cp = C<_T, _OutDim, _InDim>(i, val);
	cp.derivative[i][j] = _T(1); //仅basis坐标对应的梯度中关于subscript的偏导数为1，其余偏导数均为0
	return cp;
}

/**
 * @brief 将一个分量乘以其基底，以和其他分量相加共同运算。
 * @param i 该自变量的值作为系数对应的基底下标。
 */
template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> ei(size_t i, const ad_point<_T, 1, _InDim>& val)
{
	ad_point<_T, _OutDim, _InDim> cp = C<_T, _OutDim, _InDim>(i, val.value[0]);
	cp.derivative[i] = val.derivative[0]; //该分量的偏导数保留
	return cp;
}

//函数都会映射所有下标的分量
template<typename _T, size_t _OutDim, size_t _InDim, typename _Derived>
struct mapping
{
	//必须声明为__attribute__((pure))以告知编译器该函数不修改任何变量，从而允许编译器优化this指针带来的开销，将其运行时调用代价降低到static函数相当
	template<typename ..._ExtraParams>
	__attribute__((const)) inline static _T value(_T x, _ExtraParams ...extra_params);

	template<typename ..._ExtraParams>
	__attribute__((pure)) inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad, _ExtraParams ...extra_params);

	template<typename ..._ExtraParams>
	__attribute__((always_inline, pure)) inline static ad_point<_T, _OutDim, _InDim> value(const ad_point<_T, _OutDim, _InDim>& x, _ExtraParams ...extra_params)
	{
		ad_point<_T, _OutDim, _InDim> result;
		for(size_t i = 0; i < _OutDim; ++i)
		{
			_T result_value = _Derived::value(x.value[i], extra_params...);
			result.value[i] = result_value;
			result.derivative[i] = _Derived::derivative(x.value[i], result_value, x.derivative[i], extra_params...);
		}
		return result;
	}
};

/**
 * @brief 以e为底数的指数函数
 */
template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> exp(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __exp_impl: public mapping<_T, _OutDim, _InDim, __exp_impl>
	{
		__attribute__((always_inline, const)) inline static _T value(_T x)
		{
			return exp(x);
		}

		__attribute__((always_inline, pure)) inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad * fx;
		}

		using mapping<_T, _OutDim, _InDim, __exp_impl >::value;
	};
	return __exp_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> ln(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __ln_impl: public mapping<_T, _OutDim, _InDim, __ln_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return ln(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / x;
		}

		using mapping<_T, _OutDim, _InDim, __ln_impl >::value;
	};
	return __ln_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> lg(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __lg_impl: public mapping<_T, _OutDim, _InDim, __lg_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return lg(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / (x * constant<_T>::ln10);
		}

		using mapping<_T, _OutDim, _InDim, __lg_impl >::value;
	};
	return __lg_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> sqrt(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __sqrt_impl: public mapping<_T, _OutDim, _InDim, __sqrt_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return sqrt(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / (_T(2) * fx);
		}

		using mapping<_T, _OutDim, _InDim, __sqrt_impl >::value;
	};
	return __sqrt_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, 1, _InDim> ad_point<_T, _OutDim, _InDim>::norm() const
{
	return sqrt((*this) * (*this));
}

//幂函数：x^n，n是常数
template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> pow(const ad_point<_T, _OutDim, _InDim>& x, const _T n)
{
	struct __pow_impl: public mapping<_T, _OutDim, _InDim, __pow_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x, const _T n)
		{
			return pow(x, n);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad, const _T n)
		{
			return grad * (n * pow(x, n - _T(1)));
		}

		using mapping<_T, _OutDim, _InDim, __pow_impl >::value;
	};
	return __pow_impl::value(x, n);
}

//指数函数：a^x，a是常数
template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> pow(const _T a, const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __pow_exp_impl: public mapping<_T, _OutDim, _InDim, __pow_exp_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x, const _T a)
		{
			return pow(a, x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad, const _T a)
		{
			return grad * ((ln(a) * fx));
		}

		using mapping<_T, _OutDim, _InDim, __pow_exp_impl >::value;
	};
	return __pow_exp_impl::value(x, a);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> abs(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __abs_impl: public mapping<_T, _OutDim, _InDim, __abs_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return abs(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad * (x >= _T(0) ? _T(1) : _T(-1));
		}

		using mapping<_T, _OutDim, _InDim, __abs_impl >::value;
	};
	return __abs_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> sin(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __sin_impl: public mapping<_T, _OutDim, _InDim, __sin_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return sin(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad * cos(x);
		}

		using mapping<_T, _OutDim, _InDim, __sin_impl >::value;
	};
	return __sin_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> cos(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __cos_impl: public mapping<_T, _OutDim, _InDim, __cos_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return cos(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad * (-sin(x));
		}

		using mapping<_T, _OutDim, _InDim, __cos_impl >::value;
	};
	return __cos_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> tan(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __tan_impl: public mapping<_T, _OutDim, _InDim, __tan_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return tan(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			_T cos_val = cos(x);
			return grad / (cos_val * cos_val);
		}

		using mapping<_T, _OutDim, _InDim, __tan_impl >::value;
	};
	return __tan_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> arcsin(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __arcsin_impl: public mapping<_T, _OutDim, _InDim, __arcsin_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arcsin(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / sqrt(_T(1) - x * x);
		}

		using mapping<_T, _OutDim, _InDim, __arcsin_impl >::value;
	};
	return __arcsin_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> arccos(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __arccos_impl: public mapping<_T, _OutDim, _InDim, __arccos_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arccos(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / sqrt(_T(1) - x * x);
		}

		using mapping<_T, _OutDim, _InDim, __arccos_impl >::value;
	};
	return __arccos_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> arctan(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __arctan_impl: public mapping<_T, _OutDim, _InDim, __arctan_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arctan(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / (_T(1) + x * x);
		}

		using mapping<_T, _OutDim, _InDim, __arctan_impl >::value;
	};
	return __arctan_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> sinh(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __sinh_impl: public mapping<_T, _OutDim, _InDim, __sinh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return sinh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad * cosh(x);
		}

		using mapping<_T, _OutDim, _InDim, __sinh_impl >::value;
	};
	return __sinh_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> cosh(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __cosh_impl: public mapping<_T, _OutDim, _InDim, __cosh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return cosh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad * sinh(x);
		}

		using mapping<_T, _OutDim, _InDim, __cosh_impl >::value;
	};
	return __cosh_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> tanh(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __tanh_impl: public mapping<_T, _OutDim, _InDim, __tanh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return tanh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad * (_T(1) - fx * fx);
		}

		using mapping<_T, _OutDim, _InDim, __tanh_impl >::value;
	};
	return __tanh_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> arcsinh(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __arcsinh_impl: public mapping<_T, _OutDim, _InDim, __arcsinh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arcsinh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / sqrt(x * x + _T(1));
		}

		using mapping<_T, _OutDim, _InDim, __arcsinh_impl >::value;
	};
	return __arcsinh_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> arccosh(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __arccosh_impl: public mapping<_T, _OutDim, _InDim, __arccosh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arccosh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / sqrt(x * x - _T(1));
		}

		using mapping<_T, _OutDim, _InDim, __arccosh_impl >::value;
	};
	return __arccosh_impl::value(x);
}

template<typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, _OutDim, _InDim> arctanh(const ad_point<_T, _OutDim, _InDim>& x)
{
	struct __arctanh_impl: public mapping<_T, _OutDim, _InDim, __arctanh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arctanh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_T, _InDim> derivative(_T x, _T fx, const vector<_T, _InDim>& grad)
		{
			return grad / (_T(1) - x * x);
		}

		using mapping<_T, _OutDim, _InDim, __arctanh_impl >::value;
	};
	return __arctanh_impl::value(x);
}

}

#endif//_MATH_AUTODIFF
