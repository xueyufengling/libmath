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
/**
 * @brief 计算中间值，包括函数值和偏导数值。
 * 表示_InDim元_OutDim维向量函数。当_OutDim==1时退化为标量函数。
 * 每行都是各自独立的，只与对应下标的分量有关。
 */
template<size_t _OutDim, size_t _InDim, typename _T>
struct ad_point
{
	vector<_OutDim, _T> value; //多元函数值，即坐标，列向量
	matrix<_OutDim, _InDim, _T> derivative; //多元向量函数导数，即Jacobi矩阵，每个行向量都是对应坐标的梯度

	inline ad_point() = default;

	inline ad_point(const vector<_OutDim, _T>& c) :
			value(c)
	{
	}

	inline ad_point(const matrix<_OutDim, _InDim, _T>& d) :
			derivative(d)
	{
	}

	inline ad_point(const vector<_OutDim, _T>& c, const matrix<_OutDim, _InDim, _T>& d) :
			value(c), derivative(d)
	{
	}

	/**
	 * @brief 隐式转换到坐标
	 */
	inline operator vector<_OutDim, _T>&()
	{
		return value;
	}

	inline operator const vector<_OutDim, _T>&() const
	{
		return value;
	}

	/**
	 * @brief 强制转换赋值
	 */
	template<typename _T2>
	inline ad_point<_OutDim, _InDim, _T>& operator=(const ad_point<_OutDim, _InDim, _T2>& other)
	{
		value = (vector<_OutDim, _T> )other.value;
		derivative = (matrix<_OutDim, _InDim, _T> )other.derivative;
		return *this;
	}

	/**
	 * @brief 提取分量
	 */
	inline ad_point<1, _InDim, _T> operator[](size_t i) const
	{
		return ad_point<1, _InDim, _T>({value[i]}, (matrix<1, _InDim, _T> )derivative[i]);
	}

	/**
	 * @brief 提取指定坐标的梯度
	 */
	inline vector<_InDim, _T>& grad(size_t i)
	{
		return derivative[i];
	}

	inline const vector<_InDim, _T>& grad(size_t i) const
	{
		return derivative[i];
	}

	/**
	 * @brief 雅可比矩阵值
	 */
	inline matrix<_OutDim, _InDim, _T>& jacobi(size_t i)
	{
		return derivative;
	}

	inline const matrix<_OutDim, _InDim, _T>& jacobi(size_t i) const
	{
		return derivative;
	}

	/**
	 * @brief 范数
	 */
	inline ad_point<1, _InDim, _T> norm() const;

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

	inline vector<_OutDim, _T> partial(size_t j) const
	{
		return derivative.column(j);
	}

	/**
	 * Jacobi矩阵为指定列赋值
	 * 即对所有坐标的某个自变量偏导数赋值
	 */
	inline void partial_assign(size_t j, const vector<_OutDim, _T>& partial)
	{
		for(size_t i = 0; i < _OutDim; ++i)
			derivative.elem[i][j] = partial[i];
	}

	template<typename _T2>
	inline operator ad_point<_OutDim, _InDim, _T2>() const
	{
		return ad_point<_OutDim, _InDim, _T2>((vector<_OutDim, _T2> )value, (matrix<_OutDim, _InDim, _T2> )derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline ad_point<_OutDim, _InDim, _Result> operator+(const ad_point<_OutDim, _InDim, _T2>& other) const
	{
		// ∇(u+v)=∇u+∇v
		return ad_point<_OutDim, _InDim, _Result>(value + other.value, derivative + other.derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline ad_point<_OutDim, _InDim, _Result> operator-(const ad_point<_OutDim, _InDim, _T2>& other) const
	{
		// ∇(u-v)=∇u-∇v
		return ad_point<_OutDim, _InDim, _Result>(value - other.value, derivative - other.derivative);
	}

	inline ad_point<_OutDim, _InDim, _T> operator-() const
	{
		// ∇(-u)=-∇u
		return ad_point<_OutDim, _InDim, _T>(-value, -derivative);
	}

	/**
	 * @brief 点乘。比hadamard(adp, adp).sum()更高效
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline ad_point<1, _InDim, _Result> operator*(const ad_point<_OutDim, _InDim, _T2>& other) const
	{
		ad_point<1, _InDim, _Result> adp(vector<1, _Result> {value * other.value});
		for(size_t j = 0; j < _InDim; ++j)
		{
			adp.derivative[0][j] = 0;
			for(size_t i = 0; i < _OutDim; ++i)
			{
				// d(uv)=udv+vdu
				adp.derivative[0][j] += derivative[i][j] * other.value[i] + other.derivative[i][j] * value[i];
			}
		}
		return adp;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline ad_point<_OutDim, _InDim, _Result> operator*(_T2 t) const
	{
		// ∇(ku)=k∇u
		ad_point<_OutDim, _InDim, _Result> adp = *this;
		adp.value *= t;
		adp.derivative *= t;
		return adp;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline ad_point<_OutDim, _InDim, _Result> operator/(const ad_point<_OutDim, _InDim, _T2>& other) const
	{
		// ∇(u/v)=(v∇u-u∇v)/v²
		ad_point<_OutDim, _InDim, _Result> adp(value / other.value);
		for(size_t i = 0; i < _OutDim; ++i)
			adp.derivative[i] = (derivative[i] * other.value[i] - other.derivative[i] * value[i]) / (other.value[i] * other.value[i]);
		return adp;
	}

	/**
	 * @brief 所有坐标相加，Jacobi矩阵也是每列相加
	 */
	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline ad_point<1, _InDim, _Result> sum() const
	{
		ad_point<1, _InDim, _Result> adp(vector<1, _Result> {value.sum()});
		for(size_t j = 0; j < _InDim; ++j)
		{
			adp.derivative[0][j] = derivative.sum_column(j);
		}
		return adp;
	}

	/**
	 * @brief 向量值切片，区间[begin, begin + _Length)
	 * 		  对应的Jacobi矩阵将从(begin, begin)->(begin + _Length, begin + _Length)分块，得到新的_Length x _Length方阵
	 */
	template<size_t _Length>
	inline ad_point<_Length, _Length, _T> slice(size_t begin) const
	{
		return ad_point<_Length, _Length, _T>(value.template slice<_Length>(begin), derivative.template partition<_Length, _Length>(begin, begin));
	}

	template<typename _Result = _T>
	inline ad_point<_OutDim, _InDim, _Result> mask(size_t begin, size_t end, _T mask_value) const
	{
		return ad_point<_OutDim, _InDim, _Result>(value.template mask<_Result>(begin, end, mask_value), derivative.template mask<_Result>(begin, end, 0, _InDim, mask_value));
	}
};

template<typename _T1, size_t _OutDim, size_t _InDim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline ad_point<_OutDim, _InDim, _Result> operator*(_T1 t, const ad_point<_OutDim, _InDim, _T2>& adp)
{
	// ∇(ku)=k∇u
	ad_point<_OutDim, _InDim, _Result> mul_result;
	mul_result.value = t * adp.value;
	mul_result.derivative = t * adp.derivative;
	return mul_result;
}

template<typename _T1, size_t _OutDim, size_t _InDim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline ad_point<_OutDim, _InDim, _Result> operator/(_T1 t, const ad_point<_OutDim, _InDim, _T2>& adp)
{
	// ∇(k/u)=k∇u⁻¹=-ku⁻²∇u
	ad_point<_OutDim, _InDim, _Result> div_result;
	div_result.value = t / adp.value;
	for(size_t i = 0; i < _OutDim; ++i)
		div_result.derivative[i] = (-t / (adp.value[i] * adp.value[i])) * adp.derivative[i];
	return div_result;
}

/**
 * @brief 按位相乘，即哈达玛积
 */
template<size_t _OutDim, size_t _InDim, typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline ad_point<_OutDim, _InDim, _Result> hadamard(const ad_point<_OutDim, _InDim, _T1>& adp1, const ad_point<_OutDim, _InDim, _T2>& adp2)
{
	// ∇(uv)=u∇v+v∇u
	ad_point<_OutDim, _InDim, _Result> adp(hadamard(adp1.value, adp2.value));
	for(size_t i = 0; i < _OutDim; ++i)
		adp.derivative[i] = adp1.derivative[i] * adp2.value[i] + adp2.derivative[i] * adp1.value[i];
	return adp;
}

/**
 * @brief 拼接向量和矩阵。
 * 		  向量直接拼接，矩阵则拼接为
 * 		  [S1  0]
 * 		  [0  S2]
 * 		  其中0为0矩阵，S1、S2分别为adp1、adp2的Jacobi矩阵
 */
template<size_t _OutDim1, size_t _InDim1, size_t _OutDim2, size_t _InDim2, typename _T>
inline ad_point<_OutDim1 + _OutDim2, _InDim1 + _InDim2, _T> cat(const ad_point<_OutDim1, _InDim1, _T>& adp1, const ad_point<_OutDim2, _InDim2, _T>& adp2)
{
	ad_point<_OutDim1 + _OutDim2, _InDim1 + _InDim2, _T> adp;
	adp.value = cat(adp1.value, adp2.value);
	adp.derivative = cat(adp1.derivative, matrix<_OutDim1, _InDim2, _T>::zero(), matrix<_OutDim2, _InDim1, _T>::zero(), adp2.derivative);
	return adp;
}

template<size_t _OutDim, size_t _InDim1, size_t _InDim2, typename _T>
inline ad_point<_OutDim, _InDim1 + _InDim2, _T> cat_derivative(const ad_point<_OutDim, _InDim1, _T>& adp, const matrix<_OutDim, _InDim2, _T>& d)
{
	ad_point<_OutDim, _InDim1 + _InDim2, _T> cat_adp(adp.value);
	cat_adp.derivative = cat_column(adp.derivative, d);
	return cat_adp;
}

template<size_t _OutDim, size_t _InDim1, size_t _InDim2, typename _T>
inline ad_point<_OutDim, _InDim1 + _InDim2, _T> cat_derivative(const matrix<_OutDim, _InDim2, _T>& d, const ad_point<_OutDim, _InDim1, _T>& adp)
{
	ad_point<_OutDim, _InDim1 + _InDim2, _T> cat_adp(adp.value);
	cat_adp.derivative = cat_column(d, adp.derivative);
	return cat_adp;
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline std::ostream& operator<<(std::ostream& os, const ad_point<_OutDim, _InDim, _T>& adp)
{
	os << "value:\n" << adp.value << "\nderivative:\n" << adp.derivative;
	return os;
}

/**
 * @brief 常量。常量关于任何自变量的偏导数永远是0。在非O0调式级别下，会触发NRVO优化
 * @param basis 该常量的值作为系数对应的基底.
 */
template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> C(size_t basis, _T val)
{
	//仅对应下标的值为val，其偏导数全为1
	ad_point<_OutDim, _InDim, _T> cp(vector<_OutDim, _T>::zero(), matrix<_OutDim, _InDim, _T>::zero());
	cp.value[basis] = val;
	return cp;
}

/**
 * @brief 自变量函数。
 * @param i 该自变量的值作为系数对应的基底下标，不同基底之间各自独立，故非basis所在行的value和derivative均为0.
 * @param j 该自变量的下标
 */
template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> x(size_t i, size_t j, _T val)
{
	//仅对应下标的值为val，其偏导数全为1
	ad_point<_OutDim, _InDim, _T> cp = C<_OutDim, _InDim, _T>(i, val);
	cp.derivative[i][j] = _T(1); //仅basis坐标对应的梯度中关于subscript的偏导数为1，其余偏导数均为0
	return cp;
}

/**
 * @brief 将一个分量乘以其基底，以和其他分量相加共同运算。
 * @param i 该自变量的值作为系数对应的基底下标。
 */
template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> ei(size_t i, const ad_point<1, _InDim, _T>& val)
{
	ad_point<_OutDim, _InDim, _T> cp = C<_OutDim, _InDim, _T>(i, val.value[0]);
	cp.derivative[i] = val.derivative[0]; //该分量的偏导数保留
	return cp;
}

//函数都会映射所有下标的分量
template<size_t _OutDim, size_t _InDim, typename _T, typename _Derived>
struct mapping
{
	//必须声明为__attribute__((pure))以告知编译器该函数不修改任何变量，从而允许编译器优化this指针带来的开销，将其运行时调用代价降低到static函数相当
	template<typename ..._ExtraParams>
	__attribute__((const)) inline static _T value(_T x, _ExtraParams ...extra_params);

	template<typename ..._ExtraParams>
	__attribute__((pure)) inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad, _ExtraParams ...extra_params);

	template<typename ..._ExtraParams>
	__attribute__((always_inline, pure)) inline static ad_point<_OutDim, _InDim, _T> value(const ad_point<_OutDim, _InDim, _T>& x, _ExtraParams ...extra_params)
	{
		ad_point<_OutDim, _InDim, _T> result;
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
template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> exp(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __exp_impl: public mapping<_OutDim, _InDim, _T, __exp_impl>
	{
		__attribute__((always_inline, const)) inline static _T value(_T x)
		{
			return exp(x);
		}

		__attribute__((always_inline, pure)) inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad * fx;
		}

		using mapping<_OutDim, _InDim, _T, __exp_impl >::value;
	};
	return __exp_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> ln(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __ln_impl: public mapping<_OutDim, _InDim, _T, __ln_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return ln(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / x;
		}

		using mapping<_OutDim, _InDim, _T, __ln_impl >::value;
	};
	return __ln_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> lg(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __lg_impl: public mapping<_OutDim, _InDim, _T, __lg_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return lg(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / (x * constant<_T>::ln10);
		}

		using mapping<_OutDim, _InDim, _T, __lg_impl >::value;
	};
	return __lg_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> sqrt(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __sqrt_impl: public mapping<_OutDim, _InDim, _T, __sqrt_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return sqrt(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / (_T(2) * fx);
		}

		using mapping<_OutDim, _InDim, _T, __sqrt_impl >::value;
	};
	return __sqrt_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<1, _InDim, _T> ad_point<_OutDim, _InDim, _T>::norm() const
{
	return sqrt((*this) * (*this));
}

//幂函数：x^n，n是常数
template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> pow(const ad_point<_OutDim, _InDim, _T>& x, const _T n)
{
	struct __pow_impl: public mapping<_OutDim, _InDim, _T, __pow_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x, const _T n)
		{
			return pow(x, n);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad, const _T n)
		{
			return grad * (n * pow(x, n - _T(1)));
		}

		using mapping<_OutDim, _InDim, _T, __pow_impl >::value;
	};
	return __pow_impl::value(x, n);
}

//指数函数：a^x，a是常数
template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> pow(const _T a, const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __pow_exp_impl: public mapping<_OutDim, _InDim, _T, __pow_exp_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x, const _T a)
		{
			return pow(a, x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad, const _T a)
		{
			return grad * ((ln(a) * fx));
		}

		using mapping<_OutDim, _InDim, _T, __pow_exp_impl >::value;
	};
	return __pow_exp_impl::value(x, a);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> abs(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __abs_impl: public mapping<_OutDim, _InDim, _T, __abs_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return abs(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad * (x >= _T(0) ? _T(1) : _T(-1));
		}

		using mapping<_OutDim, _InDim, _T, __abs_impl >::value;
	};
	return __abs_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> sin(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __sin_impl: public mapping<_OutDim, _InDim, _T, __sin_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return sin(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad * cos(x);
		}

		using mapping<_OutDim, _InDim, _T, __sin_impl >::value;
	};
	return __sin_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> cos(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __cos_impl: public mapping<_OutDim, _InDim, _T, __cos_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return cos(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad * (-sin(x));
		}

		using mapping<_OutDim, _InDim, _T, __cos_impl >::value;
	};
	return __cos_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> tan(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __tan_impl: public mapping<_OutDim, _InDim, _T, __tan_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return tan(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			_T cos_val = cos(x);
			return grad / (cos_val * cos_val);
		}

		using mapping<_OutDim, _InDim, _T, __tan_impl >::value;
	};
	return __tan_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> arcsin(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __arcsin_impl: public mapping<_OutDim, _InDim, _T, __arcsin_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arcsin(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / sqrt(_T(1) - x * x);
		}

		using mapping<_OutDim, _InDim, _T, __arcsin_impl >::value;
	};
	return __arcsin_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> arccos(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __arccos_impl: public mapping<_OutDim, _InDim, _T, __arccos_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arccos(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / sqrt(_T(1) - x * x);
		}

		using mapping<_OutDim, _InDim, _T, __arccos_impl >::value;
	};
	return __arccos_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> arctan(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __arctan_impl: public mapping<_OutDim, _InDim, _T, __arctan_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arctan(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / (_T(1) + x * x);
		}

		using mapping<_OutDim, _InDim, _T, __arctan_impl >::value;
	};
	return __arctan_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> sinh(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __sinh_impl: public mapping<_OutDim, _InDim, _T, __sinh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return sinh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad * cosh(x);
		}

		using mapping<_OutDim, _InDim, _T, __sinh_impl >::value;
	};
	return __sinh_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> cosh(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __cosh_impl: public mapping<_OutDim, _InDim, _T, __cosh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return cosh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad * sinh(x);
		}

		using mapping<_OutDim, _InDim, _T, __cosh_impl >::value;
	};
	return __cosh_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> tanh(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __tanh_impl: public mapping<_OutDim, _InDim, _T, __tanh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return tanh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad * (_T(1) - fx * fx);
		}

		using mapping<_OutDim, _InDim, _T, __tanh_impl >::value;
	};
	return __tanh_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> arcsinh(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __arcsinh_impl: public mapping<_OutDim, _InDim, _T, __arcsinh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arcsinh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / sqrt(x * x + _T(1));
		}

		using mapping<_OutDim, _InDim, _T, __arcsinh_impl >::value;
	};
	return __arcsinh_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> arccosh(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __arccosh_impl: public mapping<_OutDim, _InDim, _T, __arccosh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arccosh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / sqrt(x * x - _T(1));
		}

		using mapping<_OutDim, _InDim, _T, __arccosh_impl >::value;
	};
	return __arccosh_impl::value(x);
}

template<size_t _OutDim, size_t _InDim, typename _T>
inline ad_point<_OutDim, _InDim, _T> arctanh(const ad_point<_OutDim, _InDim, _T>& x)
{
	struct __arctanh_impl: public mapping<_OutDim, _InDim, _T, __arctanh_impl>
	{
		__attribute__((always_inline, const))
		inline static _T value(_T x)
		{
			return arctanh(x);
		}

		__attribute__((always_inline, pure))
		inline static vector<_InDim, _T> derivative(_T x, _T fx, const vector<_InDim, _T>& grad)
		{
			return grad / (_T(1) - x * x);
		}

		using mapping<_OutDim, _InDim, _T, __arctanh_impl >::value;
	};
	return __arctanh_impl::value(x);
}

}

#endif//_MATH_AUTODIFF
