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
template<typename _T, size_t _OutOrder, size_t ... _Dims>
struct ad_point;

template<template<typename, size_t, size_t...> typename _AdTpl, typename _T, size_t _OutOrder, typename _Indexes, size_t ... _Dims>
struct adtpl_of_dims;

template<template<typename, size_t, size_t...> typename _AdTpl, typename _T, size_t _OutOrder, size_t ... _Indexes, size_t ... _Dims>
struct adtpl_of_dims<_AdTpl, _T, _OutOrder, tplmp::type_pack<tplmp::_size_t<_Indexes> ...>, _Dims...>
{
	typedef _AdTpl<_T, _OutOrder, tplmp::size_t_at<_Indexes, _Dims...>::value...> type;
};

template<template<typename, size_t, size_t...> typename _AdTpl, typename _T, size_t _OutOrder, size_t _ReservedOrders, size_t ... _Dims>
struct adtpl_of_dims_front
{
	typedef typename adtpl_of_dims<_AdTpl, _T, _OutOrder, tplmp::index_sequence_t<size_t, 0, _ReservedOrders>, _Dims...>::type type;
};

template<template<typename, size_t, size_t...> typename _AdTpl, typename _T, size_t _OutOrder, size_t _ReservedOrders, size_t ... _Dims>
struct adtpl_of_dims_back
{
	typedef typename adtpl_of_dims<_AdTpl, _T, _OutOrder, tplmp::index_sequence_t<size_t, sizeof...(_Dims) - _ReservedOrders, _ReservedOrders>, _Dims...>::type type;
};

template<typename _T, size_t _OutOrder, typename _Indexes, size_t ... _Dims>
struct ad_point_of_dims
{
	typedef typename adtpl_of_dims<ad_point, _T, _OutOrder, _Indexes, _Dims...>::type type;
};

/**
 * @brief AD输出张量升阶
 */
template<template<typename, size_t, size_t...> typename _AdTpl, typename _T, size_t _OutOrder, typename _AdOutDims, typename _AdInDims, size_t ..._OutOuterDims>
struct higher_out_order_adtpl;

template<template<typename, size_t, size_t...> typename _AdTpl, typename _T, size_t _OutOrder, size_t ... _AdOutDims, size_t ... _AdInDims, size_t ..._OutOuterDims>
struct higher_out_order_adtpl<_AdTpl, _T, _OutOrder, tplmp::type_pack<tplmp::_size_t<_AdOutDims>...>, tplmp::type_pack<tplmp::_size_t<_AdInDims>...>, _OutOuterDims...>
{
	typedef _AdTpl<_T, _OutOrder + sizeof...(_OutOuterDims), _OutOuterDims..., _AdOutDims..., _AdInDims...> type;
};

template<typename _Ad, size_t ... _OutOuterDims>
struct higher_out_order_ad_point
{
	typedef typename higher_out_order_adtpl<ad_point, typename _Ad::data_type, _Ad::out_order, typename _Ad::out_dim_t, typename _Ad::in_dim_t, _OutOuterDims...>::type type;
};

/**
 * @brief 张量减少一阶，去除最外层_RmOrders个维度
 */
template<template<typename, size_t, size_t...> typename _AdTpl, typename _T, size_t _OutOrder, typename _AdOutDims, typename _AdInDims, size_t _RmOrders>
struct lower_out_order_adtpl;

template<template<typename, size_t, size_t...> typename _AdTpl, typename _T, size_t _OutOrder, size_t ... _AdOutDims, size_t ... _AdInDims, size_t _RmOrders>
struct lower_out_order_adtpl<_AdTpl, _T, _OutOrder, tplmp::type_pack<tplmp::_size_t<_AdOutDims>...>, tplmp::type_pack<tplmp::_size_t<_AdInDims>...>, _RmOrders>
{
	//_OutOrder - _RmOrders为新的输出阶数
	typedef typename adtpl_of_dims_back<_AdTpl, _T, _OutOrder - _RmOrders, _OutOrder - _RmOrders + sizeof...(_AdInDims), _AdOutDims..., _AdInDims...>::type type;
};

template<typename _Ad, size_t _RmOrders>
struct lower_out_order_ad_point
{
	typedef typename lower_out_order_adtpl<ad_point, typename _Ad::data_type, _Ad::out_order, typename _Ad::out_dim_t, typename _Ad::in_dim_t, _RmOrders>::type type;
};

/**
 * @brief 计算中间值，包括函数值和偏导数值。
 * 表示in_dim[0]元out_dim[0]维向量函数。当out_dim[0]==1时退化为标量函数。
 * 每行都是各自独立的，只与对应下标的分量有关。
 */
template<typename _T, size_t _OutOrder, size_t ... _Dims>
struct ad_point
{
	static_assert(_OutOrder <= sizeof...(_Dims), "out order should be lower than dims size");

	typedef ad_point<_T, _OutOrder, _Dims...> type;

	typedef _T data_type;

	/**
	 * @brief 值类型
	 */
	typedef typename tensor_of_dims_front<_T, _OutOrder, _Dims...>::type value_type;

	typedef typename value_type::comp_type value_comp_type;

	static constexpr size_t out_order = value_type::order;

	static constexpr auto out_dim = value_type::dim;

	typedef typename value_type::dim_t out_dim_t;

	/**
	 * @brief 自变量类型
	 */
	typedef typename tensor_of_dims_back<_T, sizeof...(_Dims) - _OutOrder, _Dims...>::type independent_variable_type;

	typedef typename independent_variable_type::comp_type independent_variable_comp_type;

	static constexpr size_t in_order = independent_variable_type::order;

	static constexpr auto in_dim = independent_variable_type::dim;

	typedef typename independent_variable_type::dim_t in_dim_t;

	/**
	 * @brief 导数类型，先排列所有输出维度，再排列所有输入维度
	 */
	typedef tensor<_T, _Dims...> derivative_type;

	typedef typename derivative_type::comp_type derivative_comp_type;

	template<size_t ..._OutOuterDims>
	struct higher_out_order
	{
		typedef typename higher_out_order_adtpl<ad_point, data_type, out_order, out_dim_t, in_dim_t, _OutOuterDims...>::type type;
	};

	/**
	 * @brief 自动微分分量，输出维度提取分量，输入维度保持不变
	 */
	template<size_t _RmOrders, typename _DataType = _T>
	struct lower_out_order
	{
		typedef typename lower_out_order_adtpl<ad_point, data_type, out_order, out_dim_t, in_dim_t, _RmOrders>::type type;
	};

	typedef typename tplmp::if_else<_OutOrder>::resolve_t<
			typename lower_out_order<1, _T>::type,
			type
	>::type comp_type;

	value_type value; //张量函数值，_OutOrder阶张量
	derivative_type derivative; //张量函数导数，类比Jacobi矩阵

	inline ad_point() = default;

	inline ad_point(const value_type& c, const derivative_type& d) :
			value(c), derivative(d)
	{
	}

	/**
	 * @brief 包装为高1阶的、最外层维度为1的张量
	 */
	inline operator typename higher_out_order<1>::type&()
	{
		return *(typename higher_out_order<1>::type*)this;
	}

	inline operator const typename higher_out_order<1>::type&() const
	{
		return *(const typename higher_out_order<1>::type*)this;
	}

	/**
	 * @brief 若最外层维度为1，则降1维的张量
	 */
	template<size_t _OutermostDim = out_dim[0], typename = typename tplmp::if_else<_OutermostDim == 1>::def<> >
	inline operator typename lower_out_order<1>::type&()
	{
		return *(typename lower_out_order<1>::type*)this;
	}

	template<size_t _OutermostDim = out_dim[0], typename = typename tplmp::if_else<_OutermostDim == 1>::def<> >
	inline operator const typename lower_out_order<1>::type&() const
	{
		return *(const typename lower_out_order<1>::type*)this;
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
	inline type& operator=(const ad_point<_T2, _OutOrder, _Dims...>& rhs)
	{
		value = (value_type)rhs.value;
		derivative = (derivative_type)rhs.derivative;
		return *this;
	}

	/**
	 * @brief 提取分量
	 */
	inline comp_type operator[](size_t i) const
	{
		return comp_type(
				{value[i]}, (derivative_comp_type)derivative[i]);
	}

	/**
	 * @brief 提取指定坐标的梯度
	 */
	inline derivative_comp_type& grad(size_t i)
	{
		return derivative[i];
	}

	inline const derivative_comp_type& grad(size_t i) const
	{
		return derivative[i];
	}

	/**
	 * @brief 雅可比矩阵值
	 */
	inline derivative_type& jacobi()
	{
		return derivative;
	}

	inline const derivative_type& jacobi() const
	{
		return derivative;
	}

	template<typename _T2>
	inline operator ad_point<_T2, _OutOrder, _Dims...>() const
	{
		typedef ad_point<_T2, _OutOrder, _Dims...> result_type;
		return result_type((typename result_type::value_type)value, (typename result_type::derivative_type)derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline ad_point<_Result, _OutOrder, _Dims...> operator+(const ad_point<_T2, _OutOrder, _Dims...>& rhs) const
	{
		// ∇(u+v)=∇u+∇v
		return ad_point<_Result, _OutOrder, _Dims...>(value + rhs.value, derivative + rhs.derivative);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline ad_point<_Result, _OutOrder, _Dims...> operator-(const ad_point<_T2, _OutOrder, _Dims...>& rhs) const
	{
		// ∇(u-v)=∇u-∇v
		return ad_point<_Result, _OutOrder, _Dims...>(value - rhs.value, derivative - rhs.derivative);
	}

	inline ad_point<_T, _OutOrder, _Dims...> operator-() const
	{
		// ∇(-u)=-∇u
		return ad_point<_T, _OutOrder, _Dims...>(-value, -derivative);
	}

	/**
	 * @brief 对值求内积，对偏导数求按位乘法再相加
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline typename lower_out_order<_OutOrder, _Result>::type operator*(const ad_point<_T2, _OutOrder, _Dims...>& rhs) const
	{
		typename lower_out_order<_OutOrder, _Result>::type F;
		F.value = value * rhs.value;
		for(size_t j = 0; j < in_dim[0]; ++j)
		{
			F.derivative[j] = 0;
			for(size_t i = 0; i < out_dim[0]; ++i)
			{
				// ∇(uv)=u∇v+v∇u
				F.derivative[j] += derivative[i][j] * rhs.value[i] + rhs.derivative[i][j] * value[i];
			}
		}
		return F;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline ad_point<_Result, _OutOrder, _Dims...> operator/(const ad_point<_T2, _OutOrder, _Dims...>& rhs) const
	{
		// ∇(u/v)=(v∇u-u∇v)/v²
		ad_point<_Result, _OutOrder, _Dims...> F;
		F.value = {value / rhs.value};
		for(size_t i = 0; i < out_dim[0]; ++i)
			F.derivative[i] = (derivative[i] * rhs.value[i] - rhs.derivative[i] * value[i]) / (rhs.value[i] * rhs.value[i]);
		return F;
	}

	/**
	 * @brief 所有坐标相加，Jacobi矩阵也是每列相加
	 */
	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline typename lower_out_order<1, _Result>::type sum() const
	{
		typename lower_out_order<1, _Result>::type F;
		F.value = {value.sum()};
		for(size_t j = 0; j < in_dim[0]; ++j)
		{
			F.derivative[0][j] = derivative.sum_column(j); //TODO
		}
		return F;
	}
};

template<typename _T1, typename _T2, size_t _OutOrder, size_t ... _Dims, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline ad_point<_Result, _OutOrder, _Dims...> operator*(const _T1& lhs, const ad_point<_T2, _OutOrder, _Dims...>& rhs)
{
	// ∇(ku)=k∇u
	ad_point<_Result, _OutOrder, _Dims...> F;
	F.value = lhs * rhs.value;
	F.derivative = lhs * rhs.derivative;
	return F;
}

template<typename _T1, typename _T2, size_t _OutOrder, size_t ... _Dims, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline ad_point<_Result, _OutOrder, _Dims...> operator*(const ad_point<_T1, _OutOrder, _Dims...>& lhs, const _T2& rhs)
{
	// ∇(uk)=k∇u
	ad_point<_Result, _OutOrder, _Dims...> F;
	F.value = lhs.value * rhs;
	F.derivative = lhs.derivative * rhs;
	return F;
}

template<typename _T1, typename _T2, size_t _OutOrder, size_t ... _Dims, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline ad_point<_Result, _OutOrder, _Dims...> operator/(const _T1& lhs, const ad_point<_T2, _OutOrder, _Dims...>& rhs)
{
	// ∇(k/u)=k∇u⁻¹=-ku⁻²∇u
	typedef ad_point<_Result, _OutOrder, _Dims...> result_type;
	result_type F;
	F.value = lhs / rhs.value;
	for(size_t i = 0; i < result_type::out_dim[0]; ++i)
		F.derivative[i] = (-lhs / (rhs.value[i] * rhs.value[i])) * rhs.derivative[i];
	return F;
}

/**
 * @brief 按位相乘，即哈达玛积
 */
template<typename _T1, typename _T2, size_t _OutOrder, size_t ... _Dims, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline ad_point<_Result, _OutOrder, _Dims...> hadamard(const ad_point<_T1, _OutOrder, _Dims...>& lhs, const ad_point<_T2, _OutOrder, _Dims...>& rhs)
{
	// ∇(uv)=u∇v+v∇u
	typedef ad_point<_Result, _OutOrder, _Dims...> result_type;
	ad_point<_Result, _OutOrder, _Dims...> F(hadamard(lhs.value, rhs.value));
	for(size_t i = 0; i < result_type::out_dim[0]; ++i)
		F.derivative[i] = lhs.derivative[i] * rhs.value[i] + rhs.derivative[i] * lhs.value[i];
	return F;
}

/**
 * 多元向量函数专用
 */

/**
 * @brief 向量值切片，区间[begin, begin + _Length)
 * 		  对应的Jacobi矩阵将从(begin, begin)->(begin + _Length, begin + _Length)分块，得到新的_Length x _Length方阵
 */
template<size_t _Length, typename _T, size_t _OutDim, size_t _InDim>
inline ad_point<_T, 1, _Length, _Length> slice(const ad_point<_T, 1, _OutDim, _InDim>& ad, size_t begin)
{
	return ad_point<_T, 1, _Length, _Length>(ad.value.template slice<_Length>(begin), partition<_Length, _Length>(ad.derivative, begin, begin));
}

/**
 * @brief 拼接向量和矩阵。
 * 		  向量直接拼接，矩阵则拼接为
 * 		  [S1  0]
 * 		  [0  S2]
 * 		  其中0为0矩阵，S1、S2分别为lhs、rhs的Jacobi矩阵
 */
template<typename _T, size_t _OutDim1, size_t _InDim1, size_t _OutDim2, size_t _InDim2>
inline ad_point<_T, 1, _OutDim1 + _OutDim2, _InDim1 + _InDim2> cat(const ad_point<_T, 1, _OutDim1, _InDim1>& lhs, const ad_point<_T, 1, _OutDim2, _InDim2>& rhs)
{
	ad_point<_T, 1, _OutDim1 + _OutDim2, _InDim1 + _InDim2> F;
	F.value = cat(lhs.value, rhs.value);
	F.derivative = cat(lhs.derivative, matrix<_T, _OutDim1, _InDim2>::zero(), matrix<_T, _OutDim2, _InDim1>::zero(), rhs.derivative);
	return F;
}

template<typename _T, size_t _OutDim, size_t _InDim1, size_t _InDim2>
inline ad_point<_T, 1, _OutDim, _InDim1 + _InDim2> cat_derivative(const ad_point<_T, 1, _OutDim, _InDim1>& F, const tensor<_T, _OutDim, _InDim2>& d)
{
	ad_point<_T, 1, _OutDim, _InDim1 + _InDim2> cat_F;
	cat_F.value = F.value;
	cat_F.derivative = cat_column(F.derivative, d);
	return cat_F;
}

template<typename _T, size_t _InDim1, size_t _InDim2>
inline ad_point<_T, 0, _InDim1 + _InDim2> cat_derivative(const ad_point<_T, 0, _InDim1>& F, const tensor<_T, 1, _InDim2>& d)
{
	return (ad_point<_T, 0, _InDim1 + _InDim2> )cat_derivative((const ad_point<_T, 1, 1, _InDim1>&)F, d);
}

template<typename _T, size_t _OutDim, size_t _InDim1, size_t _InDim2>
inline ad_point<_T, 1, _OutDim, _InDim1 + _InDim2> cat_derivative(const tensor<_T, _OutDim, _InDim1>& d, const ad_point<_T, 1, _OutDim, _InDim2>& F)
{
	ad_point<_T, 1, _OutDim, _InDim1 + _InDim2> cat_F;
	cat_F.value = F.value;
	cat_F.derivative = cat_column(d, F.derivative);
	return cat_F;
}

template<typename _T, size_t _InDim1, size_t _InDim2>
inline ad_point<_T, 0, _InDim1 + _InDim2> cat_derivative(const tensor<_T, 1, _InDim1>& d, const ad_point<_T, 0, _InDim2>& F)
{
	return (ad_point<_T, 0, _InDim1 + _InDim2> )cat_derivative(d, (const ad_point<_T, 1, 1, _InDim1>&)F);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline std::ostream& operator<<(std::ostream& os, const ad_point<_T, _OutOrder, _Dims...>& F)
{
	os << "value:\n" << F.value << "\nderivative:\n" << F.derivative;
	return os;
}

/**
 * @brief 常量。常量关于任何自变量的偏导数永远是0。在非O0调式级别下，会触发NRVO优化
 * @param basis 该常量的值作为系数对应的基底.
 */
template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> C(size_t basis, const typename ad_point<_T, _OutOrder, _Dims...>::value_comp_type& val)
{
	//仅对应下标的值为val，其偏导数全为1
	typedef ad_point<_T, _OutOrder, _Dims...> result_type;
	result_type ad(result_type::value_type::zero(), result_type::derivative_type::zero());
	ad.value[basis] = val;
	return ad;
}

/**
 * @brief 自变量函数。
 * @param i 该自变量当前参与计算的输出基底下标，不同基底之间各自独立，故非basis所在行的value和derivative均为0.
 * @param j 该自变量的下标
 */
template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> x(size_t i, size_t j, const typename ad_point<_T, _OutOrder, _Dims...>::value_comp_type& val)
{
	//仅对应下标的值为val，其偏导数全为1
	ad_point<_T, _OutOrder, _Dims...> ad = C<_T, _OutOrder, _Dims...>(i, val);
	ad.derivative[i][j] = identity<_T, mul>(); //仅basis坐标对应的梯度中关于subscript的偏导数为1，其余偏导数均为0
	return ad;
}

/**
 * @brief 将一个分量乘以其基底，以和其他分量相加共同运算。
 * @param i 该自变量的值作为系数对应的基底下标。
 */
template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> ei(size_t i, const typename ad_point<_T, _OutOrder, _Dims...>::value_comp_type& val)
{
	ad_point<_T, _OutOrder, _Dims...> ad = C<_T, _OutOrder, _Dims...>(i, val.value[0]);
	ad.derivative[i] = val.derivative[0]; //该分量的偏导数保留
	return ad;
}

//函数都会映射所有下标的分量
template<typename _Derived, typename _T, size_t _OutOrder, size_t ... _Dims>
struct mapping
{
	typedef ad_point<_T, _OutOrder, _Dims...> ad_type;
	typedef typename ad_type::value_type value_type;
	typedef typename ad_type::independent_variable_type independent_variable_type;
	typedef typename ad_type::derivative_type derivative_type;
	typedef typename ad_type::derivative_comp_type derivative_comp_type;

	//必须声明为__attribute__((pure))以告知编译器该函数不修改任何变量，从而允许编译器优化this指针带来的开销，将其运行时调用代价降低到static函数相当
	template<typename ..._ExtraParams>
	__attribute__((const)) inline static _T value(const _T& x, _ExtraParams ...extra_params);

	template<typename ..._ExtraParams>
	__attribute__((pure)) inline static derivative_comp_type derivative(const _T& x, const _T& fx, const derivative_comp_type& grad, _ExtraParams ...extra_params);

	template<typename ..._ExtraParams>
	__attribute__((always_inline, pure)) inline static ad_type value(const ad_type& x, _ExtraParams ...extra_params)
	{
		ad_type result;
		for(size_t i = 0; i < ad_type::out_dim[0]; ++i)
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
template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> exp(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __exp_impl;
	typedef mapping<__exp_impl, _T, _OutOrder, _Dims...> base;
	struct __exp_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return exp(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad * fx;
		}

		using base::value;
	};
	return __exp_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> ln(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __ln_impl;
	typedef mapping<__ln_impl, _T, _OutOrder, _Dims...> base;
	struct __ln_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return ln(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / x;
		}

		using base::value;
	};
	return __ln_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> lg(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __lg_impl;
	typedef mapping<__lg_impl, _T, _OutOrder, _Dims...> base;
	struct __lg_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return lg(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / (x * constant<_T>::ln10);
		}

		using base::value;
	};
	return __lg_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> sqrt(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __sqrt_impl;
	typedef mapping<__sqrt_impl, _T, _OutOrder, _Dims...> base;
	struct __sqrt_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return sqrt(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / (_T(2) * fx);
		}

		using base::value;
	};
	return __sqrt_impl::value(x);
}

/**
 * @brief Frobenius范数
 */
template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline auto frobenius_norm(const ad_point<_T, _OutOrder, _Dims...>& T) -> decltype(sqrt(T * T))
{
	return sqrt(T * T);
}

//幂函数：x^n，n是常数
template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> pow(const ad_point<_T, _OutOrder, _Dims...>& x, const _T n)
{
	struct __pow_impl;
	typedef mapping<__pow_impl, _T, _OutOrder, _Dims...> base;
	struct __pow_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x, const _T& n)
		{
			return pow(x, n);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad, const _T& n)
		{
			return grad * (n * pow(x, n - _T(1)));
		}

		using base::value;
	};
	return __pow_impl::value(x, n);
}

//指数函数：a^x，a是常数
template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> pow(const _T a, const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __pow_exp_impl;
	typedef mapping<__pow_exp_impl, _T, _OutOrder, _Dims...> base;
	struct __pow_exp_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x, const _T& a)
		{
			return pow(a, x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad, const _T& a)
		{
			return grad * (ln(a) * fx);
		}

		using base::value;
	};
	return __pow_exp_impl::value(x, a);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> abs(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __abs_impl;
	typedef mapping<__abs_impl, _T, _OutOrder, _Dims...> base;
	struct __abs_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return abs(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad * (x >= _T(0) ? _T(1) : _T(-1));
		}

		using base::value;
	};
	return __abs_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> sin(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __sin_impl;
	typedef mapping<__sin_impl, _T, _OutOrder, _Dims...> base;
	struct __sin_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return sin(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad * cos(x);
		}

		using base::value;
	};
	return __sin_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> cos(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __cos_impl;
	typedef mapping<__cos_impl, _T, _OutOrder, _Dims...> base;
	struct __cos_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return cos(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad * (-sin(x));
		}

		using base::value;
	};
	return __cos_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> tan(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __tan_impl;
	typedef mapping<__tan_impl, _T, _OutOrder, _Dims...> base;
	struct __tan_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return tan(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			_T cos_val = cos(x);
			return grad / (cos_val * cos_val);
		}

		using base::value;
	};
	return __tan_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> arcsin(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __arcsin_impl;
	typedef mapping<__arcsin_impl, _T, _OutOrder, _Dims...> base;
	struct __arcsin_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return arcsin(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / sqrt(_T(1) - x * x);
		}

		using base::value;
	};
	return __arcsin_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> arccos(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __arccos_impl;
	typedef mapping<__arccos_impl, _T, _OutOrder, _Dims...> base;
	struct __arccos_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return arccos(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / sqrt(_T(1) - x * x);
		}

		using base::value;
	};
	return __arccos_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> arctan(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __arctan_impl;
	typedef mapping<__arctan_impl, _T, _OutOrder, _Dims...> base;
	struct __arctan_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return arctan(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / (_T(1) + x * x);
		}

		using base::value;
	};
	return __arctan_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> sinh(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __sinh_impl;
	typedef mapping<__sinh_impl, _T, _OutOrder, _Dims...> base;
	struct __sinh_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return sinh(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad * cosh(x);
		}

		using base::value;
	};
	return __sinh_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> cosh(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __cosh_impl;
	typedef mapping<__cosh_impl, _T, _OutOrder, _Dims...> base;
	struct __cosh_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return cosh(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad * sinh(x);
		}

		using base::value;
	};
	return __cosh_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> tanh(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __tanh_impl;
	typedef mapping<__tanh_impl, _T, _OutOrder, _Dims...> base;
	struct __tanh_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return tanh(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad * (_T(1) - fx * fx);
		}

		using base::value;
	};
	return __tanh_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> arcsinh(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __arcsinh_impl;
	typedef mapping<__arcsinh_impl, _T, _OutOrder, _Dims...> base;
	struct __arcsinh_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return arcsinh(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / sqrt(x * x + _T(1));
		}

		using base::value;
	};
	return __arcsinh_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> arccosh(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __arccosh_impl;
	typedef mapping<__arccosh_impl, _T, _OutOrder, _Dims...> base;
	struct __arccosh_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return arccosh(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / sqrt(x * x - _T(1));
		}

		using base::value;
	};
	return __arccosh_impl::value(x);
}

template<typename _T, size_t _OutOrder, size_t ... _Dims>
inline ad_point<_T, _OutOrder, _Dims...> arctanh(const ad_point<_T, _OutOrder, _Dims...>& x)
{
	struct __arctanh_impl;
	typedef mapping<__arctanh_impl, _T, _OutOrder, _Dims...> base;
	struct __arctanh_impl: base
	{
		__attribute__((always_inline, const)) inline static _T value(const _T& x)
		{
			return arctanh(x);
		}

		__attribute__((always_inline, pure)) inline static typename base::derivative_comp_type derivative(const _T& x, const _T& fx, const typename base::derivative_comp_type& grad)
		{
			return grad / (_T(1) - x * x);
		}

		using base::value;
	};
	return __arctanh_impl::value(x);
}

}

#endif//_MATH_AUTODIFF
