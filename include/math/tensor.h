#ifndef _MATH_TENSOR
#define _MATH_TENSOR

#include <tplmp/base.h>

#include <sstream>

#include <math/algebra.h>

#define __tensor_loop_unroll__() __loop_unroll__(9)

namespace math
{
template<typename _T, size_t ..._Dims>
struct tensor;

/**
 * 从_Dims中取指定索引的维度组成类tensor类型，_TensorTpl是模板参数与tensor一致的模板
 */
template<template<typename, size_t...> typename _TensorTpl, typename _T, typename _Indexes, size_t ... _Dims>
struct tensortpl_of_dims;

template<template<typename, size_t...> typename _TensorTpl, typename _T, size_t ... _Indexes, size_t ... _Dims>
struct tensortpl_of_dims<_TensorTpl, _T, tplmp::type_pack<tplmp::_constexpr<size_t, _Indexes> ...>, _Dims...>
{
	typedef _TensorTpl<_T, tplmp::size_t_at<_Indexes, _Dims...>::value...> type;
};

template<template<typename, size_t...> typename _TensorTpl, typename _T, size_t _Order, size_t ... _Dims>
struct tensortpl_of_dims_front
{
	typedef typename tensortpl_of_dims<_TensorTpl, _T, tplmp::index_sequence_t<size_t, 0, _Order>, _Dims...>::type type;
};

template<template<typename, size_t...> typename _TensorTpl, typename _T, size_t _Order, size_t ... _Dims>
struct tensortpl_of_dims_back
{
	typedef typename tensortpl_of_dims<_TensorTpl, _T, tplmp::index_sequence_t<size_t, sizeof...(_Dims) - _Order, _Order>, _Dims...>::type type;
};

template<typename _T, typename _Seq, size_t ... _Dims>
struct tensor_of_dims
{
	typedef typename tensortpl_of_dims<tensor, _T, _Seq, _Dims...>::type type;
};

template<typename _T, size_t _Order, size_t ... _Dims>
struct tensor_of_dims_front
{
	typedef typename tensortpl_of_dims_front<tensor, _T, _Order, _Dims...>::type type;
};

template<typename _T, size_t _Order, size_t ... _Dims>
struct tensor_of_dims_back
{
	typedef typename tensortpl_of_dims_back<tensor, _T, _Order, _Dims...>::type type;
};

/**
 * @brief 张量升阶
 * @param _OuterDims 新增的最外层维度
 */
template<template<typename, size_t...> typename _TensorTpl, typename _T, typename _TensorDims, size_t ..._OuterDims>
struct higher_order_tensortpl;

template<template<typename, size_t...> typename _TensorTpl, typename _T, size_t ... _TensorDims, size_t ..._OuterDims>
struct higher_order_tensortpl <_TensorTpl, _T, tplmp::type_pack<tplmp::_size_t<_TensorDims>...>, _OuterDims...>
{
	typedef _TensorTpl<_T, _OuterDims..., _TensorDims...> type;
};

template<typename _Tensor, size_t ..._OuterDims>
struct higher_order_tensor
{
	typedef typename higher_order_tensortpl<tensor, typename _Tensor::data_type, typename _Tensor::dim_t, _OuterDims...>::type type;
};

/**
 * @brief 张量降阶，去除最外层_RmOrders个维度
 */
template<template<typename, size_t...> typename _TensorTpl, typename _T, typename _TensorDims, size_t _RmOrders>
struct lower_order_tensortpl;

template<template<typename, size_t...> typename _TensorTpl, typename _T, size_t... _TensorDims, size_t _RmOrders>
struct lower_order_tensortpl<_TensorTpl, _T, tplmp::type_pack<tplmp::_size_t<_TensorDims>...>, _RmOrders>
{
	typedef typename tensortpl_of_dims_back<_TensorTpl, _T, sizeof...(_TensorDims) - _RmOrders, _TensorDims...>::type type;
};

template<typename _Tensor, size_t _RmOrders>
struct lower_order_tensor
{
	typedef typename lower_order_tensortpl<tensor, typename _Tensor::data_type, typename _Tensor::dim_t, _RmOrders>::type type;
};

template<typename _T>
using scalar = tensor<_T>;

template<typename _T, size_t _Dim>
using vector = tensor<_T, _Dim>;

template<typename _T, size_t _Row, size_t _Column>
using matrix = tensor<_T, _Row, _Column>;

/**
 * @brief 0阶张量，即标量
 */
template<typename _T>
struct tensor<_T>
{
	typedef tensor<_T> type;
	typedef _T data_type;
	typedef _T comp_type;

	static constexpr size_t outmost_dim = 0;

	static constexpr size_t order = 0;
	static const size_t dim[order];

	typedef tplmp::type_pack<> dim_t;

	tensor() = default;

	inline tensor(_T s) :
			comp(s)
	{
	}

	_T comp;

	inline operator _T&()
	{
		return comp;
	}

	inline operator const _T&() const
	{
		return comp;
	}

	template<typename _T2>
	inline operator scalar<_T2>() const
	{
		return
		{	(_T2)(comp)};
	}

	template<typename _T2>
	inline operator _T2() const
	{
		return (_T2)(comp);
	}

	inline comp_type& operator[](size_t i)
	{
		return comp;
	}

	inline const comp_type& operator[](size_t i) const
	{
		return comp;
	}

	template<typename _T2>
	inline scalar<_T>& operator=(const tensor<_T2>& T)
	{
		comp = (_T)T.comp;
		return *this;
	}

	template<typename _T2>
	inline scalar<_T>& operator=(const _T2& s)
	{
		comp = (_T)s;
		return *this;
	}

	inline _T operator-() const
	{
		return -comp;
	}

	//张量形式

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline scalar<_Result> operator+(const tensor<_T2>& rhs) const
	{
		return (_Result)(comp + rhs.comp);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline _Result operator-(const tensor<_T2>& rhs) const
	{
		return (_Result)(comp - rhs.comp);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const tensor<_T2>& rhs) const
	{
		return _Result(comp * rhs.comp);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline _Result operator/(const tensor<_T2>& rhs) const
	{
		return (_Result)(comp / rhs.comp);
	}

	template<typename _T2>
	inline tensor<_T>& operator+=(const tensor<_T2>& rhs)
	{
		comp += rhs.comp;
		return *this;
	}

	template<typename _T2>
	inline tensor<_T>& operator-=(const tensor<_T2>& rhs)
	{
		comp -= rhs.comp;
		return *this;
	}

	template<typename _T2>
	inline tensor<_T>& operator*=(const tensor<_T2>& rhs)
	{
		comp *= rhs.comp;
		return *this;
	}

	template<typename _T2>
	inline tensor<_T>& operator/=(const tensor<_T2>& rhs)
	{
		comp /= rhs.comp;
		return *this;
	}

	//标量形式

	template<typename _T2>
	inline tensor<_T>& operator+=(const _T2& rhs)
	{
		comp += rhs;
		return *this;
	}

	template<typename _T2>
	inline tensor<_T>& operator-=(const _T2& rhs)
	{
		comp -= rhs;
		return *this;
	}

	template<typename _T2>
	inline tensor<_T>& operator*=(const _T2& rhs)
	{
		comp *= rhs;
		return *this;
	}

	template<typename _T2>
	inline tensor<_T>& operator/=(const _T2& rhs)
	{
		comp /= rhs;
		return *this;
	}

	template<typename _T2>
	inline static type all(_T2 value)
	{
		return
		{	(_T)value};
	}

	inline static type zero()
	{
		return
		{	_T(0)};
	}

	inline static type one()
	{
		return
		{	_T(1)};
	}

	template<typename _Result = _T>
	inline scalar<_Result> unit() const
	{
		return
		{	_Result(1)};
	}
};

template<typename T>
constexpr size_t tensor<T>::dim[tensor<T>::order] = {};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() + tplmp::decl<_T2>::val())>
inline _Result operator+(const _T1& lhs, const tensor<_T2>& rhs)
{
	return (_Result)(lhs + rhs.comp);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() + tplmp::decl<_T2>::val())>
inline _Result operator+(const tensor<_T1>& lhs, const _T2& rhs)
{
	return (_Result)(lhs.comp + rhs);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() - tplmp::decl<_T2>::val())>
inline _Result operator-(const _T1& lhs, const tensor<_T2>& rhs)
{
	return (_Result)(lhs - rhs.comp);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() - tplmp::decl<_T2>::val())>
inline _Result operator-(const tensor<_T1>& lhs, const _T2& rhs)
{
	return (_Result)(lhs.comp - rhs);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline _Result operator*(const _T1& lhs, const tensor<_T2>& rhs)
{
	return (_Result)(lhs * rhs.comp);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline _Result operator*(const tensor<_T1>& lhs, const _T2& rhs)
{
	return (_Result)(lhs.comp * rhs);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline _Result operator/(const _T1& lhs, const tensor<_T2>& rhs)
{
	return (_Result)(lhs / rhs.comp);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline _Result operator/(const tensor<_T1>& lhs, const _T2& rhs)
{
	return (_Result)(lhs.comp / rhs);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline scalar<_Result> hadamard(const tensor<_T1>& lhs, const tensor<_T2>& rhs)
{
	return
	{	(_Result)(lhs.comp * rhs.comp)};
}

template<typename _T, typename _Result = _T>
inline _Result frobenius_norm(const tensor<_T>& T)
{
	return math::template sqrt<_Result>(T.comp * T.comp);
}

template<typename _T, size_t ..._Dims>
struct tensor
{
	typedef tensor<_T, _Dims...> type;
	typedef _T data_type;

	/**
	 * @brief 张量阶数
	 */
	static constexpr size_t order = sizeof...(_Dims);

	/**
	 * @brief 张量各个模的维度
	 */
	static const size_t dim[order];

	static constexpr size_t outmost_dim = dim[0];

	typedef tplmp::type_pack<tplmp::_size_t<_Dims>...> dim_t;

	template<size_t ... _OuterDims>
	struct higher_order
	{
		typedef typename higher_order_tensortpl<tensor, data_type, dim_t, _OuterDims...>::type type;
	};

	/**
	 * @brief 去除前_RmOrders个维度后剩余的维度的张量类型
	 */
	template<size_t _RmOrders, typename _DataType = _T>
	struct lower_order
	{
		typedef typename lower_order_tensortpl<tensor, _DataType, dim_t, _RmOrders>::type type;
	};

	typedef typename lower_order<1>::type comp_type;

	template<size_t _Length>
	struct slice_type
	{
		//生成{_Length, _Dims...}的索引，即0, 2, 3, ..., sizeof...(_Dims)
		typedef typename tensortpl_of_dims<tensor, _T,
				typename tplmp::index_sequence_t<size_t, 2, sizeof...(_Dims) - 1>::prepend‌<tplmp::_size_t<0> >::type,
				_Length, _Dims...>::type type;
	};

	comp_type comp[outmost_dim];

	/**
	 * @brief 包装为1维高阶张量
	 */
	inline operator typename higher_order<1>::type&()
	{
		return *(typename higher_order<1>::type*)this;
	}

	inline operator const typename higher_order<1>::type&() const
	{
		return *(const typename higher_order<1>::type*)this;
	}

	inline comp_type& operator[](size_t i)
	{
		return comp[i];
	}

	inline const comp_type& operator[](size_t i) const
	{
		return comp[i];
	}

	/**
	 * 切片，区间为[begin, begin + _Length)
	 */
	template<size_t _Length>
	inline typename slice_type<_Length>::type& slice(size_t begin)
	{
		return *(typename slice_type<_Length>::type*)(comp + begin);
	}

	template<size_t _Length>
	inline const typename slice_type<_Length>::type& slice(size_t begin) const
	{
		return *(typename slice_type<_Length>::type*)(comp + begin);
	}

	/**
	 * 类型强制转换
	 */
	template<typename _T2>
	inline operator tensor<_T2, _Dims...>() const
	{
		typedef tensor<_T2, _Dims...> result_type;
		result_type T;
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			T.comp[i] = (typename result_type::comp_type)comp[i];
		return T;
	}

	/**
	 * 强制转换赋值
	 */
	template<typename _T2>
	inline type& operator=(const tensor<_T2, _Dims...>& T)
	{
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			comp[i] = (comp_type)(T.comp[i]);
		return *this;
	}

	/**
	 * 按位加法
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline tensor<_Result, _Dims...> operator+(const tensor<_T2, _Dims...>& rhs) const
	{
		typedef tensor<_Result, _Dims...> result_type;
		result_type T;
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			T.comp[i] = (typename result_type::comp_type)(comp[i] + rhs.comp[i]);
		return T;
	}

	/**
	 * 按位减法
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline tensor<_Result, _Dims...> operator-(const tensor<_T2, _Dims...>& rhs) const
	{
		typedef tensor<_Result, _Dims...> result_type;
		result_type T;
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			T.comp[i] = (typename result_type::comp_type)(comp[i] - rhs.comp[i]);
		return T;
	}

	inline type operator-() const
	{
		type T;
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			T.comp[i] = -comp[i];
		return T;
	}

	template<typename _T2>
	inline type& operator+=(const tensor<_T2, _Dims...>& rhs)
	{
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			comp[i] += rhs.comp[i];
		return *this;
	}

	template<typename _T2>
	inline type& operator-=(const tensor<_T2, _Dims...>& rhs)
	{
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			comp[i] -= rhs.comp[i];
		return *this;
	}

	template<typename _T2>
	inline type& operator*=(const tensor<_T2>& rhs)
	{
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			comp[i] *= rhs.comp;
		return *this;
	}

	template<typename _T2>
	inline type& operator/=(const tensor<_T2>& rhs)
	{
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			comp[i] /= rhs.comp;
		return *this;
	}

	template<typename _T2>
	inline type& operator*=(const _T2& rhs)
	{
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			comp[i] *= rhs;
		return *this;
	}

	template<typename _T2>
	inline type& operator/=(const _T2& rhs)
	{
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			comp[i] /= rhs;
		return *this;
	}

	/**
	 * @brief Frobenius内积
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const tensor<_T2, _Dims...>& rhs) const
	{
		_Result s(0);
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			s += (comp[i] * rhs.comp[i]);
		return s;
	}

	/**
	 * @brief 按位除法，同哈达玛积
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline tensor<_Result, _Dims...> operator/(const tensor<_T2, _Dims...>& rhs) const
	{
		typedef tensor<_Result, _Dims...> result_type;
		result_type T;
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			T[i] = (typename result_type::comp_type)(comp[i] / rhs.comp[i]);
		return T;
	}

	/**
	 * 全部分量都相同的张量
	 */
	inline static type all(_T value)
	{
		type T;
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			T.comp[i] = comp_type::all(value);
		return T;
	}

	inline static type zero()
	{
		return all(0);
	}

	inline static type one()
	{
		return all(1);
	}

	/**
	 * 该方向单位张量
	 */
	template<typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T>::val())>
	inline tensor<_Result, _Dims...> unit() const;
};

template<typename T, size_t ... Dims>
constexpr size_t tensor<T, Dims...>::dim[tensor<T, Dims...>::order] = {Dims...};

/**
 * 标量乘法
 */
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> operator*(const _T1& lhs, const tensor<_T2, _Dims...>& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs * rhs.comp[i]);
	return T;
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> operator*(const tensor<_T1, _Dims...>& lhs, const _T2& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs.comp[i] * rhs);
	return T;
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> operator*(const tensor<_T1>& lhs, const tensor<_T2, _Dims...>& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs.comp * rhs.comp[i]);
	return T;
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> operator*(const tensor<_T1, _Dims...>& lhs, const tensor<_T2>& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs.comp[i] * rhs.comp);
	return T;
}

/**
 * @brief 矩阵乘法
 */
template<typename _T, size_t _Row, size_t _Column1, size_t _Column2>
inline tensor<_T, _Row, _Column2> operator*(const tensor<_T, _Row, _Column1>& lhs, const tensor<_T, _Column1, _Column2>& rhs)
{
	tensor<_T, _Row, _Column2> A;
	for(size_t i = 0; i < _Row; ++i)
	{
		for(size_t j = 0; j < _Column2; ++j)
		{
			_T sum(0);
			for(size_t k = 0; k < _Column1; ++k)
			{
				sum += lhs.comp[i][k] * rhs.comp[k][j];
			}
			A.comp[i][j] = sum;
		}
	}
	return A;
}

/**
 * @brief 矩阵*向量
 */
template<typename _T, size_t _Row, size_t _Column>
inline tensor<_T, _Row> operator*(const tensor<_T, _Row, _Column>& A, const tensor<_T, _Column>& x)
{
	tensor<_T, _Row> v;
	for(size_t i = 0; i < _Row; ++i)
	{
		v.comp[i] = A.comp[i] * x;
	}
	return v;
}

/**
 * 标量除法
 */
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> operator/(const tensor<_T1, _Dims...>& lhs, const _T2& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs.comp[i] / rhs);
	return T;
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> operator/(const _T1& lhs, const tensor<_T2, _Dims...>& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs / rhs.comp[i]);
	return T;
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> operator/(const tensor<_T1, _Dims...>& lhs, const tensor<_T2>& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs.comp[i] / rhs.comp);
	return T;
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> operator/(const tensor<_T1>& lhs, const tensor<_T2, _Dims...>& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs.comp / rhs.comp[i]);
	return T;
}

/**
 * 按位相乘
 */
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val()), size_t ..._Dims>
inline tensor<_Result, _Dims...> hadamard(const tensor<_T1, _Dims...>& lhs, const tensor<_T2, _Dims...>& rhs)
{
	typedef tensor<_Result, _Dims...> result_type;
	result_type T;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < result_type::outmost_dim; ++i)
		T.comp[i] = (typename result_type::comp_type)(lhs.comp[i] * rhs.comp[i]);
	return T;
}

//加法恒等元
template<typename _T, size_t ..._Dims>
struct __identity_impl<tensor<_T, _Dims...>, add>
{
	typedef tensor<_T, _Dims...> type;

	inline static type value()
	{
		return type::zero();
	}
};

//基
template<typename _T, size_t ..._Dims>
struct __basis_impl<tensor<_T, _Dims...> >
{
	typedef tensor<_T, _Dims...> type;

	inline static type value(size_t i)
	{
		type basis = type::zero();
		//basis[i] = 1;
		return basis;
	}
};

/**
 * 2维向量叉乘
 */
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline _Result operator^(const tensor<_T1, 2>& lhs, const tensor<_T2, 2>& rhs)
{
	return lhs.comp[0] * rhs.comp[1] - lhs.comp[1] * rhs.comp[0];
}

/**
 * 3维向量叉乘
 */
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<_Result, 3> operator^(const tensor<_T1, 3>& lhs, const tensor<_T2, 3>& rhs)
{
	return
	{	lhs.comp[1] * rhs.comp[2] - lhs.comp[2] * rhs.comp[1],
		lhs.comp[2] * rhs.comp[0] - lhs.comp[0] * rhs.comp[2],
		lhs.comp[0] * rhs.comp[1] - lhs.comp[1] * rhs.comp[0]};
}

template<typename _T, size_t _Dim1, size_t _Dim2>
inline vector<_T, _Dim1 + _Dim2> cat(const tensor<_T, _Dim1>& lhs, const tensor<_T, _Dim2>& rhs)
{
	vector<_T, _Dim1 + _Dim2> T;
	*(vector<_T, _Dim1>*)T.comp = lhs;
	*(vector<_T, _Dim2>*)(T.comp + _Dim1) = rhs;
	return T;
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const tensor<_T, 1, 1>& A)
{
	return A.comp[0][0];
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const tensor<_T, 2, 2>& A)
{
	return A.comp[0][0] * A.comp[1][1] - A.comp[0][1] * A.comp[1][0];
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const tensor<_T, 3, 3>& A)
{
	return A.comp[0][0] * (A.comp[1][1] * A.comp[2][2] - A.comp[1][2] * A.comp[2][1])
			+ A.comp[0][1] * (A.comp[1][2] * A.comp[2][0] - A.comp[1][0] * A.comp[2][2])
			+ A.comp[0][2] * (A.comp[1][0] * A.comp[2][1] - A.comp[1][1] * A.comp[2][0]);
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const tensor<_T, 4, 4>& A)
{
	return A.comp[0][0] * (A.comp[1][1] * A.comp[2][2] * A.comp[3][3] - A.comp[3][1] * A.comp[2][2] * A.comp[1][3])
			+ A.comp[0][1] * (A.comp[1][2] * A.comp[2][3] * A.comp[3][0] - A.comp[1][0] * A.comp[2][3] * A.comp[3][2])
			+ A.comp[0][2] * (A.comp[1][3] * A.comp[2][0] * A.comp[3][1] - A.comp[1][1] * A.comp[2][0] * A.comp[3][3])
			+ A.comp[0][3] * (A.comp[1][0] * A.comp[2][1] * A.comp[3][2] - A.comp[1][2] * A.comp[2][1] * A.comp[3][0]);
}

// ------------------ 矩阵 ------------------
/**
 * 找列主元。
 * 如果找到则返回行索引，未找到（全0）则返回_Row
 */
template<typename _T, size_t _Row, size_t _Column>
size_t partial_pivot_row(const tensor<_T, _Row, _Column>& A, size_t i, size_t j, _T eps)
{
	size_t pivot_row = _Row;
	_T max_abs = _T(0);
	for(size_t r = i; r < _Row; ++r)
	{
		_T abs_val = math::abs(A.comp[r][j]);
		if(abs_val > max_abs)
		{
			max_abs = abs_val;
			pivot_row = r;
		}
	}
	if(math::is_zero(max_abs, eps))
		return _Row;
	else
		return pivot_row;
}

/**
 * 找全主元。
 * 如果找到则返回行索引和列索引，未找到则返回row。
 */
template<typename _T, size_t _Row, size_t _Column>
vector<size_t, 2> complete_pivot_comp(const tensor<_T, _Row, _Column>& A, size_t i, size_t j, _T eps)
{
	vector<size_t, 2> pivot_comp = {_Row, _Column};
	_T max_abs = _T(0);
	for(size_t r = i; r < _Row; ++r)
	{
		for(size_t c = j; c < _Column; ++c)
		{
			_T abs_val = math::abs(A.comp[r][c]);
			if(abs_val > max_abs)
			{
				max_abs = abs_val;
				pivot_comp[0] = r;
				pivot_comp[1] = c;
			}
		}
	}
	if(math::is_zero(max_abs, eps))
		return
		{	_Row, _Column};
	else
		return pivot_comp;
}

/**
 * 交换本矩阵的行
 */
template<typename _T, size_t _Row, size_t _Column>
inline matrix<_T, _Row, _Column> swap_row(const tensor<_T, _Row, _Column>& A, size_t row1, size_t row2)
{
	vector<_T, _Column> temp = A.comp[row1];
	A.comp[row1] = A.comp[row2];
	A.comp[row2] = temp;
	return A;
}

/**
 * 交换本矩阵的列
 */
template<typename _T, size_t _Row, size_t _Column>
inline matrix<_T, _Row, _Column> swap_column(const tensor<_T, _Row, _Column>& A, size_t column1, size_t column2)
{
	__tensor_loop_unroll__()
	for(size_t i = 0; i < _Row; ++i)
	{
		_T temp = A.comp[i][column1];
		A.comp[i][column1] = A.comp[i][column2];
		A.comp[i][column2] = temp;
	}
	return A;
}

/**
 * 矩阵分块，分块区间为[_BeginRow, _EndRow)，[_BeginColumn, _EndColumn)
 */
template<size_t _PartitionRow, size_t _PartitionColumn, typename _T, size_t _Row, size_t _Column>
inline matrix<_T, _PartitionRow, _PartitionColumn> partition(const tensor<_T, _Row, _Column>& A, size_t row_begin, size_t column_begin)
{
	static_assert(_PartitionRow >= 0 && _PartitionColumn >= 0, "invalid partition size");
	matrix<_T, _PartitionRow, _PartitionColumn> partition_result;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < _PartitionRow; ++i)
	{
		partition_result.comp[i] = A.comp[row_begin + i].template slice<_PartitionColumn>(column_begin);
	}
	return partition_result;
}

template<typename _T, size_t _Row, size_t _Column1, size_t _Column2>
inline matrix<_T, _Row, _Column1 + _Column2> cat_column(const tensor<_T, _Row, _Column1>& A1, const tensor<_T, _Row, _Column2>& A2)
{
	matrix<_T, _Row, _Column1 + _Column2> A;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < _Row; ++i)
		A.comp[i] = cat(A1.comp[i], A2.comp[i]);
	return A;
}

/**
 * @brief 4个分块矩阵拼接
 */
template<typename _T, size_t _Row1, size_t _Column1, size_t _Row2, size_t _Column2>
inline matrix<_T, _Row1 + _Row2, _Column1 + _Column2> cat(const tensor<_T, _Row1, _Column1>& left_up, const tensor<_T, _Row1, _Column2>& right_up, const tensor<_T, _Row2, _Column1>& left_down, const tensor<_T, _Row2, _Column2>& right_down)
{
	constexpr size_t cat_row = _Row1 + _Row2;
	constexpr size_t cat_column = _Column1 + _Column2;
	matrix<_T, cat_row, cat_column> A;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < _Row1; ++i)
		A.comp[i] = cat(left_up.comp[i], right_up.comp[i]);
	__tensor_loop_unroll__()
	for(size_t i = _Row1; i < cat_row; ++i)
		A.comp[i] = cat(left_down.comp[i], right_down.comp[i]);
	return A;
}

/**
 * Row Echelon Form，即变换为行阶梯形
 * 列主元高斯消元法计算
 */
template<typename _T, size_t _Row, size_t _Column, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
matrix<_Result, _Row, _Column> ref(const tensor<_T, _Row, _Column>& A, _T eps, size_t* mat_rank = nullptr, size_t* swap_count = nullptr)
{
	if(swap_count)
		*swap_count = 0;
	size_t rank = 0; //当前有主元的行数，即秩
	matrix<_Result, _Row, _Column> ref = A;
	for(size_t current_column = 0; current_column < _Column && rank < _Row; ++current_column)
	{
		//找列主元
		size_t pivot_row = ref.partial_pivot_row(rank, current_column, eps);
		if(pivot_row < _Row)
		{
			if(pivot_row != rank)
			{
				ref.swap_row(rank, pivot_row); //当前主元所在行不同则交换当前行与主元行
				if(swap_count)
					++*swap_count;
			}
			_T pivot_val = ref[rank][current_column];
			for(size_t r = rank + 1; r < _Row; ++r)
			{
				//消去左下方元素得到上三角矩阵
				_T factor = ref[r][current_column] / pivot_val;
				for(size_t c = current_column; c < _Column; ++c)
				{
					ref[r][c] -= (_Result)(factor * ref[current_column][c]);
					ref[r][current_column] = _Result(0); //直接置零避免浮点数误差
				}
			}
			++rank; //找到主元才增加行号
		}
	}
	if(mat_rank)
		*mat_rank = rank;
	return ref;
}

/**
 * @brief REF回代法求解线性方程组
 * 		  即在已经矩阵A、向量b的情况下，求满足Ax=b的x向量
 */
template<typename _T, size_t _Order, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
vector<_Result, _Order> solve_linear_system(const tensor<_T, _Order, _Order>& coeff_mat, const tensor<_T, _Order>& vals, bool* success, _T eps = 0)
{
	//将系数矩阵和值向量横向拼接成增广矩阵
	matrix<_T, _Order, _Order + 1> augmented;
	for(size_t i = 0; i < _Order; ++i)
	{
		for(size_t j = 0; j < _Order; ++j)
			augmented[i][j] = coeff_mat[i][j];
		augmented[i][_Order] = vals[i];
	}
	vector<_Result, _Order> solution;
	size_t rank;
	auto ref_mat = ref(augmented, eps, &rank);
	//有解时增广矩阵必定满秩，且REF是上梯形矩阵且对角线不为0，直接从最后一行回代
	if(rank < _Order)
	{
		//无穷多解
		goto FAILED;
	}
	else
	{
		for(size_t i = _Order; i--;)
		{
			_Result sum = ref_mat[i][_Order]; //ref_mat是增广矩阵，维度为_Order x _Order+1
			for(size_t j = i + 1; j < _Order; ++j)
			{
				sum -= ref_mat[i][j] * solution[j];
			}
			if(is_zero(ref_mat[i][i], eps))
			{
				//rank<_Order时可能无解或无穷多解
				goto FAILED;
			}
			solution[i] = sum / ref_mat[i][i];
		}
	}
	if(success)
		*success = true;
	return solution;
	FAILED:
	if(success)
		*success = false;
	return solution;
}

/**
 * 秩
 * 相对于ref()的优化
 */
template<typename _T, size_t _Row, size_t _Column, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
size_t rank(const tensor<_T, _Row, _Column>& A, _T eps = 0)
{
	matrix<_Result, _Row, _Column> ref = A;
	size_t rank = 0;
	for(size_t current_column = 0; current_column < _Column; ++current_column)
	{
		if(rank >= _Row)
			break; //满秩直接返回
		size_t pivot_row = ref.partial_pivot_row(rank, current_column, eps);
		if(pivot_row >= _Row)
			continue;  //该列全零，跳过
		if(pivot_row != rank)
		{
			ref.swap_row(rank, pivot_row);
		}
		_T pivot_val = ref[rank][current_column];
		for(size_t r = rank + 1; r < _Row; ++r)
		{
			_T factor = ref[r][current_column] / pivot_val;
			for(size_t c = current_column; c < _Column; ++c)
			{
				ref[r][c] -= factor * ref[rank][c];
			}
		}
		++rank;
	}
	return rank;
}

/**
 * 行列式
 * 相对于ref()的优化
 */
template<typename _T, size_t _Order, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
_Result det(const tensor<_T, _Order, _Order>& A, _T eps = 0)
{
	matrix<_T, _Order, _Order> temp = A;
	_Result det = _Result(1);
	for(size_t i = 0; i < _Order; ++i)
	{
		//找列主元
		size_t pivot_row = temp.partial_pivot_row(i, i, eps);
		if(pivot_row < _Order)
		{
			if(pivot_row != i)
			{
				temp.swap_row(i, pivot_row); //当前主元所在行不同则交换当前行与主元行
				det = -det; //交换一次行列则行列式取反
			}
			_T pivot_val = temp[i][i];
			for(size_t r = i + 1; r < _Order; ++r)
			{
				//消去左下方元素得到上三角矩阵
				_T factor = temp[r][i] / pivot_val;
				for(size_t c = i; c < _Order; ++c)
				{
					temp[r][c] -= factor * temp[i][c]; //由于只有对角线元素参与最终计算，下三角的微小误差不需要再置零
				}
			}
			//累乘对角线
			det *= temp[i][i];
		}
		else
		{
			//任意一列没找到主元则行列式恒等于0
			return _Result(0);
		}
	}
	return det;
}

//方阵特有

/**
 * 迹
 */
template<typename _T, size_t _Order, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
inline _Result tr(const tensor<_T, _Order, _Order>& A)
{
	_Result tr = 0;
	__tensor_loop_unroll__()
	for(size_t i = 0; i < _Order; ++i)
		tr += A.comp[i][i];
	return tr;
}

/**
 * 恒等元
 */

//乘法恒等元：恒等矩阵
template<typename _T, size_t _Order>
struct __identity_impl<matrix<_T, _Order, _Order>, mul>
{
	inline static matrix<_T, _Order, _Order> value()
	{
		matrix<_T, _Order, _Order> I = matrix<_T, _Order, _Order>::zero();
		__tensor_loop_unroll__()
		for(size_t i = 0; i < _Order; ++i)
			I[i][i] = 1;
		return I;
	}
};

/*
 template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val()), size_t ..._Dims>
 inline _Result contract(const tensor<_T1, _Dims...>& lhs, const tensor<_T2, _Dims...>& rhs)
 {
 _T result(0);
 __tensor_loop_unroll__()
 for(size_t i = 0; i < outmost_dim; ++i)
 result += comp[i];
 return result;
 }
 */

/**
 * Frobenius范数，对向量而言数值为向量长度
 */
template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val()), size_t ...Dims>
inline _Result frobenius_norm(const tensor<_T, Dims...>& T)
{
	return math::template sqrt<_Result>(T * T);
}

/**
 * 单位张量
 */
template<typename _T, size_t ..._Dims>
template<typename _Result>
inline tensor<_Result, _Dims...> tensor<_T, _Dims...>::unit() const
{
	return operator/<_T, _T, _Result, _Dims...>(*this, frobenius_norm<_T, _Result>(*this));
}

/**
 * 字符串转换及打印
 */

template<typename _T>
inline std::string to_string(const tensor<_T>& T)
{
	std::ostringstream oss;
	oss << T.comp;
	return oss.str();
}

template<typename _T, size_t _Dim>
inline std::string to_string(const tensor<_T, _Dim>& T)
{
	std::ostringstream oss;
	oss << '(';
	for(size_t i = 0; i < _Dim; ++i)
	{
		oss << (_T)T.comp[i];
		if(i < _Dim - 1)
			oss << ", ";
	}
	oss << ')';
	return oss.str();
}

template<typename _T, size_t _Row, size_t _Column>
inline std::string to_string(const tensor<_T, _Row, _Column>& T)
{
	std::ostringstream oss;
	for(size_t i = 0; i < _Row; ++i)
	{
		oss << '[';
		for(size_t j = 0; j < _Column; ++j)
		{
			oss << (_T)T.comp[i][j];
			if(j < _Column - 1)
				oss << '\t';
		}
		oss << ']';
		if(i < _Row - 1)
			oss << '\n';
	}
	return oss.str();
}

template<typename _T, size_t ..._Dims>
inline std::ostream& operator<<(std::ostream& os, const tensor<_T, _Dims...>& T)
{
	os << to_string(T);
	return os;
}

/**
 * @brief 张量子视图，可以直接在其上进行读写
 */
template<typename _T, size_t ..._Dims>
class tensor_view
{
protected:
	_T* comp = nullptr;

public:
	typedef typename tensor<_T, _Dims...>::type type;
	typedef tensor_view<_T, _Dims...> view_type;

	typedef typename tensor<_T, _Dims...>::data_type data_type;
	typedef typename tensor<_T, _Dims...>::comp_type comp_type;

	static const size_t order = type::oder;
	static constexpr size_t dim[order] = type::dim;

	static constexpr size_t outmost_dim = type::outmost_dim;

	tensor_view() = default;

	template<size_t ..._Dims2>
	tensor_view(tensor<_T, _Dims2...>& target, size_t row_begin, size_t column_begin)
	{
		static_assert(tensor<_T, _Dims2...>::outmost_dim >= outmost_dim && sizeof...(_Dims) == sizeof...(_Dims2), "tensor view dim should be less than target");
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
		{
			//必须确保column_begin+_Column<=_Column2
			comp[i] = target.comp[row_begin + i].comp + column_begin;
		}
	}

	template<typename _T2>
	inline view_type& operator=(const tensor<_T2, _Dims...>& T)
	{
		__tensor_loop_unroll__()
		for(size_t i = 0; i < outmost_dim; ++i)
			*(comp_type*)comp[i] = T.comp[i];
		return *this;
	}

	inline operator type() const
	{
		type T;
		for(size_t i = 0; i < outmost_dim; ++i)
			T.comp[i] = *(comp_type*)comp[i];
		return T;
	}
};

template<size_t ..._Dims>
using tensorf = tensor<float, _Dims...>;
template<size_t ..._Dims>
using tensord = tensor<double, _Dims...>;
template<size_t ..._Dims>
using tensorq = tensor<__float128, _Dims...>;

template<typename _T>
using vector1 = vector<_T, 1>;
template<typename _T>
using vector2 = vector<_T, 2>;
template<typename _T>
using vector3 = vector<_T, 3>;
template<typename _T>
using vector4 = vector<_T, 4>;
template<typename _T>
using vector5 = vector<_T, 5>;
template<typename _T>
using vector6 = vector<_T, 6>;
template<typename _T>
using vector7 = vector<_T, 7>;
template<typename _T>
using vector8 = vector<_T, 8>;
template<typename _T>
using vector9 = vector<_T, 9>;

using vector1f = vector1<float>;
using vector2f = vector2<float>;
using vector3f = vector3<float>;
using vector4f = vector4<float>;
using vector5f = vector5<float>;
using vector6f = vector6<float>;
using vector7f = vector7<float>;
using vector8f = vector8<float>;
using vector9f = vector9<float>;

using vector1d = vector1<double>;
using vector2d = vector2<double>;
using vector3d = vector3<double>;
using vector4d = vector4<double>;
using vector5d = vector5<double>;
using vector6d = vector6<double>;
using vector7d = vector7<double>;
using vector8d = vector8<double>;
using vector9d = vector9<double>;

using vector1q = vector1<__float128>;
using vector2q = vector2<__float128>;
using vector3q = vector3<__float128>;
using vector4q = vector4<__float128>;
using vector5q = vector5<__float128>;
using vector6q = vector6<__float128>;
using vector7q = vector7<__float128>;
using vector8q = vector8<__float128>;
using vector9q = vector9<__float128>;

template<typename _T>
using matrix1x1 = matrix<_T, 1, 1>;
template<typename _T>
using matrix2x2 = matrix<_T, 2, 2>;
template<typename _T>
using matrix3x3 = matrix<_T, 3, 3>;
template<typename _T>
using matrix4x4 = matrix<_T, 4, 4>;
template<typename _T>
using matrix5x5 = matrix<_T, 5, 5>;
template<typename _T>
using matrix6x6 = matrix<_T, 6, 6>;
template<typename _T>
using matrix7x7 = matrix<_T, 7, 7>;
template<typename _T>
using matrix8x8 = matrix<_T, 8, 8>;
template<typename _T>
using matrix9x9 = matrix<_T, 9, 9>;

using matrix1x1f = matrix1x1<float>;
using matrix2x2f = matrix2x2<float>;
using matrix3x3f = matrix3x3<float>;
using matrix4x4f = matrix4x4<float>;
using matrix5x5f = matrix5x5<float>;
using matrix6x6f = matrix6x6<float>;
using matrix7x7f = matrix7x7<float>;
using matrix8x8f = matrix8x8<float>;
using matrix9x9f = matrix9x9<float>;

using matrix1x1d = matrix1x1<double>;
using matrix2x2d = matrix2x2<double>;
using matrix3x3d = matrix3x3<double>;
using matrix4x4d = matrix4x4<double>;
using matrix5x5d = matrix5x5<double>;
using matrix6x6d = matrix6x6<double>;
using matrix7x7d = matrix7x7<double>;
using matrix8x8d = matrix8x8<double>;
using matrix9x9d = matrix9x9<double>;

using matrix1x1q = matrix1x1<__float128>;
using matrix2x2q = matrix2x2<__float128>;
using matrix3x3q = matrix3x3<__float128>;
using matrix4x4q = matrix4x4<__float128>;
using matrix5x5q = matrix5x5<__float128>;
using matrix6x6q = matrix6x6<__float128>;
using matrix7x7q = matrix7x7<__float128>;
using matrix8x8q = matrix8x8<__float128>;
using matrix9x9q = matrix9x9<__float128>;
}

#endif//_MATH_TENSOR
