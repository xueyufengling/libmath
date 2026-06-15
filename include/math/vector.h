#ifndef _MATH_VECTOR
#define _MATH_VECTOR

#include <tplmp/base.h>

#include <string.h>
#include <sstream>

#include <math/algebra.h>

#define __vector_loop_unroll__() __loop_unroll__(9)

namespace math
{
template<size_t _Row, size_t _Column, typename _T>
class matrix;

/**
 * 列向量，初始化使用vector vec={x,y,z,w...};
 */
template<size_t _Dim, typename _T>
class vector
{
public:
	_T coord[_Dim];
	typedef _T type;
	static const size_t dim = _Dim;

	inline operator _T*()
	{
		return coord;
	}
	template<size_t _Dim2>
	inline operator math::vector<_Dim2, _T>&()
	{
		return *(math::vector<_Dim2, _T>*)this;
	}

	/**
	 * 坐标
	 */
	inline _T& operator[](size_t i)
	{
		return coord[i];
	}

	inline const _T& operator[](size_t i) const
	{
		return coord[i];
	}

	/**
	 * 向量长度
	 */
	template<typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
	inline _Result norm() const
	{
		return math::sqrt<_Result>(((*this)) * (*this));
	}

	/**
	 * 单位向量
	 */
	template<typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T>::val())>
	inline vector<_Dim, _Result> unit() const
	{
		return (*this) / norm();
	}

	/**
	 * 视为1 x _Dim或_Dim x 1矩阵
	 */
	template<int _D = _Dim>
	operator typename tplmp::if_else<(_D != 1)>::def<math::matrix<_D, 1, _T>&>() const
	{
		return *(math::matrix<_D, 1, _T>*)this;
	}

	template<int _D = _Dim>
	operator typename tplmp::if_else<(_D != 1)>::def<math::matrix<1, _D, _T>&>() const
	{
		return *(math::matrix<1, _D, _T>*)this;
	}

	template<int _D = _Dim>
	operator typename tplmp::if_else<(_D == 1)>::def<math::matrix<_D, 1, _T>&>() const
	{
		return *(math::matrix<1, 1, _T>*)this;
	}

	/**
	 * 向量切片，区间为[begin, begin + _Length)
	 */
	template<size_t _Length>
	inline vector<_Length, _T>& slice(size_t begin)
	{
		return *(vector<_Length, _T>*)(coord + begin);
	}

	template<size_t _Length>
	inline const vector<_Length, _T>& slice(size_t begin) const
	{
		return *(const vector<_Length, _T>*)(coord + begin);
	}

	/**
	 * 向量切片，区间为[_Begin, _End)，剩余值以指定值m填充
	 */
	template<typename _Result = _T>
	inline vector<_Dim, _Result> mask(size_t begin, size_t end, _Result mask_value) const
	{
		vector<_Dim, _Result> mask_vec;
		for(size_t i = 0; i < begin; ++i)
			mask_vec.coord[i] = mask_value;
		for(size_t i = begin; i < end; ++i)
			mask_vec.coord[i] = coord[i];
		for(size_t i = end; i < _Dim; ++i)
			mask_vec.coord[i] = mask_value;
		return mask_vec;
	}

	/**
	 * 类型强制转换
	 */
	template<typename _T2>
	inline operator vector<_Dim, _T2>() const
	{
		vector<_Dim, _T2> cast_result;
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			cast_result[i] = (_T2)coord[i];
		return cast_result;
	}

	/**
	 * 强制转换赋值
	 */
	template<typename _T2>
	inline vector<_Dim, _T>& operator=(const vector<_Dim, _T2>& vec)
	{
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			coord[i] = (_T)(vec.coord[i]);
		return *this;
	}

	/**
	 * 向量加法
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline vector<_Dim, _Result> operator+(const vector<_Dim, _T2>& vec) const
	{
		vector<_Dim, _Result> add_result;
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			add_result[i] = (_Result)(coord[i] + vec.coord[i]);
		return add_result;
	}

	/**
	 * 向量减法
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline vector<_Dim, _Result> operator-(const vector<_Dim, _T2>& vec) const
	{
		vector<_Dim, _Result> sub_result;
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			sub_result[i] = (_Result)(coord[i] - vec.coord[i]);
		return sub_result;
	}

	template<typename _T2>
	inline vector<_Dim, _T>& operator+=(const vector<_Dim, _T2>& vec)
	{
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			coord[i] += vec.coord[i];
		return *this;
	}

	template<typename _T2>
	inline vector<_Dim, _T>& operator-=(const vector<_Dim, _T2>& vec)
	{
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			coord[i] -= vec.coord[i];
		return *this;
	}

	/**
	 * 向量取反
	 */
	inline vector<_Dim, _T> operator-() const
	{
		vector<_Dim, _T> inv;
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			inv[i] = -coord[i];
		return inv;
	}

	template<typename _T2>
	inline vector<_Dim, _T>& operator*=(const _T2& t)
	{
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			coord[i] *= t;
		return *this;
	}

	template<typename _T2>
	inline vector<_Dim, _T>& operator/=(const _T2& t)
	{
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			coord[i] /= t;
		return *this;
	}

	/**
	 * 向量点乘
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<_Dim, _T2> vec) const
	{
		_Result dot_result = 0;
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			dot_result += coord[i] * vec.coord[i];
		return dot_result;
	}

	/**
	 * 按位除法，同哈达玛积
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<_Dim, _Result> operator/(const vector<_Dim, _T2> vec) const
	{
		vector<_Dim, _Result> div_result;
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			div_result[i] = (_Result)(coord[i] / vec.coord[i]);
		return div_result;
	}

	/**
	 * 全部分量都相同的向量
	 */
	inline static vector<_Dim, _T> all(_T value)
	{
		vector<_Dim, _T> vec;
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			vec[i] = value;
		return vec;
	}

	inline static vector<_Dim, _T> zero()
	{
		return all(0);
	}

	inline static vector<_Dim, _T> one()
	{
		return all(1);
	}

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum() const
	{
		_T result(0);
		__vector_loop_unroll__()
		for(size_t i = 0; i < _Dim; ++i)
			result += coord[i];
		return result;
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(';
		for(size_t i = 0; i < _Dim; ++i)
		{
			oss << coord[i];
			if(i < _Dim - 1)
				oss << ", ";
		}
		oss << ')';
		return oss.str();
	}
};

template<size_t _Dim1, size_t _Dim2, typename _T>
inline vector<_Dim1 + _Dim2, _T> cat(const vector<_Dim1, _T>& vec1, const vector<_Dim2, _T>& vec2)
{
	vector<_Dim1 + _Dim2, _T> vec;
	*(vector<_Dim1, _T>*)vec.coord = vec1;
	*(vector<_Dim2, _T>*)(vec.coord + _Dim1) = vec2;
	return vec;
}

template<size_t _Dim, typename _T>
inline std::ostream& operator<<(std::ostream& os, const vector<_Dim, _T>& vec)
{
	os << ((std::string)vec);
	return os;
}

/**
 * 标量乘法
 */
template<typename _T1, size_t _Dim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> operator*(const vector<_Dim, _T1>& vec, _T2 t)
{
	vector<_Dim, _Result> mul_result;
	__vector_loop_unroll__()
	for(size_t i = 0; i < _Dim; ++i)
		mul_result[i] = (_Result)(vec.coord[i] * t);
	return mul_result;
}

/**
 * 标量除法
 */
template<typename _T1, size_t _Dim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> operator/(const vector<_Dim, _T1>& vec, _T2 t)
{
	vector<_Dim, _Result> div_result;
	__vector_loop_unroll__()
	for(size_t i = 0; i < _Dim; ++i)
		div_result[i] = (_Result)(vec.coord[i] / t);
	return div_result;
}

//交换律
template<typename _T1, size_t _Dim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> operator*(_T1 t, const vector<_Dim, _T2>& vec)
{
	vector<_Dim, _Result> mul_result;
	__vector_loop_unroll__()
	for(size_t i = 0; i < _Dim; ++i)
		mul_result[i] = (_Result)(t * vec.coord[i]);
	return mul_result;
}

template<typename _T1, size_t _Dim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> operator/(_T1 t, const vector<_Dim, _T2>& vec)
{
	vector<_Dim, _Result> div_result;
	__vector_loop_unroll__()
	for(size_t i = 0; i < _Dim; ++i)
		div_result[i] = (_Result)(t / vec.coord[i]);
	return div_result;
}

template<size_t _Dim, typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> hadamard(const vector<_Dim, _T1>& vec1, const vector<_Dim, _T2>& vec2)
{
	vector<_Dim, _Result> hadamard_result;
	__vector_loop_unroll__()
	for(size_t i = 0; i < _Dim; ++i)
		hadamard_result[i] = vec1[i] * vec2[i];
	return hadamard_result;
}

//加法恒等元
template<size_t _Dim, typename _T>
struct __identity_impl<vector<_Dim, _T>, add>
{
	inline static vector<_Dim, _T> value()
	{
		return vector<_Dim, _T>::zero();
	}
};

//基向量
template<size_t _Dim, typename _T>
struct __basis_impl<vector<_Dim, _T> >
{
	static vector<_Dim, _T> value(size_t i)
	{
		vector<_Dim, _T> basis = vector<_Dim, _T>::zero();
		basis[i] = 1;
		return basis;
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline _Result operator^(const vector<3, _T1>& vec1, const vector<3, _T2>& vec2)
{
	return vec1.coord[0] * vec2.coord[1] - vec1.coord[1] * vec2.coord[0];
}

/**
 * 3维向量叉乘
 */
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<3, _Result> operator^(const vector<3, _T1>& vec1, const vector<3, _T2>& vec2)
{
	return
	{	vec1.coord[1] * vec2.coord[2] - vec1.coord[2] * vec2.coord[1],
		vec1.coord[2] * vec2.coord[0] - vec1.coord[0] * vec2.coord[2],
		vec1.coord[0] * vec2.coord[1] - vec1.coord[1] * vec2.coord[0]};
}

template<typename _T>
using vector1 = vector<1, _T>;
template<typename _T>
using vector2 = vector<2, _T>;
template<typename _T>
using vector3 = vector<3, _T>;
template<typename _T>
using vector4 = vector<4, _T>;
template<typename _T>
using vector5 = vector<5, _T>;
template<typename _T>
using vector6 = vector<6, _T>;
template<typename _T>
using vector7 = vector<7, _T>;
template<typename _T>
using vector8 = vector<8, _T>;
template<typename _T>
using vector9 = vector<9, _T>;

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
}

#endif//_MATH_VECTOR