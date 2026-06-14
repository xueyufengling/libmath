#ifndef _MATH_VECTOR
#define _MATH_VECTOR

#include <tplmp/base.h>

#include <string.h>
#include <sstream>

#include <math/algebra.h>

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

#define __def_vector_members__(coord_type, vec_dim)\
	coord_type coord[vec_dim];\
	typedef coord_type type;\
	static const size_t dim = vec_dim;
	__def_vector_members__(_T, _Dim)

#define __def_vector_op_cast__(coord_type)\
	inline operator coord_type*()\
	{\
		return coord;\
	}\
	template<size_t _Dim2>\
	inline operator math::vector<_Dim2, elem_type>&()\
	{\
		return *(math::vector<_Dim2, elem_type>*)this;\
	}
	__def_vector_op_cast__(_T)

	/**
	 * 坐标
	 */
#define __def_vector_op_coord__(coord_type)\
	inline coord_type& operator[](size_t i)\
	{\
		return coord[i];\
	}\
	inline const coord_type& operator[](size_t i) const\
	{\
		return coord[i];\
	}
	__def_vector_op_coord__(_T)

	/**
	 * 向量长度
	 */
#define __def_vector_op_norm__(coord_type)\
	template<typename _Result = decltype(tplmp::decl<coord_type>::val() * tplmp::decl<coord_type>::val())>\
	inline _Result norm() const\
	{\
		return math::sqrt<_Result>(((*this)) * (*this));\
	}
	__def_vector_op_norm__(_T)

	/**
	 * 单位向量
	 */
#define __def_vector_op_unit__(coord_type, vec_dim)\
	template<typename _Result = decltype(tplmp::decl<coord_type>::val() / tplmp::decl<coord_type>::val())>\
	inline vector<vec_dim, _Result> unit() const\
	{\
		return (*this) / norm();\
	}
__def_vector_op_unit__	(_T, _Dim)

	/**
	 * 视为1 x _Dim或_Dim x 1矩阵
	 */
#define __def_vector_op_to_matrix__(coord_type, vec_dim)\
	operator math::matrix<vec_dim, 1, coord_type>&() const\
	{\
		return *(math::matrix<vec_dim, 1, coord_type>*)this;\
	}\
	operator math::matrix<1, vec_dim, coord_type>&() const\
	{\
		return *(math::matrix<1, vec_dim, coord_type>*)this;\
	}
	__def_vector_op_to_matrix__(_T, _Dim)

	/**
	 * 向量切片，区间为[begin, begin + _Length)
	 */
#define __def_vector_op_slice__(coord_type, vec_dim)\
	template<size_t _Length>\
	inline vector<_Length, coord_type>& slice(size_t begin)\
	{\
		return *(vector<_Length, coord_type>*)(coord + begin);\
	}\
	template<size_t _Length>\
	inline const vector<_Length, coord_type>& slice(size_t begin) const\
	{\
		return *(const vector<_Length, coord_type>*)(coord + begin);\
	}
	__def_vector_op_slice__(_T, _Dim)

	/**
	 * 向量切片，区间为[_Begin, _End)，剩余值以指定值m填充
	 */
#define __def_vector_op_mask__(coord_type, vec_dim)\
	template<typename _Result = coord_type>\
	inline vector<vec_dim, _Result> mask(size_t begin, size_t end, _Result mask_value) const\
	{\
		vector<vec_dim, _Result> mask_vec;\
		for(size_t i = 0; i < begin; ++i)\
			mask_vec.coord[i] = mask_value;\
		for(size_t i = begin; i < end; ++i)\
			mask_vec.coord[i] = coord[i];\
		for(size_t i = end; i < vec_dim; ++i)\
			mask_vec.coord[i] = mask_value;\
		return mask_vec;\
	}
	__def_vector_op_mask__(_T, _Dim)

	/**
	 * 类型强制转换
	 */
	template<typename _T2>
	inline operator vector<_Dim, _T2>() const
	{
		vector<_Dim, _T2> cast_result;
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
		for(size_t i = 0; i < _Dim; ++i)
		sub_result[i] = (_Result)(coord[i] - vec.coord[i]);
		return sub_result;
	}

	template<typename _T2>
	inline vector<_Dim, _T>& operator+=(const vector<_Dim, _T2>& vec)
	{
		for(size_t i = 0; i < _Dim; ++i)
		coord[i] += vec.coord[i];
		return *this;
	}

	template<typename _T2>
	inline vector<_Dim, _T>& operator-=(const vector<_Dim, _T2>& vec)
	{
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
		for(size_t i = 0; i < _Dim; ++i)
		inv[i] = -coord[i];
		return inv;
	}

	template<typename _T2>
	inline vector<_Dim, _T>& operator*=(const _T2& t)
	{
		for(size_t i = 0; i < _Dim; ++i)
		coord[i] *= t;
		return *this;
	}

	template<typename _T2>
	inline vector<_Dim, _T>& operator/=(const _T2& t)
	{
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
		for(size_t i = 0; i < _Dim; ++i)
		vec[i] = value;
		return vec;
	}

#define __def_vector_zero_one__(coord_type, vec_dim)\
	inline static vector<vec_dim, coord_type> zero()\
	{\
		return all(0);\
	}\
	inline static vector<vec_dim, coord_type> one()\
	{\
		return all(1);\
	}
	__def_vector_zero_one__(_T, _Dim)

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum() const
	{
		_T result(0);
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
	for(size_t i = 0; i < _Dim; ++i)
		div_result[i] = (_Result)(vec.coord[i] / t);
	return div_result;
}

//交换律
template<typename _T1, size_t _Dim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> operator*(_T1 t, const vector<_Dim, _T2>& vec)
{
	vector<_Dim, _Result> mul_result;
	for(size_t i = 0; i < _Dim; ++i)
		mul_result[i] = (_Result)(t * vec.coord[i]);
	return mul_result;
}

template<typename _T1, size_t _Dim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> operator/(_T1 t, const vector<_Dim, _T2>& vec)
{
	vector<_Dim, _Result> div_result;
	for(size_t i = 0; i < _Dim; ++i)
		div_result[i] = (_Result)(t / vec.coord[i]);
	return div_result;
}

template<size_t _Dim, typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> hadamard(const vector<_Dim, _T1>& vec1, const vector<_Dim, _T2>& vec2)
{
	vector<_Dim, _Result> hadamard_result;
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

// 1维向量特化
template<typename _T>
class vector<1, _T>
{
public:
	__def_vector_members__(_T, 1)

	__def_vector_op_cast__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_slice__(_T, 1)

	__def_vector_op_mask__(_T, 1)

template	<typename _Result = _T>
	inline _Result norm() const
	{
		return (_Result)coord[0];
	}

	template<typename _Result = _T>
	inline vector<1, _Result> unit() const
	{
		return
		{	_Result(1)};
	}

	operator math::matrix<1, 1, _T>&() const
	{
		return *(math::matrix<1, 1, _T>*)this;
	}

	template<typename _T2>
	inline operator vector<1, _T2>() const
	{
		return
		{	(_T2)coord[0]};
	}

	template<typename _T2>
	inline vector<1, _T>& operator=(const vector<1, _T2>& vec)
	{
		coord[0] = (_T)vec.coord[0];
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline vector<1, _Result> operator+(const vector<1, _T2>& vec) const
	{
		return
		{	coord[0] + vec.coord[0]};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline vector<1, _Result> operator-(const vector<1, _T2>& vec) const
	{
		return
		{	coord[0] - vec.coord[0]};
	}

	template<typename _T2>
	inline vector<1, _T>& operator+=(const vector<1, _T2>& vec)
	{
		coord[0] += vec.coord[0];
		return *this;
	}

	template<typename _T2>
	inline vector<1, _T>& operator-=(const vector<1, _T2>& vec)
	{
		coord[0] -= vec.coord[0];
		return *this;
	}

	inline vector<1, _T> operator-() const
	{
		return
		{	-coord[0]};
	}

	template<typename _T2>
	inline vector<1, _T>& operator*=(const _T2& t)
	{
		coord[0] *= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<1, _T2>& vec) const
	{
		return coord[0] * vec.coord[0];
	}

	template<typename _T2>
	inline vector<1, _T>& operator/=(const _T2& t)
	{
		coord[0] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<1, _Result> operator/(const vector<1, _T2> vec) const
	{
		return
		{	coord[0] / vec.coord[0]};
	}

	inline static vector<1, _T> all(_T value)
	{
		return
		{	value};
	}

	__def_vector_zero_one__(_T, 1)

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum() const
	{
		return coord[0];
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0] << ')';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<1, _Result> operator*(const vector<1, _T2>& vec, _T2 t)
{
	return
	{	vec.coord[0] * t};
}
template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<1, _Result> operator/(const vector<1, _T2>& vec, _T2 t)
{
	return
	{	vec.coord[0] / t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<1, _Result> operator*(_T1 t, const vector<1, _T2>& vec)
{
	return
	{	t * vec.coord[0]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<1, _Result> operator/(_T1 t, const vector<1, _T2>& vec)
{
	return
	{	t / vec.coord[0]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<1, _Result> hadamard(const vector<1, _T1>& vec1, const vector<1, _T2>& vec2)
{
	return
	{	vec1[0] * vec2[0]};
}

//2维向量特化
template<typename _T>
class vector<2, _T>
{
public:
	__def_vector_members__(_T, 2)

	__def_vector_op_cast__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_norm__(_T)

__def_vector_op_unit__	(_T, 2)

	__def_vector_op_to_matrix__(_T, 2)

	__def_vector_op_slice__(_T, 2)

	__def_vector_op_mask__(_T, 2)

	//不定义构造函数，使其成为聚合体，可使用{}直接指定数组成员初始化

	template<typename _T2>
	inline operator vector<2, _T2>() const
	{
		return
		{	(_T2)coord[0], (_T2)coord[1]};
	}

	template<typename _T2>
	inline vector<2, _T>& operator=(const vector<2, _T2>& vec)
	{
		coord[0] = (_T)vec.coord[0];
		coord[1] = (_T)vec.coord[1];
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline vector<2, _Result> operator+(const vector<2, _T2>& vec) const
	{
		return
		{	coord[0] + vec.coord[0], coord[1] + vec.coord[1]};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline vector<2, _Result> operator-(const vector<2, _T2>& vec) const
	{
		return
		{	coord[0] - vec.coord[0], coord[1] - vec.coord[1]};
	}

	inline vector<2, _T> operator-() const
	{
		return
		{	-coord[0], -coord[1]};
	}

	template<typename _T2>
	inline vector<2, _T>& operator+=(const vector<2, _T2>& vec)
	{
		coord[0] += vec.coord[0];
		coord[1] += vec.coord[1];
		return *this;
	}

	template<typename _T2>
	inline vector<2, _T>& operator-=(const vector<2, _T2>& vec)
	{
		coord[0] -= vec.coord[0];
		coord[1] -= vec.coord[1];
		return *this;
	}

	template<typename _T2>
	inline vector<2, _T>& operator*=(const _T2& t)
	{
		coord[0] *= t;
		coord[1] *= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<2, _T2>& vec) const
	{
		return coord[0] * vec.coord[0] + coord[1] * vec.coord[1];
	}

	template<typename _T2>
	inline vector<2, _T>& operator/=(const _T2& t)
	{
		coord[0] /= t;
		coord[1] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<2, _Result> operator/(const vector<2, _T2>& vec) const
	{
		return
		{
			coord[0] / vec.coord[0],
			coord[1] / vec.coord[1]
		};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator^(const vector<2, _T2>& vec) const
	{
		return coord[0] * vec.coord[1] - coord[1] * vec.coord[0];
	}

	inline static vector<2, _T> all(_T value)
	{
		return
		{	value, value};
	}
	__def_vector_zero_one__(_T, 2)

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum() const
	{
		return coord[0] + coord[1];
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0] << ", " << coord[1] << ')';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<2, _Result> operator*(const vector<2, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] * t,
		vec.coord[1] * t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<2, _Result> operator/(const vector<2, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] / t,
		vec.coord[1] / t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<2, _Result> operator*(_T1 t, const vector<2, _T2>& vec)
{
	return
	{	t * vec.coord[0],
		t * vec.coord[1]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<2, _Result> operator/(_T1 t, const vector<2, _T2>& vec)
{
	return
	{	t / vec.coord[0],
		t / vec.coord[1]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<2, _Result> hadamard(const vector<2, _T1>& vec1, const vector<2, _T2>& vec2)
{
	return
	{	vec1[0] * vec2[0],
		vec1[1] * vec2[1]};
}

//3维向量特化
template<typename _T>
class vector<3, _T>
{
public:
	__def_vector_members__(_T, 3)

	__def_vector_op_cast__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_norm__(_T)

__def_vector_op_unit__	(_T, 3)

	__def_vector_op_to_matrix__(_T, 3)

	__def_vector_op_slice__(_T, 3)

	__def_vector_op_mask__(_T, 3)

	template<typename _T2>
	inline operator vector<3, _T2>() const
	{
		return
		{	(_T2)coord[0],(_T2)coord[1],(_T2)coord[2]};
	}

	template<typename _T2>
	inline vector<3, _T>& operator=(const vector<3, _T2>& vec)
	{
		coord[0] = (_T)coord[0];
		coord[1] = (_T)coord[1];
		coord[2] = (_T)coord[2];
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline vector<3, _Result> operator+(const vector<3, _T2>& vec) const
	{
		return
		{	coord[0]+vec.coord[0],coord[1]+vec.coord[1],coord[2]+vec.coord[2]};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline vector<3, _Result> operator-(const vector<3, _T2>& vec) const
	{
		return
		{	coord[0]-vec.coord[0],coord[1]-vec.coord[1],coord[2]-vec.coord[2]};
	}

	template<typename _T2>
	inline vector<3, _T>& operator+=(const vector<3, _T2>& vec)
	{
		coord[0] += vec.coord[0];
		coord[1] += vec.coord[1];
		coord[2] += vec.coord[2];
		return *this;
	}

	template<typename _T2>
	inline vector<3, _T>& operator-=(const vector<3, _T2>& vec)
	{
		coord[0] -= vec.coord[0];
		coord[1] -= vec.coord[1];
		coord[2] -= vec.coord[2];
		return *this;
	}

	inline vector<3, _T> operator-() const
	{
		return
		{	-coord[0],-coord[1],-coord[2]};
	}

	template<typename _T2>
	inline vector<3, _T>& operator*=(const _T2& t)
	{
		coord[0] *= t;
		coord[1] *= t;
		coord[2] *= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<3, _T2> vec) const
	{
		return coord[0] * vec.coord[0] + coord[1] * vec.coord[1] + coord[2] * vec.coord[2];
	}

	template<typename _T2>
	inline vector<3, _T>& operator/=(const _T2& t)
	{
		coord[0] /= t;
		coord[1] /= t;
		coord[2] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<3, _Result> operator/(const vector<3, _T2>& vec) const
	{
		return
		{
			coord[0] / vec.coord[0],
			coord[1] / vec.coord[1],
			coord[2] / vec.coord[2]
		};
	}

	/**
	 * 3维向量叉乘
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<3, _Result> operator^(const vector<3, _T2>& vec) const
	{
		return
		{	coord[1]*vec.coord[2]-coord[2]*vec.coord[1],coord[2]*vec.coord[0]-coord[0]*vec.coord[2],coord[0]*vec.coord[1]-coord[1]*vec.coord[0]};
	}

	inline static vector<3, _T> all(_T value)
	{
		return
		{	value, value, value};
	}
	__def_vector_zero_one__(_T, 3)

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum() const
	{
		return coord[0] + coord[1] + coord[2];
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0]<< ", " << coord[1]<< ", " << coord[2] << ')';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<3, _Result> operator*(const vector<3, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] * t, vec.coord[1] * t, vec.coord[2] * t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<3, _Result> operator/(const vector<3, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] / t, vec.coord[1] / t, vec.coord[2] / t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<3, _Result> operator*(_T1 t, const vector<3, _T2>& vec)
{
	return
	{	t * vec.coord[0],
		t * vec.coord[1],
		t * vec.coord[2]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<3, _Result> operator/(_T1 t, const vector<3, _T2>& vec)
{
	return
	{	t / vec.coord[0],
		t / vec.coord[1],
		t / vec.coord[2]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<3, _Result> hadamard(const vector<3, _T1>& vec1, const vector<3, _T2>& vec2)
{
	return
	{	vec1[0] * vec2[0],
		vec1[1] * vec2[1],
		vec1[2] * vec2[2]};
}

//4维向量特化
template<typename _T>
class vector<4, _T>
{
public:
	__def_vector_members__(_T, 4)

	__def_vector_op_cast__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_norm__(_T)

__def_vector_op_unit__	(_T, 4)

	__def_vector_op_to_matrix__(_T, 4)

	__def_vector_op_slice__(_T, 4)

	__def_vector_op_mask__(_T, 4)

	template<typename _T2>
	inline operator vector<4, _T2>() const
	{
		return
		{	(_T2)coord[0], (_T2)coord[1], (_T2)coord[2], (_T2)coord[3]};
	}

	template<typename _T2>
	inline vector<4, _T>& operator=(const vector<4, _T2>& vec)
	{
		coord[0] = (_T)vec.coord[0];
		coord[1] = (_T)vec.coord[1];
		coord[2] = (_T)vec.coord[2];
		coord[3] = (_T)vec.coord[3];
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline vector<4, _Result> operator+(const vector<4, _T2>& vec) const
	{
		return
		{
			coord[0] + vec.coord[0],
			coord[1] + vec.coord[1],
			coord[2] + vec.coord[2],
			coord[3] + vec.coord[3]
		};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline vector<4, _Result> operator-(const vector<4, _T2>& vec) const
	{
		return
		{
			coord[0] - vec.coord[0],
			coord[1] - vec.coord[1],
			coord[2] - vec.coord[2],
			coord[3] - vec.coord[3]
		};
	}

	template<typename _T2>
	inline vector<4, _T>& operator+=(const vector<4, _T2>& vec)
	{
		coord[0] += vec.coord[0];
		coord[1] += vec.coord[1];
		coord[2] += vec.coord[2];
		coord[3] += vec.coord[3];
		return *this;
	}

	template<typename _T2>
	inline vector<4, _T>& operator-=(const vector<4, _T2>& vec)
	{
		coord[0] -= vec.coord[0];
		coord[1] -= vec.coord[1];
		coord[2] -= vec.coord[2];
		coord[3] -= vec.coord[3];
		return *this;
	}

	inline vector<4, _T> operator-() const
	{
		return
		{	-coord[0], -coord[1], -coord[2], -coord[3]};
	}

	template<typename _T2>
	inline vector<4, _T>& operator*=(const _T2& t)
	{
		coord[0] *= t;
		coord[1] *= t;
		coord[2] *= t;
		coord[3] *= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<4, _T2>& vec) const
	{
		return coord[0] * vec.coord[0] + coord[1] * vec.coord[1] + coord[2] * vec.coord[2] + coord[3] * vec.coord[3];
	}

	template<typename _T2>
	inline vector<4, _T>& operator/=(const _T2& t)
	{
		coord[0] /= t;
		coord[1] /= t;
		coord[2] /= t;
		coord[3] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<4, _Result> operator/(const vector<4, _T2>& vec) const
	{
		return
		{
			coord[0] / vec.coord[0],
			coord[1] / vec.coord[1],
			coord[2] / vec.coord[2],
			coord[3] / vec.coord[3]
		};
	}

	inline static vector<4, _T> all(_T value)
	{
		return
		{	value, value, value, value};
	}
	__def_vector_zero_one__(_T, 4)

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum() const
	{
		return coord[0] + coord[1] + coord[2] + coord[3];
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0] << ", " << coord[1] << ", " << coord[2] << ", " << coord[3] << ')';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<4, _Result> operator*(const vector<4, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] * t, vec.coord[1] * t, vec.coord[2] * t, vec.coord[3] * t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<4, _Result> operator/(const vector<4, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] / t, vec.coord[1] / t, vec.coord[2] / t, vec.coord[3] / t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<4, _Result> operator*(_T1 t, const vector<4, _T2>& vec)
{
	return
	{	t * vec.coord[0],
		t * vec.coord[1],
		t * vec.coord[2],
		t * vec.coord[3]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<4, _Result> operator/(_T1 t, const vector<4, _T2>& vec)
{
	return
	{	t / vec.coord[0],
		t / vec.coord[1],
		t / vec.coord[2],
		t / vec.coord[3]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<4, _Result> hadamard(const vector<4, _T1>& vec1, const vector<4, _T2>& vec2)
{
	return
	{	vec1[0] * vec2[0],
		vec1[1] * vec2[1],
		vec1[2] * vec2[2],
		vec1[3] * vec2[3]};
}

//6维向量特化
template<typename _T>
class vector<6, _T>
{
public:
	__def_vector_members__(_T, 6)

	__def_vector_op_cast__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_norm__(_T)

__def_vector_op_unit__	(_T, 6)

	__def_vector_op_to_matrix__(_T, 6)

	__def_vector_op_slice__(_T, 6)

	__def_vector_op_mask__(_T, 6)

	template<typename _T2>
	inline operator vector<6, _T2>() const
	{
		return
		{	(_T2)coord[0], (_T2)coord[1], (_T2)coord[2],
			(_T2)coord[3], (_T2)coord[4], (_T2)coord[5]};
	}

	template<typename _T2>
	inline vector<6, _T>& operator=(const vector<6, _T2>& vec)
	{
		coord[0] = (_T)vec.coord[0];
		coord[1] = (_T)vec.coord[1];
		coord[2] = (_T)vec.coord[2];
		coord[3] = (_T)vec.coord[3];
		coord[4] = (_T)vec.coord[4];
		coord[5] = (_T)vec.coord[5];
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline vector<6, _Result> operator+(const vector<6, _T2>& vec) const
	{
		return
		{
			coord[0] + vec.coord[0],
			coord[1] + vec.coord[1],
			coord[2] + vec.coord[2],
			coord[3] + vec.coord[3],
			coord[4] + vec.coord[4],
			coord[5] + vec.coord[5]
		};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline vector<6, _Result> operator-(const vector<6, _T2>& vec) const
	{
		return
		{
			coord[0] - vec.coord[0],
			coord[1] - vec.coord[1],
			coord[2] - vec.coord[2],
			coord[3] - vec.coord[3],
			coord[4] - vec.coord[4],
			coord[5] - vec.coord[5]
		};
	}

	template<typename _T2>
	inline vector<6, _T>& operator+=(const vector<6, _T2>& vec)
	{
		coord[0] += vec.coord[0];
		coord[1] += vec.coord[1];
		coord[2] += vec.coord[2];
		coord[3] += vec.coord[3];
		coord[4] += vec.coord[4];
		coord[5] += vec.coord[5];
		return *this;
	}

	template<typename _T2>
	inline vector<6, _T>& operator-=(const vector<6, _T2>& vec)
	{
		coord[0] -= vec.coord[0];
		coord[1] -= vec.coord[1];
		coord[2] -= vec.coord[2];
		coord[3] -= vec.coord[3];
		coord[4] -= vec.coord[4];
		coord[5] -= vec.coord[5];
		return *this;
	}

	inline vector<6, _T> operator-() const
	{
		return
		{	-coord[0], -coord[1], -coord[2],
			-coord[3], -coord[4], -coord[5]};
	}

	template<typename _T2>
	inline vector<6, _T>& operator*=(const _T2& t)
	{
		coord[0] *= t;
		coord[1] *= t;
		coord[2] *= t;
		coord[3] *= t;
		coord[4] *= t;
		coord[5] *= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<6, _T2>& vec) const
	{
		return coord[0] * vec.coord[0] + coord[1] * vec.coord[1] +
		coord[2] * vec.coord[2] + coord[3] * vec.coord[3] +
		coord[4] * vec.coord[4] + coord[5] * vec.coord[5];
	}

	template<typename _T2>
	inline vector<6, _T>& operator/=(const _T2& t)
	{
		coord[0] /= t;
		coord[1] /= t;
		coord[2] /= t;
		coord[3] /= t;
		coord[4] /= t;
		coord[5] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<6, _Result> operator/(const vector<6, _T2>& vec) const
	{
		return
		{
			coord[0] / vec.coord[0],
			coord[1] / vec.coord[1],
			coord[2] / vec.coord[2],
			coord[3] / vec.coord[3],
			coord[4] / vec.coord[4],
			coord[5] / vec.coord[5]
		};
	}

	inline static vector<6, _T> all(_T value)
	{
		return
		{	value, value, value, value, value, value};
	}
	__def_vector_zero_one__(_T, 6)

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum() const
	{
		return coord[0] + coord[1] + coord[2] + coord[3] + coord[4] + coord[5];
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0] << ", " << coord[1] << ", " << coord[2]
		<< ", " << coord[3] << ", " << coord[4] << ", " << coord[5] << ')';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<6, _Result> operator*(const vector<6, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] * t, vec.coord[1] * t, vec.coord[2] * t,
		vec.coord[3] * t, vec.coord[4] * t, vec.coord[5] * t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<6, _Result> operator/(const vector<6, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] / t, vec.coord[1] / t, vec.coord[2] / t,
		vec.coord[3] / t, vec.coord[4] / t, vec.coord[5] / t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<6, _Result> operator*(_T1 t, const vector<6, _T2>& vec)
{
	return
	{	t * vec.coord[0],
		t * vec.coord[1],
		t * vec.coord[2],
		t * vec.coord[3],
		t * vec.coord[4],
		t * vec.coord[5]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<6, _Result> operator/(_T1 t, const vector<6, _T2>& vec)
{
	return
	{	t / vec.coord[0],
		t / vec.coord[1],
		t / vec.coord[2],
		t / vec.coord[3],
		t / vec.coord[4],
		t / vec.coord[5]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<6, _Result> hadamard(const vector<6, _T1>& vec1, const vector<6, _T2>& vec2)
{
	return
	{	vec1[0] * vec2[0],
		vec1[1] * vec2[1],
		vec1[2] * vec2[2],
		vec1[3] * vec2[3],
		vec1[4] * vec2[4],
		vec1[5] * vec2[5]};
}

//7维向量特化
template<typename _T>
class vector<7, _T>
{
public:
	__def_vector_members__(_T, 7)

	__def_vector_op_cast__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_norm__(_T)

__def_vector_op_unit__	(_T, 7)

	__def_vector_op_to_matrix__(_T, 7)

	__def_vector_op_slice__(_T, 7)

	__def_vector_op_mask__(_T, 7)

	template<typename _T2>
	inline operator vector<7, _T2>() const
	{
		return
		{	(_T2)coord[0], (_T2)coord[1], (_T2)coord[2],
			(_T2)coord[3], (_T2)coord[4], (_T2)coord[5],
			(_T2)coord[6]};
	}

	template<typename _T2>
	inline vector<7, _T>& operator=(const vector<7, _T2>& vec)
	{
		coord[0] = (_T)vec.coord[0];
		coord[1] = (_T)vec.coord[1];
		coord[2] = (_T)vec.coord[2];
		coord[3] = (_T)vec.coord[3];
		coord[4] = (_T)vec.coord[4];
		coord[5] = (_T)vec.coord[5];
		coord[6] = (_T)vec.coord[6];
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline vector<7, _Result> operator+(const vector<7, _T2>& vec) const
	{
		return
		{
			coord[0] + vec.coord[0],
			coord[1] + vec.coord[1],
			coord[2] + vec.coord[2],
			coord[3] + vec.coord[3],
			coord[4] + vec.coord[4],
			coord[5] + vec.coord[5],
			coord[6] + vec.coord[6]
		};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline vector<7, _Result> operator-(const vector<7, _T2>& vec) const
	{
		return
		{
			coord[0] - vec.coord[0],
			coord[1] - vec.coord[1],
			coord[2] - vec.coord[2],
			coord[3] - vec.coord[3],
			coord[4] - vec.coord[4],
			coord[5] - vec.coord[5],
			coord[6] - vec.coord[6]
		};
	}

	template<typename _T2>
	inline vector<7, _T>& operator+=(const vector<7, _T2>& vec)
	{
		coord[0] += vec.coord[0];
		coord[1] += vec.coord[1];
		coord[2] += vec.coord[2];
		coord[3] += vec.coord[3];
		coord[4] += vec.coord[4];
		coord[5] += vec.coord[5];
		coord[6] += vec.coord[6];
		return *this;
	}

	template<typename _T2>
	inline vector<7, _T>& operator-=(const vector<7, _T2>& vec)
	{
		coord[0] -= vec.coord[0];
		coord[1] -= vec.coord[1];
		coord[2] -= vec.coord[2];
		coord[3] -= vec.coord[3];
		coord[4] -= vec.coord[4];
		coord[5] -= vec.coord[5];
		coord[6] -= vec.coord[6];
		return *this;
	}

	inline vector<7, _T> operator-() const
	{
		return
		{	-coord[0], -coord[1], -coord[2],
			-coord[3], -coord[4], -coord[5],
			-coord[6]};
	}

	template<typename _T2>
	inline vector<7, _T>& operator*=(const _T2& t)
	{
		coord[0] *= t;
		coord[1] *= t;
		coord[2] *= t;
		coord[3] *= t;
		coord[4] *= t;
		coord[5] *= t;
		coord[6] *= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<7, _T2>& vec) const
	{
		return coord[0] * vec.coord[0] + coord[1] * vec.coord[1] +
		coord[2] * vec.coord[2] + coord[3] * vec.coord[3] +
		coord[4] * vec.coord[4] + coord[5] * vec.coord[5] +
		coord[6] * vec.coord[6];
	}

	template<typename _T2>
	inline vector<7, _T>& operator/=(const _T2& t)
	{
		coord[0] /= t;
		coord[1] /= t;
		coord[2] /= t;
		coord[3] /= t;
		coord[4] /= t;
		coord[5] /= t;
		coord[6] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<7, _Result> operator/(const vector<7, _T2>& vec) const
	{
		return
		{
			coord[0] / vec.coord[0],
			coord[1] / vec.coord[1],
			coord[2] / vec.coord[2],
			coord[3] / vec.coord[3],
			coord[4] / vec.coord[4],
			coord[5] / vec.coord[5],
			coord[6] / vec.coord[6]
		};
	}

	inline static vector<7, _T> all(_T value)
	{
		return
		{	value, value, value, value, value, value, value};
	}
	__def_vector_zero_one__(_T, 7)

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum() const
	{
		return coord[0] + coord[1] + coord[2] + coord[3] + coord[4] + coord[5] + coord[6];
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0] << ", " << coord[1] << ", " << coord[2]
		<< ", " << coord[3] << ", " << coord[4] << ", " << coord[5]
		<< ", " << coord[6] << ')';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<7, _Result> operator*(const vector<7, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] * t, vec.coord[1] * t, vec.coord[2] * t,
		vec.coord[3] * t, vec.coord[4] * t, vec.coord[5] * t,
		vec.coord[6] * t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<7, _Result> operator/(const vector<7, _T1>& vec, _T2 t)
{
	return
	{	vec.coord[0] / t, vec.coord[1] / t, vec.coord[2] / t,
		vec.coord[3] / t, vec.coord[4] / t, vec.coord[5] / t,
		vec.coord[6] / t};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<7, _Result> operator*(_T1 t, const vector<7, _T2>& vec)
{
	return
	{	t * vec.coord[0],
		t * vec.coord[1],
		t * vec.coord[2],
		t * vec.coord[3],
		t * vec.coord[4],
		t * vec.coord[5],
		t * vec.coord[6]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<7, _Result> operator/(_T1 t, const vector<7, _T2>& vec)
{
	return
	{	t / vec.coord[0],
		t / vec.coord[1],
		t / vec.coord[2],
		t / vec.coord[3],
		t / vec.coord[4],
		t / vec.coord[5],
		t / vec.coord[6]};
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline vector<7, _Result> hadamard(const vector<7, _T1>& vec1, const vector<7, _T2>& vec2)
{
	return
	{	vec1[0] * vec2[0],
		vec1[1] * vec2[1],
		vec1[2] * vec2[2],
		vec1[3] * vec2[3],
		vec1[4] * vec2[4],
		vec1[5] * vec2[5],
		vec1[6] * vec2[6]
	};
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
using vector6 = vector<6, _T>;
template<typename _T>
using vector7 = vector<7, _T>;

using vector1f = vector1<float>;
using vector2f = vector2<float>;
using vector3f = vector3<float>;
using vector4f = vector4<float>;
using vector6f = vector6<float>;
using vector7f = vector7<float>;

using vector1d = vector1<double>;
using vector2d = vector2<double>;
using vector3d = vector3<double>;
using vector4d = vector4<double>;
using vector6d = vector6<double>;
using vector7d = vector7<double>;

using vector1q = vector1<__float128>;
using vector2q = vector2<__float128>;
using vector3q = vector3<__float128>;
using vector4q = vector4<__float128>;
using vector6q = vector6<__float128>;
using vector7q = vector7<__float128>;
}

#endif//_MATH_VECTOR