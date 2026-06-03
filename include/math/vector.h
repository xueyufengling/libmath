#ifndef _MATH_VECTOR
#define _MATH_VECTOR

#include <tplmp/tplmp.h>

#include <string.h>
#include <sstream>

#include <math/op.h>

namespace math
{
/**
 * 列向量，初始化使用vector vec={x,y,z,w...};
 */
template<size_t _Dim, typename _T = float>
class vector
{
public:

#define __def_vector_members__(coord_type, vec_dim)\
	coord_type coord[vec_dim];\
	typedef coord_type type;\
	static const size_t dim = vec_dim;
	__def_vector_members__(_T, _Dim)

#define __def_vector_op_cast_ptr__(coord_type)\
	inline operator coord_type*()\
	{\
		return (coord_type*)coord;\
	}
	__def_vector_op_cast_ptr__(_T)

	/**
	 * 坐标
	 */
#define __def_vector_op_coord__(coord_type)\
	inline coord_type& operator[](int i)\
	{\
		return coord[i];\
	}\
	inline const coord_type& operator[](int i) const\
	{\
		return coord[i];\
	}
	__def_vector_op_coord__(_T)

	/**
	 * 向量长度
	 */
#define __def_vector_op_len__(coord_type)\
	template<typename _Result = decltype(tplmp::decl<coord_type>::val() * tplmp::decl<coord_type>::val())>\
	inline _Result length() const\
	{\
		return math::sqrt<_Result>(((*this)) * (*this));\
	}
	__def_vector_op_len__(_T)

	/**
	 * 单位向量
	 */
#define __def_vector_op_unit__(coord_type, vec_dim)\
	template<typename _Result = decltype(tplmp::decl<coord_type>::val() / tplmp::decl<coord_type>::val())>\
	inline vector<vec_dim, _Result> unit() const\
	{\
		return (*this) / length();\
	}
__def_vector_op_unit__	(_T, _Dim)

	/**
	 * 类型强制转换
	 */
	template<typename _T2>
	inline operator vector<_Dim,_T2>() const
	{
		vector<_Dim, _T2> cast_result;
		for(int i = 0; i < _Dim; ++i)
		cast_result[i] = (_T2)coord[i];
		return cast_result;
	}

	/**
	 * 强制转换赋值
	 */
	template<typename _T2>
	inline vector<_Dim, _T>& operator=(const vector<_Dim, _T2>& vec)
	{
		for(int i = 0; i < _Dim; ++i)
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
		for(int i = 0; i < _Dim; ++i)
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
		for(int i = 0; i < _Dim; ++i)
		sub_result[i] = (_Result)(coord[i] - vec.coord[i]);
		return sub_result;
	}

	/**
	 * 向量取反
	 */
	inline vector<_Dim, _T> operator-() const
	{
		vector<_Dim, _T> inv;
		for(int i = 0; i < _Dim; ++i)
		inv[i] = -coord[i];
		return inv;
	}

	/**
	 * 标量乘法
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<_Dim, _Result> operator*(_T2 t) const
	{
		vector<_Dim, _Result> mul_result;
		for(size_t i = 0; i < _Dim; ++i)
		mul_result[i] = (_Result)(coord[i] * t);
		return mul_result;
	}

	/**
	 * 标量除法
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<_Dim, _Result> operator/(_T2 t) const
	{
		vector<_Dim, _Result> div_result;
		for(size_t i = 0; i < _Dim; ++i)
		div_result[i] = (_Result)(coord[i] / t);
		return div_result;
	}

	/**
	 * 向量点乘
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<_Dim, _T2> vec) const
	{
		_Result dot_result = 0;
		for(size_t i = 0; i < _Dim; ++i)
		dot_result += coord[i] * vec.coord[i];
		return dot_result;
	}

	/**
	 * 零向量
	 */
	inline static vector<_Dim, _T> zero()
	{
		vector<_Dim, _T> zero_vec;
		memset(&zero_vec.coord, 0, _Dim * sizeof(_T));
		return zero_vec;
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

template<size_t _Dim, typename _T>
inline std::ostream& operator<<(std::ostream& os, const vector<_Dim, _T>& vec)
{
	os << ((std::string)vec);
	return os;
}

//交换律
template<typename _T1, size_t _Dim, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline vector<_Dim, _Result> operator*(_T1 t, const vector<_Dim, _T2>& vec)
{
	return vec * t;
}

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

// 2维向量特化
template<typename _T>
class vector<2, _T>
{
public:
	__def_vector_members__(_T, 2)

	__def_vector_op_cast_ptr__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_len__(_T)

__def_vector_op_unit__	(_T, 2)

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

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<2, _Result> operator*(_T2 t) const
	{
		return
		{	coord[0] * t, coord[1] * t};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<2, _T2>& vec) const
	{
		return coord[0] * vec.coord[0] + coord[1] * vec.coord[1];
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<2, _Result> operator/(_T2 t) const
	{
		return
		{	coord[0] / t, coord[1] / t};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator^(const vector<2, _T2>& vec) const
	{
		return coord[0] * vec.coord[1] - coord[1] * vec.coord[0];
	}

	inline static vector<2, _T> zero()
	{
		return
		{	0, 0};
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0] << ", " << coord[1] << ')';
		return oss.str();
	}
};

// 3维向量特化
template<typename _T>
class vector<3, _T>
{
public:
	__def_vector_members__(_T, 3)

	__def_vector_op_cast_ptr__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_len__(_T)

__def_vector_op_unit__	(_T, 3)

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

	inline vector<3, _T> operator-() const
	{
		return
		{	-coord[0],-coord[1],-coord[2]};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<3, _Result> operator*(_T2 t) const
	{
		return
		{	coord[0]*t,coord[1]*t,coord[2]*t};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<3, _T2> vec) const
	{
		return coord[0] * vec.coord[0] + coord[1] * vec.coord[1] + coord[2] * vec.coord[2];
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<3, _Result> operator/(_T2 t) const
	{
		return
		{	coord[0]/t,coord[1]/t,coord[2]/t};
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

	inline static vector<3, _T> zero()
	{
		return
		{	0,0,0};
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0]<< ", " << coord[1]<< ", " << coord[2] << ')';
		return oss.str();
	}
};

// 4维向量特化
template<typename _T>
class vector<4, _T>
{
public:
	__def_vector_members__(_T, 4)

	__def_vector_op_cast_ptr__(_T)

	__def_vector_op_coord__(_T)

	__def_vector_op_len__(_T)

__def_vector_op_unit__	(_T, 4)

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

	inline vector<4, _T> operator-() const
	{
		return
		{	-coord[0], -coord[1], -coord[2], -coord[3]};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<4, _Result> operator*(_T2 t) const
	{
		return
		{	coord[0] * t, coord[1] * t, coord[2] * t, coord[3] * t};
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline _Result operator*(const vector<4, _T2>& vec) const
	{
		return coord[0] * vec.coord[0] + coord[1] * vec.coord[1] + coord[2] * vec.coord[2] + coord[3] * vec.coord[3];
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline vector<4, _Result> operator/(_T2 t) const
	{
		return
		{	coord[0] / t, coord[1] / t, coord[2] / t, coord[3] / t};
	}

	inline static vector<4, _T> zero()
	{
		return
		{	0, 0, 0, 0};
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '(' << coord[0] << ", " << coord[1] << ", " << coord[2] << ", " << coord[3] << ')';
		return oss.str();
	}
};

template<typename _T>
using vector2 = vector<2, _T>;
template<typename _T>
using vector3 = vector<3, _T>;
template<typename _T>
using vector4 = vector<4, _T>;

using vector2f = vector2<float>;
using vector3f = vector3<float>;
using vector4f = vector4<float>;

using vector2d = vector2<double>;
using vector3d = vector3<double>;
using vector4d = vector4<double>;

}

#endif//_MATH_VECTOR