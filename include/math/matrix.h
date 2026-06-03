#ifndef _MATH_MATRIX
#define _MATH_MATRIX

#include <tplmp/tplmp.h>

#include <string.h>
#include <sstream>

#include <math/vector.h>

namespace math
{
/**
 * 矩阵。
 * 采用行主序以贴合矩阵乘向量的数学含义
 */
template<size_t _Row, size_t _Column, typename _T = float>
class matrix
{
public:
	//索引行、列。储存Row个Column维行向量
#define __def_matrix_members__(elem_type, row, column)\
	math::vector<column, elem_type> elem[row];\
	typedef elem_type type;\
	static const size_t row_num = row;\
	static const size_t column_num = column;
	__def_matrix_members__(_T, _Row, _Column)

#define __def_matrix_op_cast_ptr__(elem_type)\
	inline operator elem_type*()\
	{\
		return (elem_type*)elem;\
	}
	__def_matrix_op_cast_ptr__(_T)

#define __def_matrix_op_cast_vec_ptr__(elem_type, column)\
	inline operator math::vector<column, elem_type>*()\
	{\
		return (math::vector<column, elem_type>*)elem;\
	}
	__def_matrix_op_cast_vec_ptr__(_T, _Column)

	/**
	 * 提取行向量
	 */
#define __def_matrix_op_row_vec__(elem_type, column)\
	inline math::vector<column, elem_type>& operator[](size_t row)\
	{\
		return elem[row];\
	}\
	inline const math::vector<column, elem_type>& operator[](size_t row) const\
	{\
		return elem[row];\
	}\
	inline math::vector<column, _T> row(int row)\
	{\
		return elem[row];\
	}
	__def_matrix_op_row_vec__(_T, _Column)

	/**
	 * 提取列向量
	 */
	inline math::vector<_Row, _T> column(int column)
	{
		vector<_Row, _T> row_vec;
		for(size_t i = 0; i < _Row; ++i)
			row_vec[i] = elem[i][column];
		return row_vec;
	}

	template<typename _T2>
	inline operator matrix<_Row, _Column , _T2>() const
	{
		matrix<_Row, _Column, _T2> cast_result;
		for(size_t i = 0; i < _Row; ++i)
			cast_result[i] = (_T2)elem[i];
		return cast_result;
	}

	template<typename _T2>
	inline matrix<_Row, _Column, _T>& operator=(const matrix<_Row, _Column, _T2>& mat)
	{
		for(size_t i = 0; i < _Row; ++i)
			elem[i] = mat.elem[i];
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator+(const matrix<_Row, _Column, _T2>& mat) const
	{
		matrix<_Row, _Column, _Result> add_result;
		for(size_t i = 0; i < _Row; ++i)
			add_result[i] = elem[i] + mat.elem[i];
		return add_result;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator-(const matrix<_Row, _Column, _T2>& mat) const
	{
		matrix<_Row, _Column, decltype(tplmp::decl<_T>::val()-tplmp::decl<_T2>::val())> sub_result;
		for(size_t i = 0; i < _Row; ++i)
			sub_result[i] = elem[i] - mat.elem[i];
		return sub_result;
	}

	/**
	 * 数乘
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator*(_T2 t) const
	{
		matrix<_Row, _Column, _Result> mul_result;
		for(size_t i = 0; i < _Row; ++i)
			mul_result[i] = elem[i] * t;
		return mul_result;
	}

	template<size_t _Row2, typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Row2, _Result> operator*(const matrix<_Column, _Row2, _T2>& mat) const
	{
		matrix<_Row, _Row2, _Result> mul_result;
		for(size_t i = 0; i < _Row; ++i)
		{
			for(size_t j = 0; j < _Row2; ++j)
				mul_result[i][j] = row(i) * mat.column(j);
		}
		return mul_result;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline math::vector<_Row, _Result> operator*(const math::vector<_Column, _T2>& vec) const
	{
		vector<_Row, _Result> mul_result;
		for(size_t i = 0; i < _Row; ++i)
			mul_result[i] = row(i) * vec;
		return mul_result;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator/(_T2 t) const
	{
		matrix<_Row, _Column, _Result> div_result;
		for(size_t i = 0; i < _Row; ++i)
			div_result[i] = elem[i] / t;
		return div_result;
	}

	/**
	 * 转置
	 */
	inline matrix<_Column, _Row, _T> operator~() const
	{
		matrix<_Column, _Row, _T> transpose_result;
		for(size_t i = 0; i < _Column; ++i)
			for(size_t j = 0; j < _Row; ++j)
				transpose_result[i][j] = elem[j][i];
		return transpose_result;
	}

	inline static matrix<_Row, _Column, _T> zero()
	{
		matrix<_Row, _Column, _T> zero_mat;
		memset(&zero_mat.elem, 0, _Row * _Column * sizeof(_T));
		return zero_mat;
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		for(size_t i = 0; i < _Row; ++i)
		{
			oss << '[';
			for(size_t j = 0; j < _Column; ++j)
			{
				oss << elem[i][j];
				if(j < _Column - 1)
					oss << ' ';
			}
			oss << ']';
			if(i < _Row - 1)
				oss << '\n';
		}
		return oss.str();
	}
};

template<size_t _Row, size_t _Column, typename _T>
inline std::ostream& operator<<(std::ostream& os, const matrix<_Row, _Column, _T>& mat)
{
	os << ((std::string)mat);
	return os;
}

//交换律
template<typename _T1, size_t _Row, size_t _Column, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<_Row, _Column, _Result> operator*(_T1 t, const matrix<_Row, _Column, _T2>& mat)
{
	return mat * t;
}

/**
 * 按位相乘的哈达玛积
 */
template<size_t _Row, size_t _Column, typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<_Row, _Column, _Result> hadamard(const matrix<_Row, _Column, _T1>& mat1, const matrix<_Row, _Column, _T2>& mat2)
{
	matrix<_Row, _Column, _Result> hadamard_result;
	for(size_t i = 0; i < _Row; ++i)
		for(size_t j = 0; j < _Column; ++j)
			hadamard_result[i][j] = mat1.elem[i][j] * mat2.elem[i][j];
	return hadamard_result;
}

//方阵特有

/**
 * 迹
 */
template<size_t _Order, typename _T, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
inline _Result tr(const matrix<_Order, _Order, _T>& mat)
{
	_Result tr = 0;
	for(size_t i = 0; i < _Order; ++i)
		tr += mat.elem[i][i];
	return tr;
}

/**
 * 恒等元
 */
template<size_t _Order, typename _T>
struct __identity_impl<matrix<_Order, _Order, _T> >
{
	static matrix<_Order, _Order, _T> value()
	{
		matrix<_Order, _Order, _T> identity = matrix<_Order, _Order, _T>::zero();
		for(size_t i = 0; i < _Order; ++i)
			identity[i][i] = 1;
		return identity;
	}
};

//2x2方阵特化
template<typename _T>
class matrix<2, 2, _T>
{
public:
	__def_matrix_members__(_T, 2, 2)

	__def_matrix_op_cast_ptr__(_T)

	__def_matrix_op_cast_vec_ptr__(_T, 2)

	__def_matrix_op_row_vec__(_T, 2)

	template<typename _T2>
	inline operator matrix<2, 2, _T2>()
	{
		_T2 m_elem[] =
				{
						(_T2)elem[0][0], (_T2)elem[0][1],
						(_T2)elem[1][0], (_T2)elem[1][1]
				};
		return tplmp::cast<matrix<2, 2, _T2> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<2, 2, _Result> operator+(const matrix<2, 2, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] + mat.elem[0][0], elem[0][1] + mat.elem[0][1],
						elem[1][0] + mat.elem[1][0], elem[1][1] + mat.elem[1][1]
				};
		return tplmp::cast<matrix<2, 2, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<2, 2, _Result> operator-(const matrix<2, 2, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] - mat.elem[0][0], elem[0][1] - mat.elem[0][1],
						elem[1][0] - mat.elem[1][0], elem[1][1] - mat.elem[1][1]
				};
		return tplmp::cast<matrix<2, 2, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<2, 2, _Result> operator*(_T2 t) const
	{
		_Result m_elem[] =
				{
						elem[0][0] * t, elem[0][1] * t,
						elem[1][0] * t, elem[1][1] * t
				};
		return tplmp::cast<matrix<2, 2, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<2, 2, _Result> operator*(const matrix<2, 2, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] * mat.elem[0][0] + elem[0][1] * mat.elem[1][0],
						elem[0][0] * mat.elem[0][1] + elem[0][1] * mat.elem[1][1],

						elem[1][0] * mat.elem[0][0] + elem[1][1] * mat.elem[1][0],
						elem[1][0] * mat.elem[0][1] + elem[1][1] * mat.elem[1][1]
				};
		return tplmp::cast<matrix<2, 2, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<2, _Result> operator*(const vector<2, _T2>& vec) const
	{
		_Result v_elem[] =
				{
						elem[0][0] * vec[0] + elem[0][1] * vec[1],
						elem[1][0] * vec[0] + elem[1][1] * vec[1]
				};
		return tplmp::cast<vector<2, _Result> >(v_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<2, 2, _Result> operator/(_T2 t) const
	{
		_Result m_elem[] =
				{
						elem[0][0] / t, elem[0][1] / t,
						elem[1][0] / t, elem[1][1] / t
				};
		return tplmp::cast<matrix<2, 2, _Result> >(m_elem);
	}

	inline vector<2, _T> column(int column) const
	{
		_T v_elem[] =
				{
						elem[0][column], elem[1][column]
				};
		return tplmp::cast<vector<2, _T> >(v_elem);
	}

	inline matrix<2, 2, _T> operator~() const
	{
		_T m_elem[] =
				{
						elem[0][0], elem[1][0],
						elem[0][1], elem[1][1]
				};
		return tplmp::cast<matrix<2, 2, _T> >(m_elem);
	}

	inline static matrix<2, 2, _T> zero()
	{
		_T m_elem[] =
				{
						0, 0,
						0, 0
				};
		return tplmp::cast<matrix<2, 2, _T> >(m_elem);
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << ' ' << elem[0][1] << "]\n"
				<< '[' << elem[1][0] << ' ' << elem[1][1] << ']';
		return oss.str();
	}
};

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const matrix<2, 2, _T>& mat)
{
	return mat.elem[0][0] * mat.elem[1][1] - mat.elem[0][1] * mat.elem[1][0];
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<2, 2, _Result> hadamard(const matrix<2, 2, _T1>& mat1, const matrix<2, 2, _T2>& mat2)
{
	_Result m_elem[] =
			{
					mat1.elem[0][0] * mat2.elem[0][0], mat1.elem[0][1] * mat2.elem[0][1],
					mat1.elem[1][0] * mat2.elem[1][0], mat1.elem[1][1] * mat2.elem[1][1]
			};
	return tplmp::cast<matrix<2, 2, _Result> >(m_elem);
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
inline _Result tr(const matrix<2, 2, _T>& mat)
{
	return mat.elem[0][0] + mat.elem[1][1];
}

template<typename _T>
struct __identity_impl<matrix<2, 2, _T>>
{
	static matrix<2, 2, _T> value()
	{
		_T m_elem[] =
				{
						1, 0,
						0, 1
				};
		return tplmp::cast<matrix<2, 2, _T> >(m_elem);
	}
};

//3x3方阵特化
template<typename _T>
class matrix<3, 3, _T>
{
public:
	__def_matrix_members__(_T, 3, 3)

	__def_matrix_op_cast_ptr__(_T)

	__def_matrix_op_cast_vec_ptr__(_T, 3)

	__def_matrix_op_row_vec__(_T, 3)

	template<typename _T2>
	inline operator matrix<3, 3, _T2>()
	{
		_T2 m_elem[] =
				{
						(_T2)elem[0][0], (_T2)elem[0][1], (_T2)elem[0][2],
						(_T2)elem[1][0], (_T2)elem[1][1], (_T2)elem[1][2],
						(_T2)elem[2][0], (_T2)elem[2][1], (_T2)elem[2][2]
				};
		return tplmp::cast<matrix<3, 3, _T2> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<3, 3, _Result> operator+(const matrix<3, 3, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] + mat.elem[0][0], elem[0][1] + mat.elem[0][1], elem[0][2] + mat.elem[0][2],
						elem[1][0] + mat.elem[1][0], elem[1][1] + mat.elem[1][1], elem[1][2] + mat.elem[1][2],
						elem[2][0] + mat.elem[2][0], elem[2][1] + mat.elem[2][1], elem[2][2] + mat.elem[2][2]
				};
		return tplmp::cast<matrix<3, 3, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<3, 3, _Result> operator-(const matrix<3, 3, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] - mat.elem[0][0], elem[0][1] - mat.elem[0][1], elem[0][2] - mat.elem[0][2],
						elem[1][0] - mat.elem[1][0], elem[1][1] - mat.elem[1][1], elem[1][2] - mat.elem[1][2],
						elem[2][0] - mat.elem[2][0], elem[2][1] - mat.elem[2][1], elem[2][2] - mat.elem[2][2]
				};
		return tplmp::cast<matrix<3, 3, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<3, 3, _Result> operator*(_T2 t)
	{
		_Result m_elem[] =
				{
						elem[0][0] * t, elem[0][1] * t, elem[0][2] * t,
						elem[1][0] * t, elem[1][1] * t, elem[1][2] * t,
						elem[2][0] * t, elem[2][1] * t, elem[2][2] * t
				};
		return tplmp::cast<matrix<3, 3, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<3, 3, _Result> operator*(const matrix<3, 3, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] * mat.elem[0][0] + elem[0][1] * mat.elem[1][0] + elem[0][2] * mat.elem[2][0],
						elem[0][0] * mat.elem[0][1] + elem[0][1] * mat.elem[1][1] + elem[0][2] * mat.elem[2][1],
						elem[0][0] * mat.elem[0][2] + elem[0][1] * mat.elem[1][2] + elem[0][2] * mat.elem[2][2],

						elem[1][0] * mat.elem[0][0] + elem[1][1] * mat.elem[1][0] + elem[1][2] * mat.elem[2][0],
						elem[1][0] * mat.elem[0][1] + elem[1][1] * mat.elem[1][1] + elem[1][2] * mat.elem[2][1],
						elem[1][0] * mat.elem[0][2] + elem[1][1] * mat.elem[1][2] + elem[1][2] * mat.elem[2][2],

						elem[2][0] * mat.elem[0][0] + elem[2][1] * mat.elem[1][0] + elem[2][2] * mat.elem[2][0],
						elem[2][0] * mat.elem[0][1] + elem[2][1] * mat.elem[1][1] + elem[2][2] * mat.elem[2][1],
						elem[2][0] * mat.elem[0][2] + elem[2][1] * mat.elem[1][2] + elem[2][2] * mat.elem[2][2]
				};
		return tplmp::cast<matrix<3, 3, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<3, _Result> operator*(const vector<3, _T2>& vec) const
	{
		_Result m_elem[] =
				{
						elem[0][0] * vec[0] + elem[0][1] * vec[1] + elem[0][2] * vec[2],
						elem[1][0] * vec[0] + elem[1][1] * vec[1] + elem[1][2] * vec[2],
						elem[2][0] * vec[0] + elem[2][1] * vec[1] + elem[2][2] * vec[2]
				};
		return tplmp::cast<vector<3, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<3, 3, _Result> operator/(_T2 t)
	{
		_Result m_elem[] =
				{
						elem[0][0] / t, elem[0][1] / t, elem[0][2] / t,
						elem[1][0] / t, elem[1][1] / t, elem[1][2] / t,
						elem[2][0] / t, elem[2][1] / t, elem[2][2] / t
				};
		return tplmp::cast<matrix<3, 3, _Result> >(m_elem);
	}

	inline vector<3, _T> column(int column)
	{
		_T v_elem[] =
				{
						elem[0][column], elem[1][column], elem[2][column]
				};
		return tplmp::cast<vector<3, _T> >(v_elem);
	}

	inline matrix<3, 3, _T> operator~() const
	{
		_T m_elem[] =
				{
						elem[0][0], elem[1][0], elem[2][0],
						elem[0][1], elem[1][1], elem[2][1],
						elem[0][2], elem[1][2], elem[2][2]
				};
		return tplmp::cast<matrix<3, 3, _T> >(m_elem);
	}

	inline static matrix<3, 3, _T> zero()
	{
		_T m_elem[] =
				{
						0, 0, 0,
						0, 0, 0,
						0, 0, 0
				};
		return tplmp::cast<matrix<3, 3, _T> >(m_elem);
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << ' ' << elem[0][1] << ' ' << elem[0][2] << "]\n"
				<< '[' << elem[1][0] << ' ' << elem[1][1] << ' ' << elem[1][2] << "]\n"
				<< '[' << elem[2][0] << ' ' << elem[2][1] << ' ' << elem[2][2] << ']';
		return oss.str();
	}
};

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const matrix<3, 3, _T>& mat)
{
	return mat.elem[0][0] * (mat.elem[1][1] * mat.elem[2][2] - mat.elem[1][2] * mat.elem[2][1])
			+ mat.elem[0][1] * (mat.elem[1][2] * mat.elem[2][0] - mat.elem[1][0] * mat.elem[2][2])
			+ mat.elem[0][2] * (mat.elem[1][0] * mat.elem[2][1] - mat.elem[1][1] * mat.elem[2][0]);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<3, 3, _Result> hadamard(const matrix<3, 3, _T1>& mat1, const matrix<3, 3, _T2>& mat2)
{
	_Result m_elem[] =
			{
					mat1.elem[0][0] * mat2.elem[0][0], mat1.elem[0][1] * mat2.elem[0][1], mat1.elem[0][2] * mat2.elem[0][2],
					mat1.elem[1][0] * mat2.elem[1][0], mat1.elem[1][1] * mat2.elem[1][1], mat1.elem[1][2] * mat2.elem[1][2],
					mat1.elem[2][0] * mat2.elem[2][0], mat1.elem[2][1] * mat2.elem[2][1], mat1.elem[2][2] * mat2.elem[2][2]
			};
	return tplmp::cast<matrix<3, 3, _Result> >(m_elem);
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
inline _Result tr(const matrix<3, 3, _T>& mat)
{
	return mat.elem[0][0] + mat.elem[1][1] + mat.elem[2][2];
}

template<typename _T>
struct __identity_impl<matrix<3, 3, _T> >
{
	static matrix<3, 3, _T> value()
	{
		_T m_elem[] =
				{
						1, 0, 0,
						0, 1, 0,
						0, 0, 1
				};
		return tplmp::cast<matrix<3, 3, _T> >(m_elem);
	}
};

//4x4方阵特化
template<typename _T>
class matrix<4, 4, _T>
{
public:
	__def_matrix_members__(_T, 4, 4)

	__def_matrix_op_cast_ptr__(_T)

	__def_matrix_op_cast_vec_ptr__(_T, 4)

	__def_matrix_op_row_vec__(_T, 4)

	template<typename _T2>
	inline operator matrix<4, 4, _T2>()
	{
		_T2 m_elem[] =
				{
						(_T2)elem[0][0], (_T2)elem[0][1], (_T2)elem[0][2], (_T2)elem[0][3],
						(_T2)elem[1][0], (_T2)elem[1][1], (_T2)elem[1][2], (_T2)elem[1][3],
						(_T2)elem[2][0], (_T2)elem[2][1], (_T2)elem[2][2], (_T2)elem[2][3],
						(_T2)elem[3][0], (_T2)elem[3][1], (_T2)elem[3][2], (_T2)elem[3][3]
				};
		return tplmp::cast<matrix<4, 4, _T2> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<4, 4, _Result> operator+(const matrix<4, 4, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] + mat.elem[0][0], elem[0][1] + mat.elem[0][1], elem[0][2] + mat.elem[0][2], elem[0][3] + mat.elem[0][3],
						elem[1][0] + mat.elem[1][0], elem[1][1] + mat.elem[1][1], elem[1][2] + mat.elem[1][2], elem[1][3] + mat.elem[1][3],
						elem[2][0] + mat.elem[2][0], elem[2][1] + mat.elem[2][1], elem[2][2] + mat.elem[2][2], elem[2][3] + mat.elem[2][3],
						elem[3][0] + mat.elem[3][0], elem[3][1] + mat.elem[3][1], elem[3][2] + mat.elem[3][2], elem[3][3] + mat.elem[3][3]
				};
		return tplmp::cast<matrix<4, 4, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<4, 4, _Result> operator-(const matrix<4, 4, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] - mat.elem[0][0], elem[0][1] - mat.elem[0][1], elem[0][2] - mat.elem[0][2], elem[0][3] - mat.elem[0][3],
						elem[1][0] - mat.elem[1][0], elem[1][1] - mat.elem[1][1], elem[1][2] - mat.elem[1][2], elem[1][3] - mat.elem[1][3],
						elem[2][0] - mat.elem[2][0], elem[2][1] - mat.elem[2][1], elem[2][2] - mat.elem[2][2], elem[2][3] - mat.elem[2][3],
						elem[3][0] - mat.elem[3][0], elem[3][1] - mat.elem[3][1], elem[3][2] - mat.elem[3][2], elem[3][3] - mat.elem[3][3]
				};
		return tplmp::cast<matrix<4, 4, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<4, 4, _Result> operator*(_T2 t) const
	{
		_Result m_elem[] =
				{
						elem[0][0] * t, elem[0][1] * t, elem[0][2] * t, elem[0][3] * t,
						elem[1][0] * t, elem[1][1] * t, elem[1][2] * t, elem[1][3] * t,
						elem[2][0] * t, elem[2][1] * t, elem[2][2] * t, elem[2][3] * t,
						elem[3][0] * t, elem[3][1] * t, elem[3][2] * t, elem[3][3] * t
				};
		return tplmp::cast<matrix<4, 4, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<4, 4, _Result> operator*(const matrix<4, 4, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] * mat.elem[0][0] + elem[0][1] * mat.elem[1][0] + elem[0][2] * mat.elem[2][0] + elem[0][3] * mat.elem[3][0],
						elem[0][0] * mat.elem[0][1] + elem[0][1] * mat.elem[1][1] + elem[0][2] * mat.elem[2][1] + elem[0][3] * mat.elem[3][1],
						elem[0][0] * mat.elem[0][2] + elem[0][1] * mat.elem[1][2] + elem[0][2] * mat.elem[2][2] + elem[0][3] * mat.elem[3][2],
						elem[0][0] * mat.elem[0][3] + elem[0][1] * mat.elem[1][3] + elem[0][2] * mat.elem[2][3] + elem[0][3] * mat.elem[3][3],

						elem[1][0] * mat.elem[0][0] + elem[1][1] * mat.elem[1][0] + elem[1][2] * mat.elem[2][0] + elem[1][3] * mat.elem[3][0],
						elem[1][0] * mat.elem[0][1] + elem[1][1] * mat.elem[1][1] + elem[1][2] * mat.elem[2][1] + elem[1][3] * mat.elem[3][1],
						elem[1][0] * mat.elem[0][2] + elem[1][1] * mat.elem[1][2] + elem[1][2] * mat.elem[2][2] + elem[1][3] * mat.elem[3][2],
						elem[1][0] * mat.elem[0][3] + elem[1][1] * mat.elem[1][3] + elem[1][2] * mat.elem[2][3] + elem[1][3] * mat.elem[3][3],

						elem[2][0] * mat.elem[0][0] + elem[2][1] * mat.elem[1][0] + elem[2][2] * mat.elem[2][0] + elem[2][3] * mat.elem[3][0],
						elem[2][0] * mat.elem[0][1] + elem[2][1] * mat.elem[1][1] + elem[2][2] * mat.elem[2][1] + elem[2][3] * mat.elem[3][1],
						elem[2][0] * mat.elem[0][2] + elem[2][1] * mat.elem[1][2] + elem[2][2] * mat.elem[2][2] + elem[2][3] * mat.elem[3][2],
						elem[2][0] * mat.elem[0][3] + elem[2][1] * mat.elem[1][3] + elem[2][2] * mat.elem[2][3] + elem[2][3] * mat.elem[3][3],

						elem[3][0] * mat.elem[0][0] + elem[3][1] * mat.elem[1][0] + elem[3][2] * mat.elem[2][0] + elem[3][3] * mat.elem[3][0],
						elem[3][0] * mat.elem[0][1] + elem[3][1] * mat.elem[1][1] + elem[3][2] * mat.elem[2][1] + elem[3][3] * mat.elem[3][1],
						elem[3][0] * mat.elem[0][2] + elem[3][1] * mat.elem[1][2] + elem[3][2] * mat.elem[2][2] + elem[3][3] * mat.elem[3][2],
						elem[3][0] * mat.elem[0][3] + elem[3][1] * mat.elem[1][3] + elem[3][2] * mat.elem[2][3] + elem[3][3] * mat.elem[3][3]
				};
		return tplmp::cast<matrix<4, 4, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<4, _Result> operator*(const vector<4, _T2>& vec) const
	{
		_Result v_elem[] =
				{
						elem[0][0] * vec[0] + elem[0][1] * vec[1] + elem[0][2] * vec[2] + elem[0][3] * vec[3],
						elem[1][0] * vec[0] + elem[1][1] * vec[1] + elem[1][2] * vec[2] + elem[1][3] * vec[3],
						elem[2][0] * vec[0] + elem[2][1] * vec[1] + elem[2][2] * vec[2] + elem[2][3] * vec[3],
						elem[3][0] * vec[0] + elem[3][1] * vec[1] + elem[3][2] * vec[2] + elem[3][3] * vec[3]
				};
		return tplmp::cast<vector<4, _Result> >(v_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<4, 4, _Result> operator/(_T2 t) const
	{
		_Result m_elem[] =
				{
						elem[0][0] / t, elem[0][1] / t, elem[0][2] / t, elem[0][3] / t,
						elem[1][0] / t, elem[1][1] / t, elem[1][2] / t, elem[1][3] / t,
						elem[2][0] / t, elem[2][1] / t, elem[2][2] / t, elem[2][3] / t,
						elem[3][0] / t, elem[3][1] / t, elem[3][2] / t, elem[3][3] / t
				};
		return tplmp::cast<matrix<4, 4, _Result> >(m_elem);
	}

	inline vector<4, _T> column(int column) const
	{
		_T v_elem[] =
				{
						elem[0][column], elem[1][column], elem[2][column], elem[3][column]
				};
		return tplmp::cast<vector<4, _T> >(v_elem);
	}

	inline matrix<4, 4, _T> operator~() const
	{
		_T m_elem[] =
				{
						elem[0][0], elem[1][0], elem[2][0], elem[3][0],
						elem[0][1], elem[1][1], elem[2][1], elem[3][1],
						elem[0][2], elem[1][2], elem[2][2], elem[3][2],
						elem[0][3], elem[1][3], elem[2][3], elem[3][3]
				};
		return tplmp::cast<matrix<4, 4, _T> >(m_elem);
	}

	inline static matrix<4, 4, _T> zero()
	{
		_T m_elem[] =
				{
						0, 0, 0, 0,
						0, 0, 0, 0,
						0, 0, 0, 0,
						0, 0, 0, 0
				};
		return tplmp::cast<matrix<4, 4, _T> >(m_elem);
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << ' ' << elem[0][1] << ' ' << elem[0][2] << ' ' << elem[0][3] << "]\n"
				<< '[' << elem[1][0] << ' ' << elem[1][1] << ' ' << elem[1][2] << ' ' << elem[1][3] << "]\n"
				<< '[' << elem[2][0] << ' ' << elem[2][1] << ' ' << elem[2][2] << ' ' << elem[2][3] << "]\n"
				<< '[' << elem[3][0] << ' ' << elem[3][1] << ' ' << elem[3][2] << ' ' << elem[3][3] << ']';
		return oss.str();
	}
};

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const matrix<4, 4, _T>& mat)
{
	return mat.elem[0][0] * (mat.elem[1][1] * mat.elem[2][2] * mat.elem[3][3] - mat.elem[3][1] * mat.elem[2][2] * mat.elem[1][3])
			+ mat.elem[0][1] * (mat.elem[1][2] * mat.elem[2][3] * mat.elem[3][0] - mat.elem[1][0] * mat.elem[2][3] * mat.elem[3][2])
			+ mat.elem[0][2] * (mat.elem[1][3] * mat.elem[2][0] * mat.elem[3][1] - mat.elem[1][1] * mat.elem[2][0] * mat.elem[3][3])
			+ mat.elem[0][3] * (mat.elem[1][0] * mat.elem[2][1] * mat.elem[3][2] - mat.elem[1][2] * mat.elem[2][1] * mat.elem[3][0]);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<4, 4, _Result> hadamard(const matrix<4, 4, _T1>& mat1, const matrix<4, 4, _T2>& mat2)
{
	_Result m_elem[] =
			{
					mat1.elem[0][0] * mat2.elem[0][0], mat1.elem[0][1] * mat2.elem[0][1], mat1.elem[0][2] * mat2.elem[0][2], mat1.elem[0][3] * mat2.elem[0][3],
					mat1.elem[1][0] * mat2.elem[1][0], mat1.elem[1][1] * mat2.elem[1][1], mat1.elem[1][2] * mat2.elem[1][2], mat1.elem[1][3] * mat2.elem[1][3],
					mat1.elem[2][0] * mat2.elem[2][0], mat1.elem[2][1] * mat2.elem[2][1], mat1.elem[2][2] * mat2.elem[2][2], mat1.elem[2][3] * mat2.elem[2][3],
					mat1.elem[3][0] * mat2.elem[3][0], mat1.elem[3][1] * mat2.elem[3][1], mat1.elem[3][2] * mat2.elem[3][2], mat1.elem[3][3] * mat2.elem[3][3]
			};
	return tplmp::cast<matrix<4, 4, _Result> >(m_elem);
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
inline _Result tr(const matrix<4, 4, _T>& mat)
{
	return mat.elem[0][0] + mat.elem[1][1] + mat.elem[2][2] + mat.elem[3][3];
}

template<typename _T>
struct __identity_impl<matrix<4, 4, _T> >
{
	static matrix<4, 4, _T> value()
	{
		_T m_elem[] =
				{
						1, 0, 0, 0,
						0, 1, 0, 0,
						0, 0, 1, 0,
						0, 0, 0, 1
				};
		return tplmp::cast<matrix<4, 4, _T> >(m_elem);
	}
};

template<typename _T>
using matrix2x2 = matrix<2, 2, _T>;
template<typename _T>
using matrix3x3 = matrix<3, 3, _T>;
template<typename _T>
using matrix4x4 = matrix<4, 4, _T>;

using matrix2x2f = matrix2x2<float>;
using matrix3x3f = matrix3x3<float>;
using matrix4x4f = matrix4x4<float>;

using matrix2x2d = matrix2x2<double>;
using matrix3x3d = matrix3x3<double>;
using matrix4x4d = matrix4x4<double>;
}
#endif//_MATH_MATRIX
