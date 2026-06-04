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
	inline math::vector<column, elem_type> row(size_t row)\
	{\
		return elem[row];\
	}
	__def_matrix_op_row_vec__(_T, _Column)

	/**
	 * 找列主元。
	 * 如果找到则返回行索引，未找到（全0）则返回row
	 */
#define __def_matrix_op_partial_pivot_row__(elem_type, row)\
	size_t partial_pivot_row(size_t i, size_t j, elem_type eps)\
	{\
		size_t pivot_row = row;\
		elem_type max_abs = elem_type(0);\
		for(size_t r = i; r < row; ++r)\
		{\
			elem_type abs_val = math::abs(elem[r][j]);\
			if(abs_val > max_abs)\
			{\
				max_abs = abs_val;\
				pivot_row = r;\
			}\
		}\
		if(math::is_zero(max_abs, eps))\
			return row;\
		else\
			return pivot_row;\
	}
	__def_matrix_op_partial_pivot_row__(_T, _Row)

	/**
	 * 找全主元。
	 * 如果找到则返回行索引和列索引，未找到则返回row。
	 */
#define __def_matrix_op_complete_pivot_coord__(elem_type, row, column)\
	vector2<size_t> complete_pivot_coord(size_t i, size_t j, elem_type eps)\
	{\
		vector2<size_t> pivot_coord = {row, column};\
		elem_type max_abs = elem_type(0);\
		for(size_t r = i; r < row; ++r)\
		{\
			for(size_t c = j; c < column; ++c)\
			{\
				elem_type abs_val = math::abs(elem[r][c]);\
				if(abs_val > max_abs)\
				{\
					max_abs = abs_val;\
					pivot_coord[0] = r;\
					pivot_coord[1] = c;\
				}\
			}\
		}\
		if(math::is_zero(max_abs, eps))\
			return {row, column};\
		else\
			return pivot_coord;\
	}
	__def_matrix_op_complete_pivot_coord__(_T, _Row, _Column)

	/**
	 * 交换本矩阵的行
	 */
#define __def_matrix_op_swap_row__(elem_type, row, column)\
	inline matrix<row, column, elem_type> swap_row(size_t row1, size_t row2)\
	{\
		vector<column, elem_type> temp = elem[row1];\
		elem[row1] = elem[row2];\
		elem[row2] = temp;\
		return *this;\
	}
	__def_matrix_op_swap_row__(_T, _Row, _Column)

	/**
	 * 交换本矩阵的列
	 */
#define __def_matrix_op_swap_column__(elem_type, row, column)\
	inline matrix<row, column, elem_type> swap_column(size_t column1, size_t column2)\
	{\
		for(size_t i = 0; i < row; ++i)\
		{\
			elem_type temp = elem[i][column1];\
			elem[i][column1] = elem[i][column2];\
			elem[i][column2] = temp;\
		}\
		return *this;\
	}
			__def_matrix_op_swap_column__(_T, _Row, _Column)

			/**
			 * 提取列向量
			 */
			inline math::vector<_Row,
	_T> column(size_t column)
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
					oss << '\t';
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
	return mat.template operator*<_T1, _Result>(t);
}

//基
template<size_t _Row, size_t _Column, typename _T>
struct __basis_impl<matrix<_Row, _Column, _T> >
{
	static matrix<_Row, _Column, _T> value(size_t i, size_t j)
	{
		matrix<_Row, _Column, _T> basis = matrix<_Row, _Column, _T>::zero();
		basis[i][j] = 1;
		return basis;
	}
};

/**
 * Row Echelon Form，即变换为行阶梯形
 * 列主元高斯消元法计算
 */
template<size_t _Row, size_t _Column, typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
matrix<_Row, _Column, _Result> ref(const matrix<_Row, _Column, _T>& mat, _T eps, size_t* mat_rank = nullptr, size_t* swap_count = nullptr)
{
	if(swap_count)
		*swap_count = 0;
	size_t rank = 0; //当前有主元的行数，即秩
	matrix<_Row, _Column, _Result> ref = mat;
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
 * 秩
 * 相对于ref()的优化
 */
template<size_t _Row, size_t _Column, typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
size_t rank(const matrix<_Row, _Column, _T>& mat, _T eps)
{
	matrix<_Row, _Column, _Result> ref = mat;
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
template<size_t _Order, typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
_Result det(const matrix<_Order, _Order, _T>& mat, _T eps)
{
	matrix<_Order, _Order, _T> temp = mat;
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

	__def_matrix_op_partial_pivot_row__(_T, 2)

	__def_matrix_op_complete_pivot_coord__(_T, 2, 2)

	__def_matrix_op_swap_row__(_T, 2, 2)

	//交换列
			inline matrix<2,
	2, _T>& swap_column(size_t column1, size_t column2)
	{
		_T temp = elem[0][column1];
		elem[0][column1] = elem[0][column2];
		elem[0][column2] = temp;
		temp = elem[1][column1];
		elem[1][column1] = elem[1][column2];
		elem[1][column2] = temp;
		return *this;
	}

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

	inline vector<2, _T> column(size_t column) const
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
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << "]\n"
				<< '[' << elem[1][0] << '\t' << elem[1][1] << ']';
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

	__def_matrix_op_partial_pivot_row__(_T, 3)

	__def_matrix_op_complete_pivot_coord__(_T, 3, 3)

	//交换行
	__def_matrix_op_swap_row__(_T, 3, 3)

	//交换列
			inline matrix<3,
	3, _T>& swap_column(size_t column1, size_t column2)
	{
		_T temp = elem[0][column1];
		elem[0][column1] = elem[0][column2];
		elem[0][column2] = temp;
		temp = elem[1][column1];
		elem[1][column1] = elem[1][column2];
		elem[1][column2] = temp;
		temp = elem[2][column1];
		elem[2][column1] = elem[2][column2];
		elem[2][column2] = temp;
		return *this;
	}

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

	inline vector<3, _T> column(size_t column)
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
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << '\t' << elem[0][2] << "]\n"
				<< '[' << elem[1][0] << '\t' << elem[1][1] << '\t' << elem[1][2] << "]\n"
				<< '[' << elem[2][0] << '\t' << elem[2][1] << '\t' << elem[2][2] << ']';
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

	__def_matrix_op_partial_pivot_row__(_T, 4)

	__def_matrix_op_complete_pivot_coord__(_T, 4, 4)

	//交换行
	__def_matrix_op_swap_row__(_T, 4, 4)

	//交换列
			inline matrix<4,
	4, _T>& swap_column(size_t column1, size_t column2)
	{
		_T temp = elem[0][column1];
		elem[0][column1] = elem[0][column2];
		elem[0][column2] = temp;
		temp = elem[1][column1];
		elem[1][column1] = elem[1][column2];
		elem[1][column2] = temp;
		temp = elem[2][column1];
		elem[2][column1] = elem[2][column2];
		elem[2][column2] = temp;
		temp = elem[3][column1];
		elem[3][column1] = elem[3][column2];
		elem[3][column2] = temp;
		return *this;
	}

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

	inline vector<4, _T> column(size_t column) const
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
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << '\t' << elem[0][2] << '\t' << elem[0][3] << "]\n"
				<< '[' << elem[1][0] << '\t' << elem[1][1] << '\t' << elem[1][2] << '\t' << elem[1][3] << "]\n"
				<< '[' << elem[2][0] << '\t' << elem[2][1] << '\t' << elem[2][2] << '\t' << elem[2][3] << "]\n"
				<< '[' << elem[3][0] << '\t' << elem[3][1] << '\t' << elem[3][2] << '\t' << elem[3][3] << ']';
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

using matrix2x2q = matrix2x2<__float128>;
using matrix3x3q = matrix3x3<__float128>;
using matrix4x4q = matrix4x4<__float128>;
}
#endif//_MATH_MATRIX
