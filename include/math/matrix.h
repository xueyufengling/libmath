#ifndef _MATH_MATRIX
#define _MATH_MATRIX

#include <tplmp/base.h>
#include <tplmp/array.h>

#include <string.h>
#include <sstream>

#include <math/vector.h>

namespace math
{
/**
 * 矩阵。
 * 采用行主序以贴合矩阵乘向量的数学含义
 */
template<size_t _Row, size_t _Column, typename _T>
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
		return elem;\
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
		if(math::algebra::is_zero(max_abs, eps))\
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
		if(math::algebra::is_zero(max_abs, eps))\
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
			 * 矩阵分块，分块区间为[_BeginRow, _EndRow)，[_BeginColumn, _EndColumn)
			 */
#define __def_matrix_op_partition__(elem_type, row, column)\
	template<size_t _PartitionRow, size_t _PartitionColumn>\
	inline matrix<_PartitionRow, _PartitionColumn, elem_type> partition(size_t row_begin, size_t column_begin) const\
	{\
		static_assert(_PartitionRow >= 0 && _PartitionColumn >= 0, "invalid partition size");\
		matrix<_PartitionRow, _PartitionColumn, elem_type> partition_result;\
		for(size_t i = 0; i < _PartitionRow; ++i)\
		{\
			partition_result.elem[i] = elem[row_begin + i].template slice<_PartitionColumn>(column_begin);\
		}\
		return partition_result;\
	}

			__def_matrix_op_partition__(_T, _Row, _Column)

#define __def_matrix_op_mask__(elem_type, row, column)\
	template<typename _Result = _T>\
	inline matrix<row, column, _Result> mask(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end, _Result mask_value) const\
	{\
		matrix<row, column, _Result> mask_mat;\
		for(size_t i = 0; i < row_begin; ++i)\
		{\
			mask_mat.elem[i] = vector<column, _Result>::all(mask_value);\
		}\
		for(size_t i = row_begin; i < row_end; ++i)\
		{\
			mask_mat.elem[i] = elem[i].template mask<_Result>(column_begin, column_end, mask_value);\
		}\
		for(size_t i = row_end; i < row; ++i)\
		{\
			mask_mat.elem[i] = vector<column, _Result>::all(mask_value);\
		}\
		return mask_mat;\
	}

			__def_matrix_op_mask__(_T, _Row, _Column)

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

	/**
	 * 视为_Row x _Column维向量
	 */
	operator vector<_Row * _Column, _T>&() const
	{
		return *(vector<_Row * _Column, _T>*)this;
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
	 * 取反
	 */
	inline matrix<_Row, _Column, _T> operator-() const
	{
		matrix<_Row, _Column, _T> inv;
		for(size_t i = 0; i < _Row; ++i)
			inv[i] = -elem[i];
		return inv;
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

	template<typename _T2>
	inline matrix<_Row, _Column, _T>& operator*=(const _T2& t)
	{
		for(size_t i = 0; i < _Row; ++i)
			elem[i] *= t;
		return *this;
	}

	/**
	 * 矩阵乘法
	 */
#define __def_matrix_op_matrix_mul__(elem_type, row_, column_)\
	template<size_t _Row2, typename _T2, typename _Result = decltype(tplmp::decl<elem_type>::val() * tplmp::decl<_T2>::val())>\
	inline matrix<row_, _Row2, _Result> operator*(const matrix<column_, _Row2, _T2>& mat) const\
	{\
		matrix<row_, _Row2, _Result> mul_result;\
		for(size_t i = 0; i < row_; ++i)\
		{\
			for(size_t j = 0; j < _Row2; ++j)\
				mul_result[i][j] = this->row(i) * mat.column(j);\
		}\
		return mul_result;\
	}
	__def_matrix_op_matrix_mul__(_T, _Row, _Column)

	/**
	 * 向量乘法
	 */
#define __def_matrix_op_vector_mul__(elem_type, row_, column_)\
	template<typename _T2, typename _Result = decltype(tplmp::decl<elem_type>::val() * tplmp::decl<_T2>::val())>\
	inline math::vector<row_, _Result> operator*(const math::vector<column_, _T2>& vec) const\
	{\
		vector<row_, _Result> mul_result;\
		for(size_t i = 0; i < row_; ++i)\
		mul_result[i] = row(i) * vec;\
		return mul_result;\
	}
			__def_matrix_op_vector_mul__(_T, _Row, _Column)

			template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator/(_T2 t) const
	{
		matrix<_Row, _Column, _Result> div_result;
		for(size_t i = 0; i < _Row; ++i)
		div_result[i] = elem[i] / t;
		return div_result;
	}

	/**
	 * 按位相除，同哈达玛积
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator/(const matrix<_Row, _Column, _T2>& mat) const
	{
		matrix<_Row, _Column, _Result> div_result;
		for(size_t i = 0; i < _Row; ++i)
		div_result[i] = elem[i] / mat.elem[i];
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

	inline static matrix<_Row, _Column, _T> all(_T value)
	{
		matrix<_Row, _Column, _T> mat;
		mat.elem[0] = vector<_Column, _T>::all(value);
		for(size_t i = 1; i < _Row; ++i)
		mat.elem[i] = mat.elem[0];
		return mat;
	}

#define __def_matrix_zero_one__(elem_type, row_, column_)\
	inline static matrix<row_, column_, elem_type> zero()\
	{\
		return all(0);\
	}\
	inline static matrix<row_, column_, elem_type> one()\
	{\
		return all(1);\
	}
	__def_matrix_zero_one__(_T, _Row, _Column)

#define __def_matrix_op_sum__(elem_type, row_, column_)\
	template<typename _Result = decltype(tplmp::decl<elem_type>::val() + tplmp::decl<elem_type>::val())>\
	inline _Result sum_row(size_t i) const\
	{\
		_T sum(0);\
		for(size_t j = 0; j < column_; ++j)\
		{\
			sum += elem[i][j];\
		}\
		return sum;\
	}\
	template<typename _Result = decltype(tplmp::decl<elem_type>::val() + tplmp::decl<elem_type>::val())>\
	inline _Result sum_column(size_t j) const\
	{\
		_T sum(0);\
		for(size_t i = 0; i < row_; ++i)\
		{\
			sum += elem[i][j];\
		}\
		return sum;\
	}
	__def_matrix_op_sum__(_T, _Row, _Column)

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
	matrix<_Row, _Column, _Result> mul_result;
	for(size_t i = 0; i < _Row; ++i)
		mul_result[i] = t * mat.elem[i];
	return mul_result;
}

template<typename _T1, size_t _Row, size_t _Column, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<_Row, _Column, _Result> operator/(_T1 t, const matrix<_Row, _Column, _T2>& mat)
{
	matrix<_Row, _Column, _Result> div_result;
	for(size_t i = 0; i < _Row; ++i)
		div_result[i] = t / mat.elem[i];
	return div_result;
}

template<size_t _Row, size_t _Column1, size_t _Column2, typename _T>
inline matrix<_Row, _Column1 + _Column2, _T> cat_column(const matrix<_Row, _Column1, _T>& mat1, const matrix<_Row, _Column2, _T>& mat2)
{
	matrix<_Row, _Column1 + _Column2, _T> mat;
	for(size_t i = 0; i < _Row; ++i)
		mat.elem[i] = cat(mat1.elem[i], mat2.elem[i]);
	return mat;
}

/**
 * @brief 4个分块矩阵拼接
 */
template<size_t _Row1, size_t _Column1, size_t _Row2, size_t _Column2, typename _T>
inline matrix<_Row1 + _Row2, _Column1 + _Column2, _T> cat(const matrix<_Row1, _Column1, _T>& left_up, const matrix<_Row1, _Column2, _T>& right_up, const matrix<_Row2, _Column1, _T>& left_down, const matrix<_Row2, _Column2, _T>& right_down)
{
	constexpr size_t cat_row = _Row1 + _Row2;
	constexpr size_t cat_column = _Column1 + _Column2;
	matrix<cat_row, cat_column, _T> mat;
	for(size_t i = 0; i < _Row1; ++i)
		mat.elem[i] = cat(left_up.elem[i], right_up.elem[i]);
	for(size_t i = _Row1; i < cat_row; ++i)
		mat.elem[i] = cat(left_down.elem[i], right_down.elem[i]);
	return mat;
}

namespace algebra
{
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
}

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
namespace algebra
{
//加法恒等元：零矩阵
template<size_t _Row, size_t _Column, typename _T>
struct __identity_impl<matrix<_Row, _Column, _T>, algebra::add>
{
	inline static matrix<_Row, _Column, _T> value()
	{
		return matrix<_Row, _Column, _T>::zero();
	}
};
//乘法恒等元：恒等矩阵
template<size_t _Order, typename _T>
struct __identity_impl<matrix<_Order, _Order, _T>, algebra::mul>
{
	inline static matrix<_Order, _Order, _T> value()
	{
		matrix<_Order, _Order, _T> identity = matrix<_Order, _Order, _T>::zero();
		for(size_t i = 0; i < _Order; ++i)
			identity[i][i] = 1;
		return identity;
	}
};
//乘法零元：零矩阵
template<size_t _Row, size_t _Column, typename _T>
struct __zero_impl<matrix<_Row, _Column, _T>, algebra::mul>
{
	inline static matrix<_Row, _Column, _T> value()
	{
		return matrix<_Row, _Column, _T>::zero();
	}
};
}

// 1x1 方阵特化
template<typename _T>
class matrix<1, 1, _T>
{
public:
	__def_matrix_members__(_T, 1, 1)

	__def_matrix_op_cast_ptr__(_T)

	__def_matrix_op_cast_vec_ptr__(_T, 1)

	__def_matrix_op_row_vec__(_T, 1)

	__def_matrix_op_partial_pivot_row__(_T, 1)

	__def_matrix_op_complete_pivot_coord__(_T, 1, 1)

	//交换行
	__def_matrix_op_swap_row__(_T, 1, 1)

	__def_matrix_op_sum__(_T, 1, 1)

	__def_matrix_op_partition__(_T, 1, 1)

	__def_matrix_op_mask__(_T, 1, 1)

	inline matrix<1,
	1, _T>& swap_column(size_t column1, size_t column2)
	{
		return *this;
	}

	template<typename _T2>
	inline operator matrix<1, 1, _T2>()
	{
		_T2 m_elem[] =
				{
						(_T2)elem[0][0]
				};
		return tplmp::cast<matrix<1, 1, _T2> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<1, 1, _Result> operator+(const matrix<1, 1, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] + mat.elem[0][0]
				};
		return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<1, 1, _Result> operator-(const matrix<1, 1, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] - mat.elem[0][0]
				};
		return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
	}

	inline matrix<1, 1, _T> operator-() const
	{
		_T m_elem[] =
				{
						-elem[0][0]
				};
		return tplmp::cast<matrix<1, 1, _T> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<1, 1, _Result> operator*(_T2 t) const
	{
		_Result m_elem[] =
				{
						elem[0][0] * t
				};
		return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<1, 1, _T>& operator*=(const _T2& t)
	{
		elem[0] *= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<1, 1, _Result> operator*(const matrix<1, 1, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] * mat.elem[0][0]
				};
		return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline vector<1, _Result> operator*(const vector<1, _T2>& vec) const
	{
		_Result v_elem[] =
				{
						elem[0][0] * vec[0]
				};
		return tplmp::cast<vector<1, _Result> >(v_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<1, 1, _Result> operator/(_T2 t) const
	{
		_Result m_elem[] =
				{
						elem[0][0] / t
				};
		return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<1, 1, _T>& operator/=(const _T2& t)
	{
		elem[0] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<1, 1, _Result> operator/(const matrix<1, 1, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] / mat.elem[0][0]
				};
		return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
	}

	inline vector<1, _T> column(size_t column) const
	{
		_T v_elem[] =
				{
						elem[0][column]
				};
		return tplmp::cast<vector<1, _T> >(v_elem);
	}

	inline matrix<1, 1, _T> operator~() const
	{
		_T m_elem[] =
				{
						elem[0][0]
				};
		return tplmp::cast<matrix<1, 1, _T> >(m_elem);
	}

	inline static matrix<1, 1, _T> all(_T value)
	{
		_T m_elem[] = {value};
		return tplmp::cast<matrix<1, 1, _T> >(m_elem);
	}
	__def_matrix_zero_one__(_T, 1, 1)

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << ']';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<1, 1, _Result> operator*(_T1 t, const matrix<1, 1, _T2>& mat)
{
	_Result m_elem[] =
			{
					t * mat.elem[0][0]
			};
	return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<1, 1, _Result> operator/(_T1 t, const matrix<1, 1, _T2>& mat)
{
	_Result m_elem[] =
			{
					t / mat.elem[0][0]
			};
	return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const matrix<1, 1, _T>& mat)
{
	return mat.elem[0][0];
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<1, 1, _Result> hadamard(const matrix<1, 1, _T1>& mat1, const matrix<1, 1, _T2>& mat2)
{
	_Result m_elem[] =
			{
					mat1.elem[0][0] * mat2.elem[0][0]
			};
	return tplmp::cast<matrix<1, 1, _Result> >(m_elem);
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
inline _Result tr(const matrix<1, 1, _T>& mat)
{
	return mat.elem[0][0];
}

namespace algebra
{
template<typename _T>
struct __identity_impl<matrix<1, 1, _T>, algebra::mul>
{
	static matrix<1, 1, _T> value()
	{
		_T m_elem[] = {1};
		return tplmp::cast<matrix<1, 1, _T> >(m_elem);
	}
};
}

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

	__def_matrix_op_sum__(_T, 2, 2)

	__def_matrix_op_partition__(_T, 2, 2)

	__def_matrix_op_mask__(_T, 2, 2)

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

	inline matrix<2, 2, _T> operator-() const
	{
		_T m_elem[] =
				{
						-elem[0][0], -elem[0][1],
						-elem[1][0], -elem[1][1]
				};
		return tplmp::cast<matrix<2, 2, _T> >(m_elem);
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

	template<typename _T2>
	inline matrix<2, 2, _T>& operator*=(const _T2& t)
	{
		elem[0] *= t;
		elem[1] *= t;
		return *this;
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

	template<typename _T2>
	inline matrix<2, 2, _T>& operator/=(const _T2& t)
	{
		elem[0] /= t;
		elem[1] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<2, 2, _Result> operator/(const matrix<2, 2, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] / mat.elem[0][0], elem[0][1] / mat.elem[0][1],
						elem[1][0] / mat.elem[1][0], elem[1][1] / mat.elem[1][1]
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

	inline static matrix<2, 2, _T> all(_T value)
	{
		_T m_elem[] = {
				value, value,
				value, value
		};
		return tplmp::cast<matrix<2, 2, _T> >(m_elem);
	}
	__def_matrix_zero_one__(_T, 2, 2)

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << "]\n"
				<< '[' << elem[1][0] << '\t' << elem[1][1] << ']';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<2, 2, _Result> operator*(_T1 t, const matrix<2, 2, _T2>& mat)
{
	_Result m_elem[] = {
			t * mat.elem[0][0], t * mat.elem[0][1],
			t * mat.elem[1][0], t * mat.elem[1][1]
	};
	return tplmp::cast<matrix<2, 2, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<2, 2, _Result> operator/(_T1 t, const matrix<2, 2, _T2>& mat)
{
	_Result m_elem[] = {
			t / mat.elem[0][0], t / mat.elem[0][1],
			t / mat.elem[1][0], t / mat.elem[1][1]
	};
	return tplmp::cast<matrix<2, 2, _Result> >(m_elem);
}

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

namespace algebra
{
template<typename _T>
struct __identity_impl<matrix<2, 2, _T>, algebra::mul>
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
}

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

	__def_matrix_op_sum__(_T, 3, 3)

	__def_matrix_op_partition__(_T, 3, 3)

	__def_matrix_op_mask__(_T, 3, 3)

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

	inline matrix<3, 3, _T> operator-() const
	{
		_T m_elem[] =
				{
						-elem[0][0], -elem[0][1], -elem[0][2],
						-elem[1][0], -elem[1][1], -elem[1][2],
						-elem[2][0], -elem[2][1], -elem[2][2]
				};
		return tplmp::cast<matrix<3, 3, _T> >(m_elem);
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

	template<typename _T2>
	inline matrix<3, 3, _T>& operator*=(const _T2& t)
	{
		elem[0] *= t;
		elem[1] *= t;
		elem[2] *= t;
		return *this;
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

	template<typename _T2>
	inline matrix<3, 3, _T>& operator/=(const _T2& t)
	{
		elem[0] /= t;
		elem[1] /= t;
		elem[2] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<3, 3, _Result> operator/(const matrix<3, 3, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] / mat.elem[0][0], elem[0][1] / mat.elem[0][1], elem[0][2] / mat.elem[0][2],
						elem[1][0] / mat.elem[1][0], elem[1][1] / mat.elem[1][1], elem[1][2] / mat.elem[1][2],
						elem[2][0] / mat.elem[2][0], elem[2][1] / mat.elem[2][1], elem[2][2] / mat.elem[2][2]
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

	inline static matrix<3, 3, _T> all(_T value)
	{
		_T m_elem[] = {
				value, value, value,
				value, value, value,
				value, value, value
		};
		return tplmp::cast<matrix<3, 3, _T> >(m_elem);
	}
	__def_matrix_zero_one__(_T, 3, 3)

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << '\t' << elem[0][2] << "]\n"
				<< '[' << elem[1][0] << '\t' << elem[1][1] << '\t' << elem[1][2] << "]\n"
				<< '[' << elem[2][0] << '\t' << elem[2][1] << '\t' << elem[2][2] << ']';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<3, 3, _Result> operator*(_T1 t, const matrix<3, 3, _T2>& mat)
{
	_Result m_elem[] = {
			t * mat.elem[0][0], t * mat.elem[0][1], t * mat.elem[0][2],
			t * mat.elem[1][0], t * mat.elem[1][1], t * mat.elem[1][2],
			t * mat.elem[2][0], t * mat.elem[2][1], t * mat.elem[2][2]
	};
	return tplmp::cast<matrix<3, 3, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<3, 3, _Result> operator/(_T1 t, const matrix<3, 3, _T2>& mat)
{
	_Result m_elem[] = {
			t / mat.elem[0][0], t / mat.elem[0][1], t / mat.elem[0][2],
			t / mat.elem[1][0], t / mat.elem[1][1], t / mat.elem[1][2],
			t / mat.elem[2][0], t / mat.elem[2][1], t / mat.elem[2][2]
	};
	return tplmp::cast<matrix<3, 3, _Result> >(m_elem);
}

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

namespace algebra
{
template<typename _T>
struct __identity_impl<matrix<3, 3, _T>, algebra::mul>
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
}

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

	__def_matrix_op_sum__(_T, 4, 4)

	__def_matrix_op_partition__(_T, 4, 4)

	__def_matrix_op_mask__(_T, 4, 4)

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

	inline matrix<4, 4, _T> operator-() const
	{
		_T m_elem[] =
				{
						-elem[0][0], -elem[0][1], -elem[0][2], -elem[0][3],
						-elem[1][0], -elem[1][1], -elem[1][2], -elem[1][3],
						-elem[2][0], -elem[2][1], -elem[2][2], -elem[2][3],
						-elem[3][0], -elem[3][1], -elem[3][2], -elem[3][3]
				};
		return tplmp::cast<matrix<4, 4, _T> >(m_elem);
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

	template<typename _T2>
	inline matrix<4, 4, _T>& operator*=(const _T2& t)
	{
		elem[0] *= t;
		elem[1] *= t;
		elem[2] *= t;
		elem[3] *= t;
		return *this;
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

	template<typename _T2>
	inline matrix<4, 4, _T>& operator/=(const _T2& t)
	{
		elem[0] /= t;
		elem[1] /= t;
		elem[2] /= t;
		elem[3] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<4, 4, _Result> operator/(const matrix<4, 4, _T2>& mat) const
	{
		_Result m_elem[] =
				{
						elem[0][0] / mat.elem[0][0], elem[0][1] / mat.elem[0][1], elem[0][2] / mat.elem[0][2], elem[0][3] / mat.elem[0][3],
						elem[1][0] / mat.elem[1][0], elem[1][1] / mat.elem[1][1], elem[1][2] / mat.elem[1][2], elem[1][3] / mat.elem[1][3],
						elem[2][0] / mat.elem[2][0], elem[2][1] / mat.elem[2][1], elem[2][2] / mat.elem[2][2], elem[2][3] / mat.elem[2][3],
						elem[3][0] / mat.elem[3][0], elem[3][1] / mat.elem[3][1], elem[3][2] / mat.elem[3][2], elem[3][3] / mat.elem[3][3]
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

	inline static matrix<4, 4, _T> all(_T value)
	{
		_T m_elem[] = {
				value, value, value, value,
				value, value, value, value,
				value, value, value, value,
				value, value, value, value
		};
		return tplmp::cast<matrix<4, 4, _T> >(m_elem);
	}
	__def_matrix_zero_one__(_T, 4, 4)

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

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<4, 4, _Result> operator*(_T1 t, const matrix<4, 4, _T2>& mat)
{
	_Result m_elem[] = {
			t * mat.elem[0][0], t * mat.elem[0][1], t * mat.elem[0][2], t * mat.elem[0][3],
			t * mat.elem[1][0], t * mat.elem[1][1], t * mat.elem[1][2], t * mat.elem[1][3],
			t * mat.elem[2][0], t * mat.elem[2][1], t * mat.elem[2][2], t * mat.elem[2][3],
			t * mat.elem[3][0], t * mat.elem[3][1], t * mat.elem[3][2], t * mat.elem[3][3]
	};
	return tplmp::cast<matrix<4, 4, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<4, 4, _Result> operator/(_T1 t, const matrix<4, 4, _T2>& mat)
{
	_Result m_elem[] = {
			t / mat.elem[0][0], t / mat.elem[0][1], t / mat.elem[0][2], t / mat.elem[0][3],
			t / mat.elem[1][0], t / mat.elem[1][1], t / mat.elem[1][2], t / mat.elem[1][3],
			t / mat.elem[2][0], t / mat.elem[2][1], t / mat.elem[2][2], t / mat.elem[2][3],
			t / mat.elem[3][0], t / mat.elem[3][1], t / mat.elem[3][2], t / mat.elem[3][3]
	};
	return tplmp::cast<matrix<4, 4, _Result> >(m_elem);
}

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

namespace algebra
{
template<typename _T>
struct __identity_impl<matrix<4, 4, _T>, algebra::mul>
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
}

//6x6矩阵特化
//用于分析力学相空间计算
template<typename _T>
class matrix<6, 6, _T>
{
public:
	__def_matrix_members__(_T, 6, 6)

	__def_matrix_op_cast_ptr__(_T)

	__def_matrix_op_cast_vec_ptr__(_T, 36)

	__def_matrix_op_row_vec__(_T, 6)

	__def_matrix_op_partial_pivot_row__(_T, 6)

	__def_matrix_op_complete_pivot_coord__(_T, 6, 6)

	//交换行
	__def_matrix_op_swap_row__(_T, 6, 6)

	__def_matrix_op_sum__(_T, 6, 6)

	__def_matrix_op_partition__(_T, 6, 6)

	__def_matrix_op_mask__(_T, 6, 6)

	inline matrix<6,
	6, _T>& swap_column(size_t column1, size_t column2)
	{
		return *this;
	}

	template<typename _T2>
	inline operator matrix<6, 6, _T2>()
	{
		_T2 m_elem[] = {
				(_T2)elem[0][0], (_T2)elem[0][1], (_T2)elem[0][2], (_T2)elem[0][3], (_T2)elem[0][4], (_T2)elem[0][5],
				(_T2)elem[1][0], (_T2)elem[1][1], (_T2)elem[1][2], (_T2)elem[1][3], (_T2)elem[1][4], (_T2)elem[1][5],
				(_T2)elem[2][0], (_T2)elem[2][1], (_T2)elem[2][2], (_T2)elem[2][3], (_T2)elem[2][4], (_T2)elem[2][5],
				(_T2)elem[3][0], (_T2)elem[3][1], (_T2)elem[3][2], (_T2)elem[3][3], (_T2)elem[3][4], (_T2)elem[3][5],
				(_T2)elem[4][0], (_T2)elem[4][1], (_T2)elem[4][2], (_T2)elem[4][3], (_T2)elem[4][4], (_T2)elem[4][5],
				(_T2)elem[5][0], (_T2)elem[5][1], (_T2)elem[5][2], (_T2)elem[5][3], (_T2)elem[5][4], (_T2)elem[5][5]
		};
		return tplmp::cast<matrix<6, 6, _T2> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<6, 6, _Result> operator+(const matrix<6, 6, _T2>& mat) const
	{
		_Result m_elem[] = {
				elem[0][0] + mat.elem[0][0], elem[0][1] + mat.elem[0][1], elem[0][2] + mat.elem[0][2], elem[0][3] + mat.elem[0][3], elem[0][4] + mat.elem[0][4], elem[0][5] + mat.elem[0][5],
				elem[1][0] + mat.elem[1][0], elem[1][1] + mat.elem[1][1], elem[1][2] + mat.elem[1][2], elem[1][3] + mat.elem[1][3], elem[1][4] + mat.elem[1][4], elem[1][5] + mat.elem[1][5],
				elem[2][0] + mat.elem[2][0], elem[2][1] + mat.elem[2][1], elem[2][2] + mat.elem[2][2], elem[2][3] + mat.elem[2][3], elem[2][4] + mat.elem[2][4], elem[2][5] + mat.elem[2][5],
				elem[3][0] + mat.elem[3][0], elem[3][1] + mat.elem[3][1], elem[3][2] + mat.elem[3][2], elem[3][3] + mat.elem[3][3], elem[3][4] + mat.elem[3][4], elem[3][5] + mat.elem[3][5],
				elem[4][0] + mat.elem[4][0], elem[4][1] + mat.elem[4][1], elem[4][2] + mat.elem[4][2], elem[4][3] + mat.elem[4][3], elem[4][4] + mat.elem[4][4], elem[4][5] + mat.elem[4][5],
				elem[5][0] + mat.elem[5][0], elem[5][1] + mat.elem[5][1], elem[5][2] + mat.elem[5][2], elem[5][3] + mat.elem[5][3], elem[5][4] + mat.elem[5][4], elem[5][5] + mat.elem[5][5]
		};
		return tplmp::cast<matrix<6, 6, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<6, 6, _Result> operator-(const matrix<6, 6, _T2>& mat) const
	{
		_Result m_elem[] = {
				elem[0][0] - mat.elem[0][0], elem[0][1] - mat.elem[0][1], elem[0][2] - mat.elem[0][2], elem[0][3] - mat.elem[0][3], elem[0][4] - mat.elem[0][4], elem[0][5] - mat.elem[0][5],
				elem[1][0] - mat.elem[1][0], elem[1][1] - mat.elem[1][1], elem[1][2] - mat.elem[1][2], elem[1][3] - mat.elem[1][3], elem[1][4] - mat.elem[1][4], elem[1][5] - mat.elem[1][5],
				elem[2][0] - mat.elem[2][0], elem[2][1] - mat.elem[2][1], elem[2][2] - mat.elem[2][2], elem[2][3] - mat.elem[2][3], elem[2][4] - mat.elem[2][4], elem[2][5] - mat.elem[2][5],
				elem[3][0] - mat.elem[3][0], elem[3][1] - mat.elem[3][1], elem[3][2] - mat.elem[3][2], elem[3][3] - mat.elem[3][3], elem[3][4] - mat.elem[3][4], elem[3][5] - mat.elem[3][5],
				elem[4][0] - mat.elem[4][0], elem[4][1] - mat.elem[4][1], elem[4][2] - mat.elem[4][2], elem[4][3] - mat.elem[4][3], elem[4][4] - mat.elem[4][4], elem[4][5] - mat.elem[4][5],
				elem[5][0] - mat.elem[5][0], elem[5][1] - mat.elem[5][1], elem[5][2] - mat.elem[5][2], elem[5][3] - mat.elem[5][3], elem[5][4] - mat.elem[5][4], elem[5][5] - mat.elem[5][5]
		};
		return tplmp::cast<matrix<6, 6, _Result> >(m_elem);
	}

	inline matrix<6, 6, _T> operator-() const
	{
		_T m_elem[] = {
				-elem[0][0], -elem[0][1], -elem[0][2], -elem[0][3], -elem[0][4], -elem[0][5],
				-elem[1][0], -elem[1][1], -elem[1][2], -elem[1][3], -elem[1][4], -elem[1][5],
				-elem[2][0], -elem[2][1], -elem[2][2], -elem[2][3], -elem[2][4], -elem[2][5],
				-elem[3][0], -elem[3][1], -elem[3][2], -elem[3][3], -elem[3][4], -elem[3][5],
				-elem[4][0], -elem[4][1], -elem[4][2], -elem[4][3], -elem[4][4], -elem[4][5],
				-elem[5][0], -elem[5][1], -elem[5][2], -elem[5][3], -elem[5][4], -elem[5][5]
		};
		return tplmp::cast<matrix<6, 6, _T> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<6, 6, _Result> operator*(_T2 t) const
	{
		_Result m_elem[] = {
				elem[0][0] * t, elem[0][1] * t, elem[0][2] * t, elem[0][3] * t, elem[0][4] * t, elem[0][5] * t,
				elem[1][0] * t, elem[1][1] * t, elem[1][2] * t, elem[1][3] * t, elem[1][4] * t, elem[1][5] * t,
				elem[2][0] * t, elem[2][1] * t, elem[2][2] * t, elem[2][3] * t, elem[2][4] * t, elem[2][5] * t,
				elem[3][0] * t, elem[3][1] * t, elem[3][2] * t, elem[3][3] * t, elem[3][4] * t, elem[3][5] * t,
				elem[4][0] * t, elem[4][1] * t, elem[4][2] * t, elem[4][3] * t, elem[4][4] * t, elem[4][5] * t,
				elem[5][0] * t, elem[5][1] * t, elem[5][2] * t, elem[5][3] * t, elem[5][4] * t, elem[5][5] * t
		};
		return tplmp::cast<matrix<6, 6, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<6, 6, _T>& operator*=(const _T2& t)
	{
		elem[0] *= t;
		elem[1] *= t;
		elem[2] *= t;
		elem[3] *= t;
		elem[4] *= t;
		elem[5] *= t;
		return *this;
	}

	__def_matrix_op_matrix_mul__(_T, 6, 6)

	__def_matrix_op_vector_mul__(_T, 6, 6)

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<6, 6, _Result> operator/(_T2 t) const
	{
		_Result m_elem[] =
		{
			elem[0][0] / t, elem[0][1] / t, elem[0][2] / t, elem[0][3] / t, elem[0][4] / t, elem[0][5] / t,
			elem[1][0] / t, elem[1][1] / t, elem[1][2] / t, elem[1][3] / t, elem[1][4] / t, elem[1][5] / t,
			elem[2][0] / t, elem[2][1] / t, elem[2][2] / t, elem[2][3] / t, elem[2][4] / t, elem[2][5] / t,
			elem[3][0] / t, elem[3][1] / t, elem[3][2] / t, elem[3][3] / t, elem[3][4] / t, elem[3][5] / t,
			elem[4][0] / t, elem[4][1] / t, elem[4][2] / t, elem[4][3] / t, elem[4][4] / t, elem[4][5] / t,
			elem[5][0] / t, elem[5][1] / t, elem[5][2] / t, elem[5][3] / t, elem[5][4] / t, elem[5][5] / t
		};
		return tplmp::cast<matrix<6, 6, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<6, 6, _T>& operator/=(const _T2& t)
	{
		elem[0] /= t;
		elem[1] /= t;
		elem[2] /= t;
		elem[3] /= t;
		elem[4] /= t;
		elem[5] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<6, 6, _Result> operator/(const matrix<6, 6, _T2>& mat) const
	{
		_Result m_elem[] =
		{
			elem[0][0] / mat.elem[0][0], elem[0][1] / mat.elem[0][1], elem[0][2] / mat.elem[0][2], elem[0][3] / mat.elem[0][3], elem[0][4] / mat.elem[0][4], elem[0][5] / mat.elem[0][5],
			elem[1][0] / mat.elem[1][0], elem[1][1] / mat.elem[1][1], elem[1][2] / mat.elem[1][2], elem[1][3] / mat.elem[1][3], elem[1][4] / mat.elem[1][4], elem[1][5] / mat.elem[1][5],
			elem[2][0] / mat.elem[2][0], elem[2][1] / mat.elem[2][1], elem[2][2] / mat.elem[2][2], elem[2][3] / mat.elem[2][3], elem[2][4] / mat.elem[2][4], elem[2][5] / mat.elem[2][5],
			elem[3][0] / mat.elem[3][0], elem[3][1] / mat.elem[3][1], elem[3][2] / mat.elem[3][2], elem[3][3] / mat.elem[3][3], elem[3][4] / mat.elem[3][4], elem[3][5] / mat.elem[3][5],
			elem[4][0] / mat.elem[4][0], elem[4][1] / mat.elem[4][1], elem[4][2] / mat.elem[4][2], elem[4][3] / mat.elem[4][3], elem[4][4] / mat.elem[4][4], elem[4][5] / mat.elem[4][5],
			elem[5][0] / mat.elem[5][0], elem[5][1] / mat.elem[5][1], elem[5][2] / mat.elem[5][2], elem[5][3] / mat.elem[5][3], elem[5][4] / mat.elem[5][4], elem[5][5] / mat.elem[5][5]
		};
		return tplmp::cast<matrix<6, 6, _Result> >(m_elem);
	}

	inline vector<6, _T> column(size_t column) const
	{
		_T v_elem[] =
		{
			elem[0][column],
			elem[1][column],
			elem[2][column],
			elem[3][column],
			elem[4][column],
			elem[5][column]
		};
		return tplmp::cast<vector<6, _T> >(v_elem);
	}

	inline matrix<6, 6, _T> operator~() const
	{
		_T m_elem[] =
		{
			elem[0][0], elem[1][0], elem[2][0], elem[3][0], elem[4][0], elem[5][0],
			elem[0][1], elem[1][1], elem[2][1], elem[3][1], elem[4][1], elem[5][1],
			elem[0][2], elem[1][2], elem[2][2], elem[3][2], elem[4][2], elem[5][2],
			elem[0][3], elem[1][3], elem[2][3], elem[3][3], elem[4][3], elem[5][3],
			elem[0][4], elem[1][4], elem[2][4], elem[3][4], elem[4][4], elem[5][4],
			elem[0][5], elem[1][5], elem[2][5], elem[3][5], elem[4][5], elem[5][5]
		};
		return tplmp::cast<matrix<6, 6, _T> >(m_elem);
	}

	inline static matrix<6, 6, _T> all(_T value)
	{
		_T m_elem[] =
		{
			value, value, value, value, value, value,
			value, value, value, value, value, value,
			value, value, value, value, value, value,
			value, value, value, value, value, value,
			value, value, value, value, value, value,
			value, value, value, value, value, value
		};
		return tplmp::cast<matrix<6, 6, _T> >(m_elem);
	}
	__def_matrix_zero_one__(_T, 6, 6)

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << '\t' << elem[0][2] << '\t' << elem[0][3] << '\t' << elem[0][4] << '\t' << elem[0][5] << "]\n"
		<< '[' << elem[1][0] << '\t' << elem[1][1] << '\t' << elem[1][2] << '\t' << elem[1][3] << '\t' << elem[1][4] << '\t' << elem[1][5] << "]\n"
		<< '[' << elem[2][0] << '\t' << elem[2][1] << '\t' << elem[2][2] << '\t' << elem[2][3] << '\t' << elem[2][4] << '\t' << elem[2][5] << "]\n"
		<< '[' << elem[3][0] << '\t' << elem[3][1] << '\t' << elem[3][2] << '\t' << elem[3][3] << '\t' << elem[3][4] << '\t' << elem[3][5] << "]\n"
		<< '[' << elem[4][0] << '\t' << elem[4][1] << '\t' << elem[4][2] << '\t' << elem[4][3] << '\t' << elem[4][4] << '\t' << elem[4][5] << "]\n"
		<< '[' << elem[5][0] << '\t' << elem[5][1] << '\t' << elem[5][2] << '\t' << elem[5][3] << '\t' << elem[5][4] << '\t' << elem[5][5] << ']';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<6, 6, _Result> operator*(_T1 t, const matrix<6, 6, _T2>& mat)
{
	_Result m_elem[] = {
			t * mat.elem[0][0], t * mat.elem[0][1], t * mat.elem[0][2], t * mat.elem[0][3], t * mat.elem[0][4], t * mat.elem[0][5],
			t * mat.elem[1][0], t * mat.elem[1][1], t * mat.elem[1][2], t * mat.elem[1][3], t * mat.elem[1][4], t * mat.elem[1][5],
			t * mat.elem[2][0], t * mat.elem[2][1], t * mat.elem[2][2], t * mat.elem[2][3], t * mat.elem[2][4], t * mat.elem[2][5],
			t * mat.elem[3][0], t * mat.elem[3][1], t * mat.elem[3][2], t * mat.elem[3][3], t * mat.elem[3][4], t * mat.elem[3][5],
			t * mat.elem[4][0], t * mat.elem[4][1], t * mat.elem[4][2], t * mat.elem[4][3], t * mat.elem[4][4], t * mat.elem[4][5],
			t * mat.elem[5][0], t * mat.elem[5][1], t * mat.elem[5][2], t * mat.elem[5][3], t * mat.elem[5][4], t * mat.elem[5][5]
	};
	return tplmp::cast<matrix<6, 6, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<6, 6, _Result> operator/(_T1 t, const matrix<6, 6, _T2>& mat)
{
	_Result m_elem[] = {
			t / mat.elem[0][0], t / mat.elem[0][1], t / mat.elem[0][2], t / mat.elem[0][3], t / mat.elem[0][4], t / mat.elem[0][5],
			t / mat.elem[1][0], t / mat.elem[1][1], t / mat.elem[1][2], t / mat.elem[1][3], t / mat.elem[1][4], t / mat.elem[1][5],
			t / mat.elem[2][0], t / mat.elem[2][1], t / mat.elem[2][2], t / mat.elem[2][3], t / mat.elem[2][4], t / mat.elem[2][5],
			t / mat.elem[3][0], t / mat.elem[3][1], t / mat.elem[3][2], t / mat.elem[3][3], t / mat.elem[3][4], t / mat.elem[3][5],
			t / mat.elem[4][0], t / mat.elem[4][1], t / mat.elem[4][2], t / mat.elem[4][3], t / mat.elem[4][4], t / mat.elem[4][5],
			t / mat.elem[5][0], t / mat.elem[5][1], t / mat.elem[5][2], t / mat.elem[5][3], t / mat.elem[5][4], t / mat.elem[5][5]
	};
	return tplmp::cast<matrix<6, 6, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<6, 6, _Result> hadamard(const matrix<6, 6, _T1>& mat1, const matrix<6, 6, _T2>& mat2)
{
	_Result m_elem[] = {
			mat1.elem[0][0] * mat2.elem[0][0], mat1.elem[0][1] * mat2.elem[0][1], mat1.elem[0][2] * mat2.elem[0][2], mat1.elem[0][3] * mat2.elem[0][3], mat1.elem[0][4] * mat2.elem[0][4], mat1.elem[0][5] * mat2.elem[0][5],
			mat1.elem[1][0] * mat2.elem[1][0], mat1.elem[1][1] * mat2.elem[1][1], mat1.elem[1][2] * mat2.elem[1][2], mat1.elem[1][3] * mat2.elem[1][3], mat1.elem[1][4] * mat2.elem[1][4], mat1.elem[1][5] * mat2.elem[1][5],
			mat1.elem[2][0] * mat2.elem[2][0], mat1.elem[2][1] * mat2.elem[2][1], mat1.elem[2][2] * mat2.elem[2][2], mat1.elem[2][3] * mat2.elem[2][3], mat1.elem[2][4] * mat2.elem[2][4], mat1.elem[2][5] * mat2.elem[2][5],
			mat1.elem[3][0] * mat2.elem[3][0], mat1.elem[3][1] * mat2.elem[3][1], mat1.elem[3][2] * mat2.elem[3][2], mat1.elem[3][3] * mat2.elem[3][3], mat1.elem[3][4] * mat2.elem[3][4], mat1.elem[3][5] * mat2.elem[3][5],
			mat1.elem[4][0] * mat2.elem[4][0], mat1.elem[4][1] * mat2.elem[4][1], mat1.elem[4][2] * mat2.elem[4][2], mat1.elem[4][3] * mat2.elem[4][3], mat1.elem[4][4] * mat2.elem[4][4], mat1.elem[4][5] * mat2.elem[4][5],
			mat1.elem[5][0] * mat2.elem[5][0], mat1.elem[5][1] * mat2.elem[5][1], mat1.elem[5][2] * mat2.elem[5][2], mat1.elem[5][3] * mat2.elem[5][3], mat1.elem[5][4] * mat2.elem[5][4], mat1.elem[5][5] * mat2.elem[5][5]
	};
	return tplmp::cast<matrix<6, 6, _Result> >(m_elem);
}

//7x7矩阵特化
//用于分析力学相空间计算
template<typename _T>
class matrix<7, 7, _T>
{
public:
	__def_matrix_members__(_T, 7, 7)

	__def_matrix_op_cast_ptr__(_T)

	__def_matrix_op_cast_vec_ptr__(_T, 49)

	__def_matrix_op_row_vec__(_T, 7)

	__def_matrix_op_partial_pivot_row__(_T, 7)

	__def_matrix_op_complete_pivot_coord__(_T, 7, 7)

	//交换行
	__def_matrix_op_swap_row__(_T, 7, 7)

	__def_matrix_op_sum__(_T, 7, 7)

	__def_matrix_op_partition__(_T, 7, 7)

	__def_matrix_op_mask__(_T, 7, 7)

	inline matrix<7,
	7, _T>& swap_column(size_t column1, size_t column2)
	{
		return *this;
	}

	template<typename _T2>
	inline operator matrix<7, 7, _T2>()
	{
		_T2 m_elem[] = {
				(_T2)elem[0][0], (_T2)elem[0][1], (_T2)elem[0][2], (_T2)elem[0][3], (_T2)elem[0][4], (_T2)elem[0][5], (_T2)elem[0][6],
				(_T2)elem[1][0], (_T2)elem[1][1], (_T2)elem[1][2], (_T2)elem[1][3], (_T2)elem[1][4], (_T2)elem[1][5], (_T2)elem[1][6],
				(_T2)elem[2][0], (_T2)elem[2][1], (_T2)elem[2][2], (_T2)elem[2][3], (_T2)elem[2][4], (_T2)elem[2][5], (_T2)elem[2][6],
				(_T2)elem[3][0], (_T2)elem[3][1], (_T2)elem[3][2], (_T2)elem[3][3], (_T2)elem[3][4], (_T2)elem[3][5], (_T2)elem[3][6],
				(_T2)elem[4][0], (_T2)elem[4][1], (_T2)elem[4][2], (_T2)elem[4][3], (_T2)elem[4][4], (_T2)elem[4][5], (_T2)elem[4][6],
				(_T2)elem[5][0], (_T2)elem[5][1], (_T2)elem[5][2], (_T2)elem[5][3], (_T2)elem[5][4], (_T2)elem[5][5], (_T2)elem[5][6],
				(_T2)elem[6][0], (_T2)elem[6][1], (_T2)elem[6][2], (_T2)elem[6][3], (_T2)elem[6][4], (_T2)elem[6][5], (_T2)elem[6][6]
		};
		return tplmp::cast<matrix<7, 7, _T2> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<7, 7, _Result> operator+(const matrix<7, 7, _T2>& mat) const
	{
		_Result m_elem[] = {
				elem[0][0] + mat.elem[0][0], elem[0][1] + mat.elem[0][1], elem[0][2] + mat.elem[0][2], elem[0][3] + mat.elem[0][3], elem[0][4] + mat.elem[0][4], elem[0][5] + mat.elem[0][5], elem[0][6] + mat.elem[0][6],
				elem[1][0] + mat.elem[1][0], elem[1][1] + mat.elem[1][1], elem[1][2] + mat.elem[1][2], elem[1][3] + mat.elem[1][3], elem[1][4] + mat.elem[1][4], elem[1][5] + mat.elem[1][5], elem[1][6] + mat.elem[1][6],
				elem[2][0] + mat.elem[2][0], elem[2][1] + mat.elem[2][1], elem[2][2] + mat.elem[2][2], elem[2][3] + mat.elem[2][3], elem[2][4] + mat.elem[2][4], elem[2][5] + mat.elem[2][5], elem[2][6] + mat.elem[2][6],
				elem[3][0] + mat.elem[3][0], elem[3][1] + mat.elem[3][1], elem[3][2] + mat.elem[3][2], elem[3][3] + mat.elem[3][3], elem[3][4] + mat.elem[3][4], elem[3][5] + mat.elem[3][5], elem[3][6] + mat.elem[3][6],
				elem[4][0] + mat.elem[4][0], elem[4][1] + mat.elem[4][1], elem[4][2] + mat.elem[4][2], elem[4][3] + mat.elem[4][3], elem[4][4] + mat.elem[4][4], elem[4][5] + mat.elem[4][5], elem[4][6] + mat.elem[4][6],
				elem[5][0] + mat.elem[5][0], elem[5][1] + mat.elem[5][1], elem[5][2] + mat.elem[5][2], elem[5][3] + mat.elem[5][3], elem[5][4] + mat.elem[5][4], elem[5][5] + mat.elem[5][5], elem[5][6] + mat.elem[5][6],
				elem[6][0] + mat.elem[6][0], elem[6][1] + mat.elem[6][1], elem[6][2] + mat.elem[6][2], elem[6][3] + mat.elem[6][3], elem[6][4] + mat.elem[6][4], elem[6][5] + mat.elem[6][5], elem[6][6] + mat.elem[6][6]
		};
		return tplmp::cast<matrix<7, 7, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<7, 7, _Result> operator-(const matrix<7, 7, _T2>& mat) const
	{
		_Result m_elem[] = {
				elem[0][0] - mat.elem[0][0], elem[0][1] - mat.elem[0][1], elem[0][2] - mat.elem[0][2], elem[0][3] - mat.elem[0][3], elem[0][4] - mat.elem[0][4], elem[0][5] - mat.elem[0][5], elem[0][6] - mat.elem[0][6],
				elem[1][0] - mat.elem[1][0], elem[1][1] - mat.elem[1][1], elem[1][2] - mat.elem[1][2], elem[1][3] - mat.elem[1][3], elem[1][4] - mat.elem[1][4], elem[1][5] - mat.elem[1][5], elem[1][6] - mat.elem[1][6],
				elem[2][0] - mat.elem[2][0], elem[2][1] - mat.elem[2][1], elem[2][2] - mat.elem[2][2], elem[2][3] - mat.elem[2][3], elem[2][4] - mat.elem[2][4], elem[2][5] - mat.elem[2][5], elem[2][6] - mat.elem[2][6],
				elem[3][0] - mat.elem[3][0], elem[3][1] - mat.elem[3][1], elem[3][2] - mat.elem[3][2], elem[3][3] - mat.elem[3][3], elem[3][4] - mat.elem[3][4], elem[3][5] - mat.elem[3][5], elem[3][6] - mat.elem[3][6],
				elem[4][0] - mat.elem[4][0], elem[4][1] - mat.elem[4][1], elem[4][2] - mat.elem[4][2], elem[4][3] - mat.elem[4][3], elem[4][4] - mat.elem[4][4], elem[4][5] - mat.elem[4][5], elem[4][6] - mat.elem[4][6],
				elem[5][0] - mat.elem[5][0], elem[5][1] - mat.elem[5][1], elem[5][2] - mat.elem[5][2], elem[5][3] - mat.elem[5][3], elem[5][4] - mat.elem[5][4], elem[5][5] - mat.elem[5][5], elem[5][6] - mat.elem[5][6],
				elem[6][0] - mat.elem[6][0], elem[6][1] - mat.elem[6][1], elem[6][2] - mat.elem[6][2], elem[6][3] - mat.elem[6][3], elem[6][4] - mat.elem[6][4], elem[6][5] - mat.elem[6][5], elem[6][6] - mat.elem[6][6]
		};
		return tplmp::cast<matrix<7, 7, _Result> >(m_elem);
	}

	inline matrix<7, 7, _T> operator-() const
	{
		_T m_elem[] = {
				-elem[0][0], -elem[0][1], -elem[0][2], -elem[0][3], -elem[0][4], -elem[0][5], -elem[0][6],
				-elem[1][0], -elem[1][1], -elem[1][2], -elem[1][3], -elem[1][4], -elem[1][5], -elem[1][6],
				-elem[2][0], -elem[2][1], -elem[2][2], -elem[2][3], -elem[2][4], -elem[2][5], -elem[2][6],
				-elem[3][0], -elem[3][1], -elem[3][2], -elem[3][3], -elem[3][4], -elem[3][5], -elem[3][6],
				-elem[4][0], -elem[4][1], -elem[4][2], -elem[4][3], -elem[4][4], -elem[4][5], -elem[4][6],
				-elem[5][0], -elem[5][1], -elem[5][2], -elem[5][3], -elem[5][4], -elem[5][5], -elem[5][6],
				-elem[6][0], -elem[6][1], -elem[6][2], -elem[6][3], -elem[6][4], -elem[6][5], -elem[6][6]
		};
		return tplmp::cast<matrix<7, 7, _T> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<7, 7, _Result> operator*(_T2 t) const
	{
		_Result m_elem[] = {
				elem[0][0] * t, elem[0][1] * t, elem[0][2] * t, elem[0][3] * t, elem[0][4] * t, elem[0][5] * t, elem[0][6] * t,
				elem[1][0] * t, elem[1][1] * t, elem[1][2] * t, elem[1][3] * t, elem[1][4] * t, elem[1][5] * t, elem[1][6] * t,
				elem[2][0] * t, elem[2][1] * t, elem[2][2] * t, elem[2][3] * t, elem[2][4] * t, elem[2][5] * t, elem[2][6] * t,
				elem[3][0] * t, elem[3][1] * t, elem[3][2] * t, elem[3][3] * t, elem[3][4] * t, elem[3][5] * t, elem[3][6] * t,
				elem[4][0] * t, elem[4][1] * t, elem[4][2] * t, elem[4][3] * t, elem[4][4] * t, elem[4][5] * t, elem[4][6] * t,
				elem[5][0] * t, elem[5][1] * t, elem[5][2] * t, elem[5][3] * t, elem[5][4] * t, elem[5][5] * t, elem[5][6] * t,
				elem[6][0] * t, elem[6][1] * t, elem[6][2] * t, elem[6][3] * t, elem[6][4] * t, elem[6][5] * t, elem[6][6] * t
		};
		return tplmp::cast<matrix<7, 7, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<7, 7, _T>& operator*=(const _T2& t)
	{
		elem[0] *= t;
		elem[1] *= t;
		elem[2] *= t;
		elem[3] *= t;
		elem[4] *= t;
		elem[5] *= t;
		elem[6] *= t;
		return *this;
	}

	__def_matrix_op_matrix_mul__(_T, 7, 7)

	__def_matrix_op_vector_mul__(_T, 7, 7)

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<7, 7, _Result> operator/(_T2 t) const
	{
		_Result m_elem[] =
		{
			elem[0][0] / t, elem[0][1] / t, elem[0][2] / t, elem[0][3] / t, elem[0][4] / t, elem[0][5] / t, elem[0][6] / t,
			elem[1][0] / t, elem[1][1] / t, elem[1][2] / t, elem[1][3] / t, elem[1][4] / t, elem[1][5] / t, elem[1][6] / t,
			elem[2][0] / t, elem[2][1] / t, elem[2][2] / t, elem[2][3] / t, elem[2][4] / t, elem[2][5] / t, elem[2][6] / t,
			elem[3][0] / t, elem[3][1] / t, elem[3][2] / t, elem[3][3] / t, elem[3][4] / t, elem[3][5] / t, elem[3][6] / t,
			elem[4][0] / t, elem[4][1] / t, elem[4][2] / t, elem[4][3] / t, elem[4][4] / t, elem[4][5] / t, elem[4][6] / t,
			elem[5][0] / t, elem[5][1] / t, elem[5][2] / t, elem[5][3] / t, elem[5][4] / t, elem[5][5] / t, elem[5][6] / t,
			elem[6][0] / t, elem[6][1] / t, elem[6][2] / t, elem[6][3] / t, elem[6][4] / t, elem[6][5] / t, elem[6][6] / t
		};
		return tplmp::cast<matrix<7, 7, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<7, 7, _T>& operator/=(const _T2& t)
	{
		elem[0] /= t;
		elem[1] /= t;
		elem[2] /= t;
		elem[3] /= t;
		elem[4] /= t;
		elem[5] /= t;
		elem[6] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<7, 7, _Result> operator/(const matrix<7, 7, _T2>& mat) const
	{
		_Result m_elem[] =
		{
			elem[0][0] / mat.elem[0][0], elem[0][1] / mat.elem[0][1], elem[0][2] / mat.elem[0][2], elem[0][3] / mat.elem[0][3], elem[0][4] / mat.elem[0][4], elem[0][5] / mat.elem[0][5], elem[0][6] / mat.elem[0][6],
			elem[1][0] / mat.elem[1][0], elem[1][1] / mat.elem[1][1], elem[1][2] / mat.elem[1][2], elem[1][3] / mat.elem[1][3], elem[1][4] / mat.elem[1][4], elem[1][5] / mat.elem[1][5], elem[1][6] / mat.elem[1][6],
			elem[2][0] / mat.elem[2][0], elem[2][1] / mat.elem[2][1], elem[2][2] / mat.elem[2][2], elem[2][3] / mat.elem[2][3], elem[2][4] / mat.elem[2][4], elem[2][5] / mat.elem[2][5], elem[2][6] / mat.elem[2][6],
			elem[3][0] / mat.elem[3][0], elem[3][1] / mat.elem[3][1], elem[3][2] / mat.elem[3][2], elem[3][3] / mat.elem[3][3], elem[3][4] / mat.elem[3][4], elem[3][5] / mat.elem[3][5], elem[3][6] / mat.elem[3][6],
			elem[4][0] / mat.elem[4][0], elem[4][1] / mat.elem[4][1], elem[4][2] / mat.elem[4][2], elem[4][3] / mat.elem[4][3], elem[4][4] / mat.elem[4][4], elem[4][5] / mat.elem[4][5], elem[4][6] / mat.elem[4][6],
			elem[5][0] / mat.elem[5][0], elem[5][1] / mat.elem[5][1], elem[5][2] / mat.elem[5][2], elem[5][3] / mat.elem[5][3], elem[5][4] / mat.elem[5][4], elem[5][5] / mat.elem[5][5], elem[5][6] / mat.elem[5][6],
			elem[6][0] / mat.elem[6][0], elem[6][1] / mat.elem[6][1], elem[6][2] / mat.elem[6][2], elem[6][3] / mat.elem[6][3], elem[6][4] / mat.elem[6][4], elem[6][5] / mat.elem[6][5], elem[6][6] / mat.elem[6][6]
		};
		return tplmp::cast<matrix<7, 7, _Result> >(m_elem);
	}

	inline vector<7, _T> column(size_t column) const
	{
		_T v_elem[] =
		{
			elem[0][column],
			elem[1][column],
			elem[2][column],
			elem[3][column],
			elem[4][column],
			elem[5][column],
			elem[6][column]
		};
		return tplmp::cast<vector<7, _T> >(v_elem);
	}

	inline matrix<7, 7, _T> operator~() const
	{
		_T m_elem[] =
		{
			elem[0][0], elem[1][0], elem[2][0], elem[3][0], elem[4][0], elem[5][0], elem[6][0],
			elem[0][1], elem[1][1], elem[2][1], elem[3][1], elem[4][1], elem[5][1], elem[6][1],
			elem[0][2], elem[1][2], elem[2][2], elem[3][2], elem[4][2], elem[5][2], elem[6][2],
			elem[0][3], elem[1][3], elem[2][3], elem[3][3], elem[4][3], elem[5][3], elem[6][3],
			elem[0][4], elem[1][4], elem[2][4], elem[3][4], elem[4][4], elem[5][4], elem[6][4],
			elem[0][5], elem[1][5], elem[2][5], elem[3][5], elem[4][5], elem[5][5], elem[6][5],
			elem[0][6], elem[1][6], elem[2][6], elem[3][6], elem[4][6], elem[5][6], elem[6][6]
		};
		return tplmp::cast<matrix<7, 7, _T> >(m_elem);
	}

	inline static matrix<7, 7, _T> all(_T value)
	{
		_T m_elem[] =
		{
			value, value, value, value, value, value, value,
			value, value, value, value, value, value, value,
			value, value, value, value, value, value, value,
			value, value, value, value, value, value, value,
			value, value, value, value, value, value, value,
			value, value, value, value, value, value, value,
			value, value, value, value, value, value, value
		};
		return tplmp::cast<matrix<7, 7, _T> >(m_elem);
	}
	__def_matrix_zero_one__(_T, 7, 7)

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << '\t' << elem[0][2] << '\t' << elem[0][3] << '\t' << elem[0][4] << '\t' << elem[0][5] << '\t' << elem[0][6] << "]\n"
		<< '[' << elem[1][0] << '\t' << elem[1][1] << '\t' << elem[1][2] << '\t' << elem[1][3] << '\t' << elem[1][4] << '\t' << elem[1][5] << '\t' << elem[1][6] << "]\n"
		<< '[' << elem[2][0] << '\t' << elem[2][1] << '\t' << elem[2][2] << '\t' << elem[2][3] << '\t' << elem[2][4] << '\t' << elem[2][5] << '\t' << elem[2][6] << "]\n"
		<< '[' << elem[3][0] << '\t' << elem[3][1] << '\t' << elem[3][2] << '\t' << elem[3][3] << '\t' << elem[3][4] << '\t' << elem[3][5] << '\t' << elem[3][6] << "]\n"
		<< '[' << elem[4][0] << '\t' << elem[4][1] << '\t' << elem[4][2] << '\t' << elem[4][3] << '\t' << elem[4][4] << '\t' << elem[4][5] << '\t' << elem[4][6] << "]\n"
		<< '[' << elem[5][0] << '\t' << elem[5][1] << '\t' << elem[5][2] << '\t' << elem[5][3] << '\t' << elem[5][4] << '\t' << elem[5][5] << '\t' << elem[5][6] << "]\n"
		<< '[' << elem[6][0] << '\t' << elem[6][1] << '\t' << elem[6][2] << '\t' << elem[6][3] << '\t' << elem[6][4] << '\t' << elem[6][5] << '\t' << elem[6][6] << ']';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<7, 7, _Result> operator*(_T1 t, const matrix<7, 7, _T2>& mat)
{
	_Result m_elem[] = {
			t * mat.elem[0][0], t * mat.elem[0][1], t * mat.elem[0][2], t * mat.elem[0][3], t * mat.elem[0][4], t * mat.elem[0][5], t * mat.elem[0][6],
			t * mat.elem[1][0], t * mat.elem[1][1], t * mat.elem[1][2], t * mat.elem[1][3], t * mat.elem[1][4], t * mat.elem[1][5], t * mat.elem[1][6],
			t * mat.elem[2][0], t * mat.elem[2][1], t * mat.elem[2][2], t * mat.elem[2][3], t * mat.elem[2][4], t * mat.elem[2][5], t * mat.elem[2][6],
			t * mat.elem[3][0], t * mat.elem[3][1], t * mat.elem[3][2], t * mat.elem[3][3], t * mat.elem[3][4], t * mat.elem[3][5], t * mat.elem[3][6],
			t * mat.elem[4][0], t * mat.elem[4][1], t * mat.elem[4][2], t * mat.elem[4][3], t * mat.elem[4][4], t * mat.elem[4][5], t * mat.elem[4][6],
			t * mat.elem[5][0], t * mat.elem[5][1], t * mat.elem[5][2], t * mat.elem[5][3], t * mat.elem[5][4], t * mat.elem[5][5], t * mat.elem[5][6],
			t * mat.elem[6][0], t * mat.elem[6][1], t * mat.elem[6][2], t * mat.elem[6][3], t * mat.elem[6][4], t * mat.elem[6][5], t * mat.elem[6][6]
	};
	return tplmp::cast<matrix<7, 7, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<7, 7, _Result> operator/(_T1 t, const matrix<7, 7, _T2>& mat)
{
	_Result m_elem[] = {
			t / mat.elem[0][0], t / mat.elem[0][1], t / mat.elem[0][2], t / mat.elem[0][3], t / mat.elem[0][4], t / mat.elem[0][5], t / mat.elem[0][6],
			t / mat.elem[1][0], t / mat.elem[1][1], t / mat.elem[1][2], t / mat.elem[1][3], t / mat.elem[1][4], t / mat.elem[1][5], t / mat.elem[1][6],
			t / mat.elem[2][0], t / mat.elem[2][1], t / mat.elem[2][2], t / mat.elem[2][3], t / mat.elem[2][4], t / mat.elem[2][5], t / mat.elem[2][6],
			t / mat.elem[3][0], t / mat.elem[3][1], t / mat.elem[3][2], t / mat.elem[3][3], t / mat.elem[3][4], t / mat.elem[3][5], t / mat.elem[3][6],
			t / mat.elem[4][0], t / mat.elem[4][1], t / mat.elem[4][2], t / mat.elem[4][3], t / mat.elem[4][4], t / mat.elem[4][5], t / mat.elem[4][6],
			t / mat.elem[5][0], t / mat.elem[5][1], t / mat.elem[5][2], t / mat.elem[5][3], t / mat.elem[5][4], t / mat.elem[5][5], t / mat.elem[5][6],
			t / mat.elem[6][0], t / mat.elem[6][1], t / mat.elem[6][2], t / mat.elem[6][3], t / mat.elem[6][4], t / mat.elem[6][5], t / mat.elem[6][6]
	};
	return tplmp::cast<matrix<7, 7, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<7, 7, _Result> hadamard(const matrix<7, 7, _T1>& mat1, const matrix<7, 7, _T2>& mat2)
{
	_Result m_elem[] = {
			mat1.elem[0][0] * mat2.elem[0][0], mat1.elem[0][1] * mat2.elem[0][1], mat1.elem[0][2] * mat2.elem[0][2], mat1.elem[0][3] * mat2.elem[0][3], mat1.elem[0][4] * mat2.elem[0][4], mat1.elem[0][5] * mat2.elem[0][5], mat1.elem[0][6] * mat2.elem[0][6],
			mat1.elem[1][0] * mat2.elem[1][0], mat1.elem[1][1] * mat2.elem[1][1], mat1.elem[1][2] * mat2.elem[1][2], mat1.elem[1][3] * mat2.elem[1][3], mat1.elem[1][4] * mat2.elem[1][4], mat1.elem[1][5] * mat2.elem[1][5], mat1.elem[1][6] * mat2.elem[1][6],
			mat1.elem[2][0] * mat2.elem[2][0], mat1.elem[2][1] * mat2.elem[2][1], mat1.elem[2][2] * mat2.elem[2][2], mat1.elem[2][3] * mat2.elem[2][3], mat1.elem[2][4] * mat2.elem[2][4], mat1.elem[2][5] * mat2.elem[2][5], mat1.elem[2][6] * mat2.elem[2][6],
			mat1.elem[3][0] * mat2.elem[3][0], mat1.elem[3][1] * mat2.elem[3][1], mat1.elem[3][2] * mat2.elem[3][2], mat1.elem[3][3] * mat2.elem[3][3], mat1.elem[3][4] * mat2.elem[3][4], mat1.elem[3][5] * mat2.elem[3][5], mat1.elem[3][6] * mat2.elem[3][6],
			mat1.elem[4][0] * mat2.elem[4][0], mat1.elem[4][1] * mat2.elem[4][1], mat1.elem[4][2] * mat2.elem[4][2], mat1.elem[4][3] * mat2.elem[4][3], mat1.elem[4][4] * mat2.elem[4][4], mat1.elem[4][5] * mat2.elem[4][5], mat1.elem[4][6] * mat2.elem[4][6],
			mat1.elem[5][0] * mat2.elem[5][0], mat1.elem[5][1] * mat2.elem[5][1], mat1.elem[5][2] * mat2.elem[5][2], mat1.elem[5][3] * mat2.elem[5][3], mat1.elem[5][4] * mat2.elem[5][4], mat1.elem[5][5] * mat2.elem[5][5], mat1.elem[5][6] * mat2.elem[5][6],
			mat1.elem[6][0] * mat2.elem[6][0], mat1.elem[6][1] * mat2.elem[6][1], mat1.elem[6][2] * mat2.elem[6][2], mat1.elem[6][3] * mat2.elem[6][3], mat1.elem[6][4] * mat2.elem[6][4], mat1.elem[6][5] * mat2.elem[6][5], mat1.elem[6][6] * mat2.elem[6][6]
	};
	return tplmp::cast<matrix<7, 7, _Result> >(m_elem);
}

//1x6矩阵特化（行向量）
//用于分析力学相空间计算
template<typename _T>
class matrix<1, 6, _T>
{
public:
	__def_matrix_members__(_T, 1, 6)

	__def_matrix_op_cast_ptr__(_T)

	__def_matrix_op_cast_vec_ptr__(_T, 6)

	__def_matrix_op_row_vec__(_T, 6)

	__def_matrix_op_partial_pivot_row__(_T, 6)

	__def_matrix_op_complete_pivot_coord__(_T, 1, 6)

	//交换行
	__def_matrix_op_swap_row__(_T, 1, 6)

	__def_matrix_op_sum__(_T, 1, 6)

	__def_matrix_op_partition__(_T, 1, 6)

	__def_matrix_op_mask__(_T, 1, 6)

	inline matrix<1,
	6, _T>& swap_column(size_t column1, size_t column2)
	{
		return *this;
	}

	template<typename _T2>
	inline operator matrix<1, 6, _T2>()
	{
		_T2 m_elem[] = {
				(_T2)elem[0][0], (_T2)elem[0][1], (_T2)elem[0][2], (_T2)elem[0][3], (_T2)elem[0][4], (_T2)elem[0][5]
		};
		return tplmp::cast<matrix<1, 6, _T2> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<1, 6, _Result> operator+(const matrix<1, 6, _T2>& mat) const
	{
		_Result m_elem[] = {
				elem[0][0] + mat.elem[0][0],
				elem[0][1] + mat.elem[0][1],
				elem[0][2] + mat.elem[0][2],
				elem[0][3] + mat.elem[0][3],
				elem[0][4] + mat.elem[0][4],
				elem[0][5] + mat.elem[0][5]
		};
		return tplmp::cast<matrix<1, 6, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<1, 6, _Result> operator-(const matrix<1, 6, _T2>& mat) const
	{
		_Result m_elem[] = {
				elem[0][0] - mat.elem[0][0],
				elem[0][1] - mat.elem[0][1],
				elem[0][2] - mat.elem[0][2],
				elem[0][3] - mat.elem[0][3],
				elem[0][4] - mat.elem[0][4],
				elem[0][5] - mat.elem[0][5]
		};
		return tplmp::cast<matrix<1, 6, _Result> >(m_elem);
	}

	inline matrix<1, 6, _T> operator-() const
	{
		_T m_elem[] = {
				-elem[0][0],
				-elem[0][1],
				-elem[0][2],
				-elem[0][3],
				-elem[0][4],
				-elem[0][5]
		};
		return tplmp::cast<matrix<1, 6, _T> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<1, 6, _Result> operator*(_T2 t) const
	{
		_Result m_elem[] = {
				elem[0][0] * t,
				elem[0][1] * t,
				elem[0][2] * t,
				elem[0][3] * t,
				elem[0][4] * t,
				elem[0][5] * t
		};
		return tplmp::cast<matrix<1, 6, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<1, 6, _T>& operator*=(const _T2& t)
	{
		elem[0] *= t;
		return *this;
	}

	__def_matrix_op_matrix_mul__(_T, 1, 6)

	__def_matrix_op_vector_mul__(_T, 1, 6)

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<1, 6, _Result> operator/(_T2 t) const
	{
		_Result m_elem[] =
		{
			elem[0][0] / t,
			elem[0][1] / t,
			elem[0][2] / t,
			elem[0][3] / t,
			elem[0][4] / t,
			elem[0][5] / t
		};
		return tplmp::cast<matrix<1, 6, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<1, 6, _T>& operator/=(const _T2& t)
	{
		elem[0] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<1, 6, _Result> operator/(const matrix<1, 6, _T2>& mat) const
	{
		_Result m_elem[] =
		{
			elem[0][0] / mat.elem[0][0],
			elem[0][1] / mat.elem[0][1],
			elem[0][2] / mat.elem[0][2],
			elem[0][3] / mat.elem[0][3],
			elem[0][4] / mat.elem[0][4],
			elem[0][5] / mat.elem[0][5]
		};
		return tplmp::cast<matrix<1, 6, _Result> >(m_elem);
	}

	inline vector<1, _T> column(size_t column) const
	{
		_T v_elem[] =
		{
			elem[0][column]
		};
		return tplmp::cast<vector<1, _T> >(v_elem);
	}

	inline matrix<6, 1, _T> operator~() const
	{
		_T m_elem[] =
		{
			elem[0][0],
			elem[0][1],
			elem[0][2],
			elem[0][3],
			elem[0][4],
			elem[0][5]
		};
		return tplmp::cast<matrix<6, 1, _T> >(m_elem);
	}

	inline static matrix<1, 6, _T> all(_T value)
	{
		_T m_elem[] =
		{
			value, value, value, value, value, value
		};
		return tplmp::cast<matrix<1, 6, _T> >(m_elem);
	}
	__def_matrix_zero_one__(_T, 1, 6)

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << '\t' << elem[0][2] << '\t' << elem[0][3] << '\t' << elem[0][4] << '\t' << elem[0][5] << ']';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<1, 6, _Result> operator*(_T1 t, const matrix<1, 6, _T2>& mat)
{
	_Result m_elem[] = {
			t * mat.elem[0][0], t * mat.elem[0][1], t * mat.elem[0][2],
			t * mat.elem[0][3], t * mat.elem[0][4], t * mat.elem[0][5]
	};
	return tplmp::cast<matrix<1, 6, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<1, 6, _Result> operator/(_T1 t, const matrix<1, 6, _T2>& mat)
{
	_Result m_elem[] = {
			t / mat.elem[0][0], t / mat.elem[0][1], t / mat.elem[0][2],
			t / mat.elem[0][3], t / mat.elem[0][4], t / mat.elem[0][5]
	};
	return tplmp::cast<matrix<1, 6, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<1, 6, _Result> hadamard(const matrix<1, 6, _T1>& mat1, const matrix<1, 6, _T2>& mat2)
{
	_Result m_elem[] = {
			mat1.elem[0][0] * mat2.elem[0][0],
			mat1.elem[0][1] * mat2.elem[0][1],
			mat1.elem[0][2] * mat2.elem[0][2],
			mat1.elem[0][3] * mat2.elem[0][3],
			mat1.elem[0][4] * mat2.elem[0][4],
			mat1.elem[0][5] * mat2.elem[0][5]
	};
	return tplmp::cast<matrix<1, 6, _Result> >(m_elem);
}

//1x7矩阵特化（行向量）
//用于分析力学相空间计算
template<typename _T>
class matrix<1, 7, _T>
{
public:
	__def_matrix_members__(_T, 1, 7)

	__def_matrix_op_cast_ptr__(_T)

	__def_matrix_op_cast_vec_ptr__(_T, 7)

	__def_matrix_op_row_vec__(_T, 7)

	__def_matrix_op_partial_pivot_row__(_T, 7)

	__def_matrix_op_complete_pivot_coord__(_T, 1, 7)

	//交换行
	__def_matrix_op_swap_row__(_T, 1, 7)

	__def_matrix_op_sum__(_T, 1, 7)

	__def_matrix_op_partition__(_T, 1, 7)

	__def_matrix_op_mask__(_T, 1, 7)

	inline matrix<1,
	7, _T>& swap_column(size_t column1, size_t column2)
	{
		return *this;
	}

	template<typename _T2>
	inline operator matrix<1, 7, _T2>()
	{
		_T2 m_elem[] = {
				(_T2)elem[0][0], (_T2)elem[0][1], (_T2)elem[0][2], (_T2)elem[0][3], (_T2)elem[0][4], (_T2)elem[0][5], (_T2)elem[0][6]
		};
		return tplmp::cast<matrix<1, 7, _T2> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<1, 7, _Result> operator+(const matrix<1, 7, _T2>& mat) const
	{
		_Result m_elem[] = {
				elem[0][0] + mat.elem[0][0],
				elem[0][1] + mat.elem[0][1],
				elem[0][2] + mat.elem[0][2],
				elem[0][3] + mat.elem[0][3],
				elem[0][4] + mat.elem[0][4],
				elem[0][5] + mat.elem[0][5],
				elem[0][6] + mat.elem[0][6]
		};
		return tplmp::cast<matrix<1, 7, _Result> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<1, 7, _Result> operator-(const matrix<1, 7, _T2>& mat) const
	{
		_Result m_elem[] = {
				elem[0][0] - mat.elem[0][0],
				elem[0][1] - mat.elem[0][1],
				elem[0][2] - mat.elem[0][2],
				elem[0][3] - mat.elem[0][3],
				elem[0][4] - mat.elem[0][4],
				elem[0][5] - mat.elem[0][5],
				elem[0][6] - mat.elem[0][6]
		};
		return tplmp::cast<matrix<1, 7, _Result> >(m_elem);
	}

	inline matrix<1, 7, _T> operator-() const
	{
		_T m_elem[] = {
				-elem[0][0],
				-elem[0][1],
				-elem[0][2],
				-elem[0][3],
				-elem[0][4],
				-elem[0][5],
				-elem[0][6]
		};
		return tplmp::cast<matrix<1, 7, _T> >(m_elem);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<1, 7, _Result> operator*(_T2 t) const
	{
		_Result m_elem[] = {
				elem[0][0] * t,
				elem[0][1] * t,
				elem[0][2] * t,
				elem[0][3] * t,
				elem[0][4] * t,
				elem[0][5] * t,
				elem[0][6] * t
		};
		return tplmp::cast<matrix<1, 7, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<1, 7, _T>& operator*=(const _T2& t)
	{
		elem[0] *= t;
		return *this;
	}

	__def_matrix_op_matrix_mul__(_T, 1, 7)

	__def_matrix_op_vector_mul__(_T, 1, 7)

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<1, 7, _Result> operator/(_T2 t) const
	{
		_Result m_elem[] =
		{
			elem[0][0] / t,
			elem[0][1] / t,
			elem[0][2] / t,
			elem[0][3] / t,
			elem[0][4] / t,
			elem[0][5] / t,
			elem[0][6] / t
		};
		return tplmp::cast<matrix<1, 7, _Result> >(m_elem);
	}

	template<typename _T2>
	inline matrix<1, 7, _T>& operator/=(const _T2& t)
	{
		elem[0] /= t;
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<1, 7, _Result> operator/(const matrix<1, 7, _T2>& mat) const
	{
		_Result m_elem[] =
		{
			elem[0][0] / mat.elem[0][0],
			elem[0][1] / mat.elem[0][1],
			elem[0][2] / mat.elem[0][2],
			elem[0][3] / mat.elem[0][3],
			elem[0][4] / mat.elem[0][4],
			elem[0][5] / mat.elem[0][5],
			elem[0][6] / mat.elem[0][6]
		};
		return tplmp::cast<matrix<1, 7, _Result> >(m_elem);
	}

	inline vector<1, _T> column(size_t column) const
	{
		_T v_elem[] =
		{
			elem[0][column]
		};
		return tplmp::cast<vector<1, _T> >(v_elem);
	}

	inline matrix<7, 1, _T> operator~() const
	{
		_T m_elem[] =
		{
			elem[0][0],
			elem[0][1],
			elem[0][2],
			elem[0][3],
			elem[0][4],
			elem[0][5],
			elem[0][6]
		};
		return tplmp::cast<matrix<7, 1, _T> >(m_elem);
	}

	inline static matrix<1, 7, _T> all(_T value)
	{
		_T m_elem[] =
		{
			value, value, value, value, value, value, value
		};
		return tplmp::cast<matrix<1, 7, _T> >(m_elem);
	}
	__def_matrix_zero_one__(_T, 1, 7)

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << '[' << elem[0][0] << '\t' << elem[0][1] << '\t' << elem[0][2] << '\t' << elem[0][3] << '\t' << elem[0][4] << '\t' << elem[0][5] << '\t' << elem[0][6] << ']';
		return oss.str();
	}
};

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<1, 7, _Result> operator*(_T1 t, const matrix<1, 7, _T2>& mat)
{
	_Result m_elem[] = {
			t * mat.elem[0][0], t * mat.elem[0][1], t * mat.elem[0][2],
			t * mat.elem[0][3], t * mat.elem[0][4], t * mat.elem[0][5],
			t * mat.elem[0][6]
	};
	return tplmp::cast<matrix<1, 7, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<1, 7, _Result> operator/(_T1 t, const matrix<1, 7, _T2>& mat)
{
	_Result m_elem[] = {
			t / mat.elem[0][0], t / mat.elem[0][1], t / mat.elem[0][2],
			t / mat.elem[0][3], t / mat.elem[0][4], t / mat.elem[0][5],
			t / mat.elem[0][6]
	};
	return tplmp::cast<matrix<1, 7, _Result> >(m_elem);
}

template<typename _T1, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() * tplmp::decl<_T2>::val())>
inline matrix<1, 7, _Result> hadamard(const matrix<1, 7, _T1>& mat1, const matrix<1, 7, _T2>& mat2)
{
	_Result m_elem[] = {
			mat1.elem[0][0] * mat2.elem[0][0],
			mat1.elem[0][1] * mat2.elem[0][1],
			mat1.elem[0][2] * mat2.elem[0][2],
			mat1.elem[0][3] * mat2.elem[0][3],
			mat1.elem[0][4] * mat2.elem[0][4],
			mat1.elem[0][5] * mat2.elem[0][5],
			mat1.elem[0][6] * mat2.elem[0][6]
	};
	return tplmp::cast<matrix<1, 7, _Result> >(m_elem);
}

template<typename _T>
using matrix1x1 = matrix<1, 1, _T>;
template<typename _T>
using matrix1x6 = matrix<1, 6, _T>;
template<typename _T>
using matrix1x7 = matrix<1, 7, _T>;
template<typename _T>
using matrix2x2 = matrix<2, 2, _T>;
template<typename _T>
using matrix3x3 = matrix<3, 3, _T>;
template<typename _T>
using matrix4x4 = matrix<4, 4, _T>;
template<typename _T>
using matrix3x4 = matrix<3, 4, _T>;
template<typename _T>
using matrix6x6 = matrix<6, 6, _T>;
template<typename _T>
using matrix7x7 = matrix<7, 7, _T>;

using matrix1x1f = matrix1x1<float>;
using matrix1x6f = matrix1x6<float>;
using matrix1x7f = matrix1x7<float>;
using matrix2x2f = matrix2x2<float>;
using matrix3x3f = matrix3x3<float>;
using matrix3x4f = matrix3x4<float>;
using matrix4x4f = matrix4x4<float>;
using matrix6x6f = matrix6x6<float>;
using matrix7x7f = matrix7x7<float>;

using matrix1x1d = matrix1x1<double>;
using matrix1x6d = matrix1x6<double>;
using matrix1x7d = matrix1x7<double>;
using matrix2x2d = matrix2x2<double>;
using matrix3x3d = matrix3x3<double>;
using matrix3x4d = matrix3x4<double>;
using matrix4x4d = matrix4x4<double>;
using matrix6x6d = matrix6x6<double>;
using matrix7x7d = matrix7x7<double>;

using matrix1x1q = matrix1x1<__float128>;
using matrix1x6q = matrix1x6<__float128>;
using matrix1x7q = matrix1x7<__float128>;
using matrix2x2q = matrix2x2<__float128>;
using matrix3x3q = matrix3x3<__float128>;
using matrix3x4q = matrix3x4<__float128>;
using matrix4x4q = matrix4x4<__float128>;
using matrix6x6q = matrix6x6<__float128>;
using matrix7x7q = matrix7x7<__float128>;
}
#endif//_MATH_MATRIX
