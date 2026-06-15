#ifndef _MATH_MATRIX
#define _MATH_MATRIX

#include <tplmp/base.h>

#include <sstream>

#include <math/vector.h>
#include <math/algebra.h>

#define __matrix_loop_unroll__() __loop_unroll__(9)

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
	math::vector<_Column, _T> elem[_Row];
	typedef _T type;
	static const size_t row_num = _Row;
	static const size_t column_num = _Column;

	inline operator _T*()
	{
		return elem;
	}

	template<size_t _Row2, size_t _Column2>
	inline operator math::matrix<_Row2, _Column2, _T>&()
	{
		return *(math::matrix<_Row2, _Column2, _T>*)this;
	}

	inline operator math::vector<_Column, _T>*()
	{
		return (math::vector<_Column, _T>*)elem;
	}

	/**
	 * 提取行向量
	 */
	inline math::vector<_Column, _T>& operator[](size_t row)
	{
		return elem[row];
	}

	inline const math::vector<_Column, _T>& operator[](size_t row) const
	{
		return elem[row];
	}

	inline math::vector<_Column, _T> row(size_t row)
	{
		return elem[row];
	}

	/**
	 * 找列主元。
	 * 如果找到则返回行索引，未找到（全0）则返回row
	 */
	size_t partial_pivot_row(size_t i, size_t j, _T eps)
	{
		size_t pivot_row = _Row;
		_T max_abs = _T(0);
		for(size_t r = i; r < _Row; ++r)
		{
			_T abs_val = math::abs(elem[r][j]);
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
	vector2<size_t> complete_pivot_coord(size_t i, size_t j, _T eps)
	{
		vector2<size_t> pivot_coord = {_Row, _Column};
		_T max_abs = _T(0);
		for(size_t r = i; r < _Row; ++r)
		{
			for(size_t c = j; c < _Column; ++c)
			{
				_T abs_val = math::abs(elem[r][c]);
				if(abs_val > max_abs)
				{
					max_abs = abs_val;
					pivot_coord[0] = r;
					pivot_coord[1] = c;
				}
			}
		}
		if(math::is_zero(max_abs, eps))
			return
			{	_Row, _Column};
		else
			return pivot_coord;
	}

	/**
	 * 交换本矩阵的行
	 */
	inline matrix<_Row, _Column, _T> swap_row(size_t row1, size_t row2)
	{
		vector<_Column, _T> temp = elem[row1];
		elem[row1] = elem[row2];
		elem[row2] = temp;
		return *this;
	}

	/**
	 * 交换本矩阵的列
	 */
	inline matrix<_Row, _Column, _T> swap_column(size_t column1, size_t column2)
	{
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
		{
			_T temp = elem[i][column1];
			elem[i][column1] = elem[i][column2];
			elem[i][column2] = temp;
		}
		return *this;
	}

	/**
	 * 矩阵分块，分块区间为[_BeginRow, _EndRow)，[_BeginColumn, _EndColumn)
	 */
	template<size_t _PartitionRow, size_t _PartitionColumn>
	inline matrix<_PartitionRow, _PartitionColumn, _T> partition(size_t row_begin, size_t column_begin) const
	{
		static_assert(_PartitionRow >= 0 && _PartitionColumn >= 0, "invalid partition size");
		matrix<_PartitionRow, _PartitionColumn, _T> partition_result;
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _PartitionRow; ++i)
		{
			partition_result.elem[i] = elem[row_begin + i].template slice<_PartitionColumn>(column_begin);
		}
		return partition_result;
	}

	template<typename _Result = _T>
	inline matrix<_Row, _Column, _Result> mask(size_t row_begin, size_t row_end, size_t column_begin, size_t column_end, _Result mask_value) const
	{
		matrix<_Row, _Column, _Result> mask_mat;
		for(size_t i = 0; i < row_begin; ++i)
		{
			mask_mat.elem[i] = vector<_Column, _Result>::all(mask_value);
		}
		for(size_t i = row_begin; i < row_end; ++i)
		{
			mask_mat.elem[i] = elem[i].template mask<_Result>(column_begin, column_end, mask_value);
		}
		for(size_t i = row_end; i < _Row; ++i)
		{
			mask_mat.elem[i] = vector<_Column, _Result>::all(mask_value);
		}
		return mask_mat;
	}

	/**
	 * 提取列向量
	 */
	inline math::vector<_Row, _T> column(size_t column)
	{
		vector<_Row, _T> row_vec;
		__matrix_loop_unroll__()
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
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
			cast_result[i] = (_T2)elem[i];
		return cast_result;
	}

	template<typename _T2>
	inline matrix<_Row, _Column, _T>& operator=(const matrix<_Row, _Column, _T2>& mat)
	{
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
			elem[i] = mat.elem[i];
		return *this;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator+(const matrix<_Row, _Column, _T2>& mat) const
	{
		matrix<_Row, _Column, _Result> add_result;
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
			add_result[i] = elem[i] + mat.elem[i];
		return add_result;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator-(const matrix<_Row, _Column, _T2>& mat) const
	{
		matrix<_Row, _Column, decltype(tplmp::decl<_T>::val()-tplmp::decl<_T2>::val())> sub_result;
		__matrix_loop_unroll__()
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
		__matrix_loop_unroll__()
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
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
			mul_result[i] = elem[i] * t;
		return mul_result;
	}

	template<typename _T2>
	inline matrix<_Row, _Column, _T>& operator*=(const _T2& t)
	{
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
			elem[i] *= t;
		return *this;
	}

	/**
	 * 矩阵乘法
	 */
	template<size_t _Row2, typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Row2, _Result> operator*(const matrix<_Column, _Row2, _T2>& mat) const
	{
		matrix<_Row, _Row2, _Result> mul_result;
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
		{
			__matrix_loop_unroll__()
			for(size_t j = 0; j < _Row2; ++j)
				mul_result[i][j] = this->row(i) * mat.column(j);
		}
		return mul_result;
	}

	/**
	 * 向量乘法
	 */
	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	inline math::vector<_Row, _Result> operator*(const math::vector<_Column, _T2>& vec) const
	{
		vector<_Row, _Result> mul_result;
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
			mul_result[i] = row(i) * vec;
		return mul_result;
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() / tplmp::decl<_T2>::val())>
	inline matrix<_Row, _Column, _Result> operator/(_T2 t) const
	{
		matrix<_Row, _Column, _Result> div_result;
		__matrix_loop_unroll__()
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
		__matrix_loop_unroll__()
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
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Column; ++i)
		{
			__matrix_loop_unroll__()
			for(size_t j = 0; j < _Row; ++j)
				transpose_result[i][j] = elem[j][i];
		}
		return transpose_result;
	}

	inline static matrix<_Row, _Column, _T> all(_T value)
	{
		matrix<_Row, _Column, _T> mat;
		mat.elem[0] = vector<_Column, _T>::all(value);
		__matrix_loop_unroll__()
		for(size_t i = 1; i < _Row; ++i)
			mat.elem[i] = mat.elem[0];
		return mat;
	}

	inline static matrix<_Row, _Column, _T> zero()
	{
		return all(0);
	}

	inline static matrix<_Row, _Column, _T> one()
	{
		return all(1);
	}

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum_row(size_t i) const
	{
		_T sum(0);
		__matrix_loop_unroll__()
		for(size_t j = 0; j < _Column; ++j)
		{
			sum += elem[i][j];
		}
		return sum;
	}

	template<typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T>::val())>
	inline _Result sum_column(size_t j) const
	{
		_T sum(0);
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Row; ++i)
		{
			sum += elem[i][j];
		}
		return sum;
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
					oss << 't';
			}
			oss << ']';
			if(i < _Row - 1)
				oss << 'n';
		}
		return oss.str();
	}
};

/**
 * @brief 矩阵子视图，可以直接在其上进行读写
 */
template<size_t _Row, size_t _Column, typename _T>
class matrix_view
{
	_T* elem[_Row];
	typedef _T type;
	static const size_t row_num = _Row;
	static const size_t column_num = _Column;

	matrix_view() = default;

	template<size_t _Row2, size_t _Column2>
	matrix_view(const matrix<_Row2, _Column2, _T>& target, size_t row_begin, size_t column_begin)
	{
		static_assert(_Row2 >= _Row && _Column2 >= _Column, "matrix view dim should be less than target");
		for(size_t i = 0; i < _Row; ++i)
		{
			//必须确保column_begin+_Column<=_Column2
			elem[i] = target.elem[row_begin + i].coord + column_begin;
		}
	}

	template<typename _T2>
	inline matrix_view<_Row, _Column, _T>& operator=(const matrix<_Row, _Column, _T2>& mat)
	{
		for(size_t i = 0; i < _Row; ++i)
			*(vector<_Column, _T> )elem[i] = mat.elem[i];
		return *this;
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
	__matrix_loop_unroll__()
	for(size_t i = 0; i < _Row; ++i)
		mul_result[i] = t * mat.elem[i];
	return mul_result;
}

template<typename _T1, size_t _Row, size_t _Column, typename _T2, typename _Result = decltype(tplmp::decl<_T1>::val() / tplmp::decl<_T2>::val())>
inline matrix<_Row, _Column, _Result> operator/(_T1 t, const matrix<_Row, _Column, _T2>& mat)
{
	matrix<_Row, _Column, _Result> div_result;
	__matrix_loop_unroll__()
	for(size_t i = 0; i < _Row; ++i)
		div_result[i] = t / mat.elem[i];
	return div_result;
}

template<size_t _Row, size_t _Column1, size_t _Column2, typename _T>
inline matrix<_Row, _Column1 + _Column2, _T> cat_column(const matrix<_Row, _Column1, _T>& mat1, const matrix<_Row, _Column2, _T>& mat2)
{
	matrix<_Row, _Column1 + _Column2, _T> mat;
	__matrix_loop_unroll__()
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
	__matrix_loop_unroll__()
	for(size_t i = 0; i < _Row1; ++i)
		mat.elem[i] = cat(left_up.elem[i], right_up.elem[i]);
	__matrix_loop_unroll__()
	for(size_t i = _Row1; i < cat_row; ++i)
		mat.elem[i] = cat(left_down.elem[i], right_down.elem[i]);
	return mat;
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
 * @brief REF回代法求解线性方程组
 * 		  即在已经矩阵A、向量b的情况下，求满足Ax=b的x向量
 */
template<size_t _Order, typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
vector<_Order, _Result> solve_linear_system(const matrix<_Order, _Order, _T>& coeff_mat, const vector<_Order, _T>& vals, bool* success, _T eps = 0)
{
	//将系数矩阵和值向量横向拼接成增广矩阵
	matrix<_Order, _Order + 1, _T> augmented;
	for(size_t i = 0; i < _Order; ++i)
	{
		for(size_t j = 0; j < _Order; ++j)
			augmented[i][j] = coeff_mat[i][j];
		augmented[i][_Order] = vals[i];
	}
	vector<_Order, _Result> solution;
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
template<size_t _Row, size_t _Column, typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
size_t rank(const matrix<_Row, _Column, _T>& mat, _T eps = 0)
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
_Result det(const matrix<_Order, _Order, _T>& mat, _T eps = 0)
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
	__matrix_loop_unroll__()
	for(size_t i = 0; i < _Row; ++i)
	{
		__matrix_loop_unroll__()
		for(size_t j = 0; j < _Column; ++j)
			hadamard_result[i][j] = mat1.elem[i][j] * mat2.elem[i][j];
	}
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
	__matrix_loop_unroll__()
	for(size_t i = 0; i < _Order; ++i)
		tr += mat.elem[i][i];
	return tr;
}

/**
 * 恒等元
 */
//加法恒等元：零矩阵
template<size_t _Row, size_t _Column, typename _T>
struct __identity_impl<matrix<_Row, _Column, _T>, add>
{
	inline static matrix<_Row, _Column, _T> value()
	{
		return matrix<_Row, _Column, _T>::zero();
	}
};
//乘法恒等元：恒等矩阵
template<size_t _Order, typename _T>
struct __identity_impl<matrix<_Order, _Order, _T>, mul>
{
	inline static matrix<_Order, _Order, _T> value()
	{
		matrix<_Order, _Order, _T> identity = matrix<_Order, _Order, _T>::zero();
		__matrix_loop_unroll__()
		for(size_t i = 0; i < _Order; ++i)
			identity[i][i] = 1;
		return identity;
	}
};
//乘法零元：零矩阵
template<size_t _Row, size_t _Column, typename _T>
struct __zero_impl<matrix<_Row, _Column, _T>, mul>
{
	inline static matrix<_Row, _Column, _T> value()
	{
		return matrix<_Row, _Column, _T>::zero();
	}
};

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const matrix<1, 1, _T>& mat)
{
	return mat.elem[0][0];
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const matrix<2, 2, _T>& mat)
{
	return mat.elem[0][0] * mat.elem[1][1] - mat.elem[0][1] * mat.elem[1][0];
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const matrix<3, 3, _T>& mat)
{
	return mat.elem[0][0] * (mat.elem[1][1] * mat.elem[2][2] - mat.elem[1][2] * mat.elem[2][1])
			+ mat.elem[0][1] * (mat.elem[1][2] * mat.elem[2][0] - mat.elem[1][0] * mat.elem[2][2])
			+ mat.elem[0][2] * (mat.elem[1][0] * mat.elem[2][1] - mat.elem[1][1] * mat.elem[2][0]);
}

template<typename _T, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T>::val())>
inline _Result det(const matrix<4, 4, _T>& mat)
{
	return mat.elem[0][0] * (mat.elem[1][1] * mat.elem[2][2] * mat.elem[3][3] - mat.elem[3][1] * mat.elem[2][2] * mat.elem[1][3])
			+ mat.elem[0][1] * (mat.elem[1][2] * mat.elem[2][3] * mat.elem[3][0] - mat.elem[1][0] * mat.elem[2][3] * mat.elem[3][2])
			+ mat.elem[0][2] * (mat.elem[1][3] * mat.elem[2][0] * mat.elem[3][1] - mat.elem[1][1] * mat.elem[2][0] * mat.elem[3][3])
			+ mat.elem[0][3] * (mat.elem[1][0] * mat.elem[2][1] * mat.elem[3][2] - mat.elem[1][2] * mat.elem[2][1] * mat.elem[3][0]);
}

template<typename _T>
using matrix1x1 = matrix<1, 1, _T>;
template<typename _T>
using matrix2x2 = matrix<2, 2, _T>;
template<typename _T>
using matrix3x3 = matrix<3, 3, _T>;
template<typename _T>
using matrix4x4 = matrix<4, 4, _T>;
template<typename _T>
using matrix5x5 = matrix<5, 5, _T>;
template<typename _T>
using matrix6x6 = matrix<6, 6, _T>;
template<typename _T>
using matrix7x7 = matrix<7, 7, _T>;
template<typename _T>
using matrix8x8 = matrix<8, 8, _T>;
template<typename _T>
using matrix9x9 = matrix<9, 9, _T>;

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
#endif//_MATH_MATRIX
