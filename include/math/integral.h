#ifndef _MATH_INTEGRAL
#define _MATH_INTEGRAL

#include <tplmp/base.h>

#include <sstream>

#include <math/elementary_function.h>
#include <math/vector.h>
#include <math/matrix.h>
#include <math/algebra.h>
#include <math/auto_diff.h>

namespace math
{
/**
 * ODE求解，目标函数只允许有一个自变量。
 * 积分原理：给定初始点t0和初始值f(t0)，以及f(t)的导数计算函数F(t, f(t))，且满足F(t, f(t))=f'(t)。
 * f'(t)的表达式必须已知。根据情况，F只需要t或f(t)中的一个。如直接令F为使用t表达的导数代数式，则直接代入t的值计算f'(t)即可，无需使用f(t)；若使用自动微分而非直接使用含t的表达式，则根据点的值f(t)去计算导数f'(t)，无需使用t，
 * 给定值t0、导数f(t0)后计算导数F(t0, f(t0))=f'(t0)，单步积分根据初始点值和导数计算t0+Δt处的值f(t0+Δt)，再将(t0+Δt, f(t0+Δt))代入F计算，得到导数f'(t0+Δt)。
 * 在单步积分中，根据(t0, f(t0))计算出了(t0+Δt, f(t0+Δt))。重复这个单步积分步骤，得到一系列连续的点，即是原函数在t0, t0+Δt, t0+2Δt, t0+3Δt...的离散值，如果步长足够小，可以看作这一系列离散积分值是原函数f(t)。
 */

/**
 * F(t, f(t))=f'(t)
 * 适用于直接使用t表达导数的情形。
 */
template<size_t _Dim, typename _T, typename _FuncType>
struct F_t_dependent
{
	static_assert(
			tplmp::type_equal<
			typename tplmp::callable_ret<_FuncType, _T>::type,
			vector<_Dim, _T>
			>::value,
			"derivative expr must be t-dependent");

	_FuncType F;

	inline F_t_dependent(_FuncType&& direct_expr) :
			F(tplmp::forward<_FuncType>(direct_expr))
	{
	}

	inline vector<_Dim, _T> operator()(_T t, const vector<_Dim, _T>& f_t) const
	{
		return F(t);
	}
};

/**
 * F(t, f(t))=auto_diff(f(t))
 * 适用于自动微分，直接从值f(t)计算。
 */
template<size_t _Dim, typename _T, typename _FuncType>
struct F_ft_dependent
{
	static_assert(
			tplmp::type_equal<
			typename tplmp::callable_ret<_FuncType, vector<_Dim, _T> >::type,
			vector<_Dim, _T>
			>::value,
			"derivative expr must be f(t)-dependent");

	_FuncType F;

	inline F_ft_dependent(_FuncType&& ad_expr) :
			F(tplmp::forward<_FuncType>(ad_expr))
	{
	}

	inline vector<_Dim, _T> operator()(_T t, const vector<_Dim, _T>& f_t) const
	{
		return F(f_t);
	}
};

/**
 * 龙格库塔法
 */
namespace runge_kutta
{
/**
 * @brief _Stage阶龙格库塔法使用的Butcher Tableau，即算法参数表
 * 		  若A为严格的下三角矩阵（对角线及上三角均为0），则属于显式RK法；若A含有对角非零元或上三角非零元，则属于隐式RK法。
 * 		  对保辛系统积分时需要使用隐式RK法，否则系统的哈密顿量会剧烈单调漂移发散。
 */
template<size_t _Stage, typename _T>
struct butcher_table
{
	vector<_Stage, _T> c;
	matrix<_Stage, _Stage, _T> A;
	vector<_Stage, _T> b;

	butcher_table() = default;

	butcher_table(const vector<_Stage, _T>& _c, const matrix<_Stage, _Stage, _T>& _A, const vector<_Stage, _T>& _b) :
			c(_c), A(_A), b(_b)
	{
	}

	/**
	 * @brief 是否是显式RK法的参数表
	 */
	inline bool is_explicit(_T eps = 0) const
	{
		//判断对角线及上三角有无非零元，如果有则是隐式RK法
		for(size_t i = 0; i < _Stage; ++i)
		{
			for(size_t j = i; j < _Stage; ++j)
			{
				if(!is_zero(A[i][j], eps))
					return false;
			}
		}
		return true;
	}

	/**
	 * @brief 是否是辛的
	 */
	bool is_symplectic(_T eps = 0) const
	{
		//bi*Aij+bj*Aji=bi*bj, ∀i,j=1,…,_Stage
		for(int i = 0; i < _Stage; i++)
		{
			for(int j = 0; j < _Stage; j++)
			{
				if(is_zero(b[i] * A[i][j] + b[j] * A[j][i] - b[i] * b[j], eps))
					return false;
			}
		}
		return true;
	}

	inline operator std::string() const
	{
		std::ostringstream oss;
		oss << "c:\n" << c << '\n';
		oss << "A:\n" << A << '\n';
		oss << "b:\n" << b;
		return oss.str();
	}
};

template<size_t _Stage, typename _T>
inline std::ostream& operator<<(std::ostream& os, const butcher_table<_Stage, _T>& table)
{
	os << ((std::string)table);
	return os;
}

/**
 * 欧拉法
 */
template<typename _T>
const butcher_table<1, _T>& erk1_euler()
{
	static const butcher_table<1, _T> table(
			{_T(0)},
			{_T(0)},
			{_T(1)}
			);
	return table;
}

/**
 * 中点法
 */
template<typename _T>
const butcher_table<2, _T>& erk2_midpoint()
{
	static const butcher_table<2, _T> table(
			{_T(0), _T(1) / _T(2)},
			{_T(0), _T(0), _T(1) / _T(2), _T(0)},
			{_T(0), _T(1)});
	return table;
}

/**
 * 改进欧拉法
 */
template<typename _T>
const butcher_table<2, _T>& erk2_heun()
{
	static const butcher_table<2, _T> table(
			{_T(0), _T(1)},
			{_T(0), _T(0), _T(1), _T(0)},
			{_T(1) / _T(2), _T(1) / _T(2)}
			);
	return table;
}

/**
 * Ralston法
 */
template<typename _T>
const butcher_table<2, _T>& erk2_ralston()
{
	static const butcher_table<2, _T> table(
			{_T(0), _T(2) / _T(3)},
			{_T(0), _T(0), _T(2) / _T(3), _T(0)},
			{_T(1) / _T(4), _T(3) / _T(4)}
			);
	return table;
}

/**
 * 经典RK3法
 */
template<typename _T>
const butcher_table<3, _T>& erk3_classic()
{
	static const butcher_table<3, _T> table(
			{_T(0), _T(1) / _T(2), _T(1)},
			{_T(0), _T(0), _T(0),
					_T(1) / _T(2), _T(0), _T(0),
					_T(-1), _T(2), _T(0)},
			{_T(1) / _T(6), _T(2) / _T(3), _T(1) / _T(6)}
			);
	return table;
}

/**
 * Heun三阶法
 */
template<typename _T>
const butcher_table<3, _T>& erk3_heun()
{
	static const butcher_table<3, _T> table(
			{_T(0), _T(1) / _T(3), _T(2) / _T(3)},
			{_T(0), _T(0), _T(0),
					_T(1) / _T(3), _T(0), _T(0),
					_T(0), _T(2) / _T(3), _T(0)},
			{_T(1) / _T(4), _T(0), _T(3) / _T(4)}
			);
	return table;
}

/**
 * Ralston三阶法
 */
template<typename _T>
const butcher_table<3, _T>& erk3_ralston()
{
	static const butcher_table<3, _T> table(
			{_T(0), _T(1) / _T(2), _T(3) / _T(4)},
			{_T(0), _T(0), _T(0),
					_T(1) / _T(2), _T(0), _T(0),
					_T(0), _T(3) / _T(4), _T(0)},
			{_T(2) / _T(9), _T(1) / _T(3), _T(4) / _T(9)}
			);
	return table;
}

/**
 * 经典RK4法
 */
template<typename _T>
const butcher_table<4, _T>& erk4_classic()
{
	static const butcher_table<4, _T> table(
			{_T(0), _T(1) / _T(2), _T(1) / _T(2), _T(1)},
			{_T(0), _T(0), _T(0), _T(0),
					_T(1) / _T(2), _T(0), _T(0), _T(0),
					_T(0), _T(1) / _T(2), _T(0), _T(0),
					_T(0), _T(0), _T(1), _T(0)},
			{_T(1) / _T(6), _T(1) / _T(3), _T(1) / _T(3), _T(1) / _T(6)}
			);
	return table;
}

/**
 * 3/8法则法
 */
template<typename _T>
const butcher_table<4, _T>& erk4_three_eighths()
{
	static const butcher_table<4, _T> table(
			{_T(0), _T(1) / _T(3), _T(2) / _T(3), _T(1)},
			{_T(0), _T(0), _T(0), _T(0),
					_T(1) / _T(3), _T(0), _T(0), _T(0),
					_T(-1) / _T(3), _T(1), _T(0), _T(0),
					_T(1), -_T(1), _T(1), _T(0)},
			{_T(1) / _T(8), _T(3) / _T(8), _T(3) / _T(8), _T(1) / _T(8)}
			);
	return table;
}

/**
 * Ralston四阶法
 */
template<typename _T>
const butcher_table<4, _T>& erk4_ralston()
{
	static const butcher_table<4, _T> table(
			{_T(0), _T(4) / _T(10), _T(6) / _T(10), _T(1)},
			{_T(0), _T(0), _T(0), _T(0),
					_T(4) / _T(10), _T(0), _T(0), _T(0),
					_T(258) / _T(1000), _T(-442) / _T(1000), _T(0), _T(0),
					_T(3) / _T(10), _T(6) / _T(10), _T(1) / _T(10), _T(0)},
			{_T(17) / _T(81), _T(46) / _T(81), _T(18) / _T(81), _T(0)}
			);
	return table;
}

/**
 * RK辛积分法使用参数，用于哈密顿系统
 */

/**
 * Gauss-Legendre 1级
 * 2阶辛方法
 */
template<typename _T>
const butcher_table<1, _T>& irk1_symplectic_gauss_legendre_1()
{
	static const butcher_table<1, _T> table(
			{_T(1) / _T(2)},
			{_T(1) / _T(2)},
			{_T(1)}
			);
	return table;
}

/**
 * Gauss-Legendre 2级
 * 4阶辛方法
 */
template<typename _T>
const butcher_table<2, _T>& irk2_symplectic_gauss_legendre_2()
{
	static const _T sqrt3 = sqrt(_T(3));
	static const butcher_table<2, _T> table(
			{_T(1) / _T(2) - sqrt3 / _T(6),
					_T(1) / _T(2) + sqrt3 / _T(6)},
			{_T(1) / _T(4), _T(1) / _T(4) - sqrt3 / _T(6),
					_T(1) / _T(4) + sqrt3 / _T(6), _T(1) / _T(4)},
			{_T(1) / _T(2), _T(1) / _T(2)}
			);
	return table;
}

/**
 * Gauss-Legendre 3级
 * 6阶辛方法
 */
template<typename _T>
const butcher_table<3, _T>& irk3_symplectic_gauss_legendre_3()
{
	static const _T sqrt15 = sqrt(_T(15));
	static const butcher_table<3, _T> table(
			{_T(1) / _T(2) - sqrt15 / _T(10),
					_T(1) / _T(2),
					_T(1) / _T(2) + sqrt15 / _T(10)},
			{_T(5) / _T(36), _T(2) / _T(9) - sqrt15 / _T(15), _T(5) / _T(36) - sqrt15 / _T(30),
					_T(5) / _T(36) + sqrt15 / _T(24), _T(2) / _T(9), _T(5) / _T(36) - sqrt15 / _T(24),
					_T(5) / _T(36) + sqrt15 / _T(30), _T(2) / _T(9) + sqrt15 / _T(15), _T(5) / _T(36)},
			{_T(5) / _T(18), _T(4) / _T(9), _T(5) / _T(18)}
			);
	return table;
}

/**
 * Lobatto IIIA
 * 4阶辛方法
 */
template<typename _T>
const butcher_table<3, _T>& irk3_symplectic_lobatto_iiia()
{
	static const butcher_table<3, _T> table(
			{_T(0),
					_T(1) / _T(2),
					_T(1)},
			{_T(0), _T(0), _T(0),
					_T(1) / _T(4), _T(1) / _T(4), _T(0),
					_T(0), _T(1), _T(0)},
			{_T(1) / _T(6), _T(2) / _T(3), _T(1) / _T(6)}
			);
	return table;
}

/**
 * @brief 根据传入的迭代过程中ki数组计算预测的原函数值。
 * @param i RK步骤序号，即stage，从0 -> _Stage-1
 */
template<size_t _Dim, size_t _Stage, typename _T>
__attribute__((always_inline)) inline vector<_Dim, _T> predicted_val(size_t i, const vector<_Dim, _T>& f_t0, _T dt, const vector<_Stage, vector<_Dim, _T> >& k, const butcher_table<_Stage, _T>& table, _T eps = 0)
{
	/**		   i-1
	 * f(t0)+dt*Σ{A[i,j]*kj})
	 * 	       j=0
	 * 其中A、b、c为butcher_table中的参数值，即各ki斜率向量的权重。
	 */
	vector<_Dim, _T> df = vector<_Dim, _T>::zero();
	for(size_t j = 0; j < i; ++j)
	{
		_T Aij = table.A[i][j];
		if(!is_zero(Aij, eps))
		{
			df += Aij * k[j];
		}
	}
	return f_t0 + df * dt;
}

/**
 * @brief 根据计算出满足精度要求的ki数组计算最终的原函数值。
 */
template<size_t _Dim, size_t _Stage, typename _T>
__attribute__((always_inline)) inline vector<_Dim, _T> final_val(const vector<_Dim, _T>& f_t0, _T dt, const vector<_Stage, vector<_Dim, _T> >& k, const butcher_table<_Stage, _T>& table, _T eps = 0)
{
	/**				  _Stage-1
	 * f(t0+dt)=f(t0)+dt*Σ{bi*ki}
	 * 					i=0
	 */
	vector<_Dim, _T> df = vector<_Dim, _T>::zero();
	for(size_t i = 0; i < _Stage; ++i)
	{
		_T bi = table.b[i];
		if(!is_zero(bi, eps))
		{
			df += bi * k[i];
		}
	}
	return f_t0 + df * dt;
}

/**
 * @brief _Stage阶显式龙格库塔法单步积分。F的调用次数为_Stage。
 * 		  给定原函数在t0时刻的初始值，计算其在t0+step_size处的值。
 * @param F 原函数的导数数值计算函数F
 * @param t0 积分起点
 * @param f_t0 原函数f(t)在t0时刻的初始值f(t0)，即F在t0的积分值
 * @param dt 单步积分步长，步长越小积分精度越高。
 * @param table 算法参数表
 */
template<size_t _Dim, size_t _Stage, typename _T, typename _Derivative>
vector<_Dim, _T> explicit_step(_Derivative F, _T t0, const vector<_Dim, _T>& f_t0, _T dt, const butcher_table<_Stage, _T>& table, _T eps = 0)
{
	vector<_Stage, vector<_Dim, _T> > k;
	k[0] = F(t0, f_t0);
	/**
	 * 计算各阶段斜率向量ki=F(t0+ci*dt, predicted_val(i)})
	 * 其中A、b、c为butcher_table中的参数值，即各ki斜率向量的权重。
	 */
	for(size_t i = 1; i < _Stage; ++i)
	{
		k[i] = F(t0 + table.c[i] * dt, predicted_val(i, f_t0, dt, k, table, eps));
	}
	return final_val(f_t0, dt, k, table, eps);
}

template<size_t _Dim, size_t _Stage, typename _T, typename _Derivative>
vector<_Dim, _T> implicit_newton_step(_Derivative F, _T t0, const vector<_Dim, _T>& f_t0, _T dt, const butcher_table<_Stage, _T>& table, _T eps = 0, size_t max_iter = 20)
{
	vector<_Stage, vector<_Dim, _T> > k;
	k[0] = F(t0, f_t0);
	//设置k数组Newton迭代初始值全部为k0
	for(size_t i = 1; i < _Stage; ++i)
	{
		k[i] = k[0];
	}
	//开始Newton迭代法
	vector<_Stage, vector<_Dim, _T> > k_residual; //ki下一步预测值与当前值之间的差值
	//@formatter:off
	/**以_Stage=3，_Dim=2为例，构建_Stage*_Dim的6x6矩阵
	 	∂k₁(s+1)/∂k₁(s)  ∂k₂(s+1)/∂k₂(s)  ∂k₃(s+1)/∂k₃(s)   ← 其中s为stage数，即下一级ki关于本级ki各分量的偏导数矩阵Ji，是_Dim x _Dim的方阵
	  ┌────────────────┬────────────────┬────────────────┐
	  │ I(_Dim)-dtA₁₁J₁│     -dtA₁₂J₁   │     -dtA₁₃J₁   │  ← 第1级，I(_Dim)为_Dim阶恒等矩阵，Aij为Butcher表中的A矩阵参数
	  ├────────────────┼────────────────┼────────────────┤
	  │     -dtA₂₁J₂   │ I(_Dim)-dtA₂₂J₂│     -dtA₂₃J₂   │  ← 第2级
	  ├────────────────┼────────────────┼────────────────┤
	  │     -dtA₃₁J₃   │     -dtA₃₂J₃   │ I(_Dim)-dtA₃₃J₃│  ← 第3级，直到第_Stage级，整个矩阵是_Stage x _Stage块，每一块是_Dim x _Dim的Jacobi矩阵
	  └────────────────┴────────────────┴────────────────┘
	 *///@formatter:on
	matrix<_Stage, _Stage, matrix<_Dim, _Dim, _T> > k_JG;	//隐式RK法的全局Jacobi矩阵
	for(size_t iter = 0; iter < max_iter; ++iter)
	{
		//单次迭代过程中求解所有_Stage级的ki残差
		for(size_t i = 1; i < _Stage; ++i)
		{
			// 导数Fi(t, fi)从数学上讲实际上只依赖于t，但这里需要知道Fi相对于fi各分量的偏导数
			//Fi.value是计算出的ki估值，Fi.derivative就是ki关于ki各分量的偏导数
			ad_point<_Dim, _Dim, _T> Fi = F(t0 + table.c[i] * dt, predicted_val(i, f_t0, dt, k, table, eps)); //计算fi预测值及偏导数
			k_residual[i] = k[i] - Fi.value;
			for(size_t j = 0; j < _Stage; ++j)
			{
				if(i == j)
				{
					k_JG[i][j] = identity<matrix<_Dim, _Dim, _T> >() - (dt * table.A[i][j]) * Fi.derivative;
				}
				else
				{
					k_JG[i][j] = -(dt * table.A[i][j]) * Fi.derivative;
				}
			}
		}
		//求解完毕后所有残差同时联立求解线性方程组
		//k_residual、k_derivative全部纵向拼接
		vector<_Dim * _Stage, _T>& KR = (vector<_Dim * _Stage, _T>&)k_residual;
		matrix<_Dim * _Stage, _Dim * _Stage, _T>& JG = (matrix<_Dim * _Stage, _Dim * _Stage, _T>&)k_JG;		//TODO 矩阵拼接转换
		vector<_Dim * _Stage, _T>& dK;
		solve_linear_system(JG, dK, KR);
	}
	return calc_integral(f_t0, dt, k, table, eps);
}

}

}

#endif//_MATH_INTEGRAL
