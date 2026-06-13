#ifndef _MATH_INTEGRAL
#define _MATH_INTEGRAL

#include <tplmp/base.h>

#include <math/elementary_function.h>
#include <math/vector.h>
#include <math/matrix.h>

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
 * @brief _Order阶龙格库塔法使用的Butcher Tableau，即算法参数表
 */
template<size_t _Order, typename _T>
struct butcher_table
{
	vector<_Order, _T> c;
	matrix<_Order, _Order, _T> A;
	vector<_Order, _T> b;

	butcher_table() = default;

	butcher_table(const vector<_Order, _T>& _c, const matrix<_Order, _Order, _T>& _A, const vector<_Order, _T>& _b) :
			c(_c), A(_A), b(_b)
	{
	}
};

/**
 * 欧拉法
 */
template<typename _T>
const butcher_table<1, _T>& rk1_euler()
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
const butcher_table<2, _T>& rk2_midpoint()
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
const butcher_table<2, _T>& rk2_heun()
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
const butcher_table<2, _T>& rk2_ralston()
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
const butcher_table<3, _T>& rk3_classic()
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
const butcher_table<3, _T>& rk3_heun()
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
const butcher_table<3, _T>& rk3_ralston()
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
const butcher_table<4, _T>& rk4_classic()
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
const butcher_table<4, _T>& rk4_three_eighths()
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
const butcher_table<4, _T>& rk4_ralston()
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
const butcher_table<1, _T>& rk1_symplectic_gauss_legendre_1()
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
const butcher_table<2, _T>& rk2_symplectic_gauss_legendre_2()
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
const butcher_table<3, _T>& rk3_symplectic_gauss_legendre_3()
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
const butcher_table<3, _T>& rk3_symplectic_lobatto_iiia()
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
 * F(t, f(t))=f'(t)
 * 适用于直接使用t表达导数的情形。
 */
template<size_t _Dim, typename _T, typename _FuncType>
struct F_t_dependent
{
	static_assert(
			tplmp::type_equal<
			decltype(tplmp::decl<_FuncType>::val()(tplmp::decl<_T>::val())),
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
			decltype(tplmp::decl<_FuncType>::val()(tplmp::decl<vector<_Dim, _T> >::val())),
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
 * @brief _Order阶显式龙格库塔法单步积分。
 * 		  给定原函数在t0时刻的初始值，计算其在t0+step_size处的值。
 * @param F 原函数的导数数值计算函数F
 * @param t0 积分起点
 * @param f_t0 原函数f(t)在t0时刻的初始值f(t0)，即F在t0的积分值
 * @param dt 单步积分步长，步长越小积分精度越高。
 * @param table 算法参数表
 */
template<size_t _Dim, size_t _Order, typename _T, typename _Derivative>
vector<_Dim, _T> runge_kutta_step(_Derivative F, _T t0, const vector<_Dim, _T>& f_t0, _T dt, const butcher_table<_Order, _T>& table)
{
	vector<_Order, vector<_Dim, _T> > k;
	k[0] = F(t0, f_t0);
	/**								   	  i-1
	 * 计算各阶段斜率向量ki=F(t0+ci*dt, f(t0)+dt*Σ{A[i,j]*kj})
	 * 									  j=0
	 * 其中A、b、c为butcher_table中的参数值，即各ki斜率向量的权重。
	 */
	for(size_t i = 1; i < _Order; ++i)
	{
		vector<_Dim, _T> df = vector<_Dim, _T>::zero();
		for(size_t j = 0; j < i; ++j)
		{
			_T Aij = table.A[i][j];
			if(Aij != 0)
			{
				df += Aij * k[j];
			}
		}
		k[i] = F(t0 + table.c[i] * dt, f_t0 + df * dt);
	}
	/**				  _Order-1
	 * f(t0+dt)=f(t0)+dt*Σ{bi*ki}
	 * 					i=0
	 */
	vector<_Dim, _T> df = vector<_Dim, _T>::zero();
	for(size_t i = 0; i < _Order; ++i)
	{
		_T bi = table.b[i];
		if(bi != 0)
		{
			df += bi * k[i];
		}
	}
	return f_t0 + df * dt;
}

}

#endif//_MATH_INTEGRAL
