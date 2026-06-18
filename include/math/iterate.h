#ifndef _MATH_ITERATE
#define _MATH_ITERATE

#include <tplmp/base.h>

#include <math/tensor.h>
#include <math/algebra.h>

/**
 * 数值迭代法
 */

namespace math
{
/**
 * @brief 牛顿迭代法求解向量方程组F(x)=0
 *		  F(x(i+1))≈F(x(i))+J(x(i))Δx，其中i为迭代次数，Δx=x(i+1)-x(i)，J=∂F/∂x，x(0)为猜测解，
 *		  假设本次迭代结果F(x(i+1))=0，则方程变为J(x(i))Δx=-F(x(i))，即可求出Δx。但实际上代入解出的Δx后，F(x(i+1))=F(x(i)+Δx)≤F(x(i))，并不一定是0。
 *		  但F(x(i))的值会随着迭代次数收敛，最终F(x(i))<eps时即可认为迭代到了满足误差范围的精确解。
 */
template<size_t _Dim, typename _T, typename _Equation, typename _Derivative>
bool newton_iterate(_Equation F, _Derivative J, vector<_T, _Dim>& x0, _T eps = precision<_T>::default_val(), size_t max_iter = 20, _T linear_eps = precision<_T>::default_val())
{
	bool solve_success = false;
	for(size_t iter = 0; iter < max_iter; ++iter)
	{
		vector<_T, _Dim> _Fxi = -F(x0);
		if(is_zero(_Fxi.norm(), eps))
			return true;
		vector<_T, _Dim> dx = solve_linear_system(J(x0), _Fxi, &solve_success, linear_eps);
		if(solve_success)
			x0 += dx;
		else
			return false;
	}
	return false;
}

}

#endif//_MATH_ITERATE
