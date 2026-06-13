#ifndef _MATH_CONSTANT
#define _MATH_CONSTANT

#include <quadmath.h>

/**
 * 数学常数定义
 */
namespace math
{
/**
 * @brief 高精度原始数值
 */
struct __constant_float128
{
	static constexpr __float128 pi = M_PIq;

	static constexpr __float128 _2pi = 2 * pi;

	static constexpr __float128 pi_2 = M_PI_2q;

	static constexpr __float128 pi_4 = M_PI_4q;

	static constexpr __float128 _1_pi = M_1_PIq;

	static constexpr __float128 e = M_Eq;

	static constexpr __float128 ln2 = M_LN2q;

	static constexpr __float128 ln10 = M_LN10q;
};

/**
 * @brief 根据模板精度转换
 */
template<typename _T>
struct constant
{
	static constexpr _T pi = (_T)__constant_float128::pi;

	static constexpr __float128 _2pi = (_T)__constant_float128::_2pi;

	static constexpr _T pi_2 = (_T)__constant_float128::pi_2;

	static constexpr _T pi_4 = (_T)__constant_float128::pi_4;

	static constexpr _T _1_pi = (_T)__constant_float128::_1_pi;

	static constexpr _T e = (_T)__constant_float128::e;

	static constexpr _T ln2 = (_T)__constant_float128::ln2;

	static constexpr _T ln10 = (_T)__constant_float128::ln10;
};
}

#endif//_MATH_CONSTANT
