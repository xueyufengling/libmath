#ifndef _MATH_QUATERNION
#define _MATH_QUATERNION

#include <tplmp/base.h>

#include <math/elementary_function.h>
#include <math/tensor.h>

namespace math
{
/**
 * 四元数
 */
template<typename _T = float>
class quaternion
{
public:
	_T w, x, y, z;

	quaternion(_T w, _T x, _T y, _T z)
	{
		this->w = w;
		this->x = x;
		this->y = y;
		this->z = z;
	}

	quaternion(_T w, vector3<_T> v)
	{
		this->w = w;
		this->x = v[0];
		this->y = v[1];
		this->z = v[2];
	}

	/**
	 * 共轭
	 */
	quaternion conj()
	{
		return quaternion(w, -x, -y, -z);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() + tplmp::decl<_T2>::val())>
	quaternion<_Result> operator+(quaternion<_T2> q)
	{
		return quaternion<_Result>(w + q.w, x + q.x, y + q.y, z + q.z);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() - tplmp::decl<_T2>::val())>
	quaternion<_Result> operator-(quaternion<_T2> q)
	{
		return quaternion<_Result>(w - q.w, x - q.x, y - q.y, z - q.z);
	}

	template<typename _T2, typename _Result = decltype(tplmp::decl<_T>::val() * tplmp::decl<_T2>::val())>
	quaternion<_Result> operator*(quaternion<_T2> q)
	{
		return quaternion<_Result>(w * q.w - x * q.x - y * q.y - z * q.z, y * q.z - q.y * z + w * q.x + q.w * x, z * q.x - q.z * x + w * q.y + q.w * y, x * q.y - q.x * y + w * q.z + q.w * z);
	}

	static vector3<_T> rotate(vector3<_T> v, vector3<_T> axis, _T rad)
	{
		auto p = quaternion(0, v);
		auto q = quaternion(math::cos(rad / 2), axis * math::sin(rad / 2));
		auto r = q * p * q.conj();
		return
		{	r.x,r.y,r.z};
	}

	//绕轴转动，四元数转换为旋转矩阵
	static inline matrix4x4<_T> R(vector3<_T> axis, _T rad)
	{
		_T m_elem[] = {math::cos(rad) + axis[0] * axis[0] * (1 - math::cos(rad)),
				axis[0] * axis[1] * (1 - math::cos(rad)) - axis[2] * math::sin(rad),
				axis[0] * axis[2] * (1 - math::cos(rad)) + axis[1] * math::sin(rad),
				0,
				axis[0] * axis[1] * (1 - math::cos(rad)) + axis[2] * math::sin(rad),
				math::cos(rad) + axis[1] * axis[1] * (1 - math::cos(rad)),
				axis[1] * axis[2] * (1 - math::cos(rad)) - axis[0] * math::sin(rad),
				0,
				axis[0] * axis[2] * (1 - math::cos(rad)) - axis[1] * math::sin(rad),
				axis[1] * axis[2] * (1 - math::cos(rad)) + axis[0] * math::sin(rad),
				math::cos(rad) + axis[2] * axis[2] * (1 - math::cos(rad)),
				0,
				0, 0, 0, 1};
		return tplmp::cast<matrix4x4<_T> >(m_elem);
	}
};

}

#endif//_MATH_QUATERNION
