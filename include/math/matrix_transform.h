#ifndef _MATH_TRANSFORM
#define _MATH_TRANSFORM

#include <math/vector.h>
#include <math/matrix.h>
#include <math/quaternion.h>

//采用右手坐标系
namespace math
{
//缩放矩阵
template<size_t _Order, typename _T = float>
inline matrix<_Order, _Order, _T> scale(_T scale_factor)
{
	return math::identity<matrix<_Order, _T> >() * scale_factor;
}

//平移矩阵，矩阵维度比矢量大1
template<size_t _Order, typename _T = float>
inline matrix<_Order + 1, _Order + 1, _T> shift(vector<_Order, _T>& shift_v)
{
	matrix<_Order + 1, _Order + 1, _T> shift_mat = math::identity<matrix<_Order + 1, _Order + 1, _T> >();
	for(int i = 0; i <= _Order; ++i)
		shift_mat[i][_Order] = shift_v[i];
	return shift_mat;
}

template<typename _T = float>
inline matrix4x4<_T> shift(vector3<_T> shift_v)
{
	_T m_elem[] =
			{
					1, 0, 0, shift_v[0],
					0, 1, 0, shift_v[1],
					0, 0, 1, shift_v[2],
					0, 0, 0, 1
			};
	return tplmp::cast<matrix4x4<_T> >(m_elem);
}

template<typename _T = float>
inline matrix4x4<_T> Rx(_T rad)
{
	_T m_elem[] =
			{
					1, 0, 0, 0,
					0, math::cos(rad), -math::sin(rad), 0,
					0, math::sin(rad), math::cos(rad), 0,
					0, 0, 0, 1
			};
	return tplmp::cast<matrix4x4<_T> >(m_elem);
}

template<typename _T = float>
inline matrix4x4<_T> Ry(_T rad)
{
	_T m_elem[] =
			{
					math::cos(rad), 0, math::sin(rad), 0,
					0, 1, 0, 0,
					-math::sin(rad), 0, math::cos(rad), 0,
					0, 0, 0, 1
			};
	return tplmp::cast<matrix4x4<_T> >(m_elem);
}

template<typename _T = float>
inline matrix4x4<_T> Rz(_T rad)
{
	_T m_elem[] =
			{
					math::cos(rad), -math::sin(rad), 0, 0,
					math::sin(rad), math::cos(rad), 0, 0,
					0, 0, 1, 0,
					0, 0, 0, 1
			};
	return tplmp::cast<matrix4x4<_T> >(m_elem);
}

//外在旋转
template<typename _T = float>
inline matrix4x4<_T> Rxyz(_T rad_x, _T rad_y, _T rad_z)
{
	return Rx(rad_z) * Ry(rad_y) * Rz(rad_x);
}

//内在旋转
template<typename _T = float>
inline matrix4x4<_T> Rzyx(_T pitch, _T yaw, _T roll)
{
	return Rz(roll) * Ry(yaw) * Rx(pitch);
}

//绕轴转动，四元数转换为旋转矩阵
template<typename _T = float>
inline matrix4x4<_T> R(vector3<_T> axis, _T rad)
{
	return quaternion<_T>::R(axis, rad);
}

//投影矩阵。OpenGL的NDC为左手坐标系，世界坐标系为右手坐标系，near和far应该取负数，向z轴负方向看去
template<typename _T = float>
inline matrix4x4<_T> frustum(_T z_near, _T z_far, _T fov, _T left, _T right, _T top, _T bottom)
{
	_T m_elem[] = {
			2 * z_near / (right - left), 0, (right + left) / (right - left), 0,
			0, 2 * z_near / (top - bottom), (top + bottom) / (top - bottom), 0,
			0, 0, -(z_far + z_near) / (z_far - z_near), -2 * z_near * z_far / (z_far - z_near),
			0, 0, -1, 0
	};
	return tplmp::cast<matrix4x4<_T> >(m_elem);
}

//投影矩阵。OpenGL的NDC为左手坐标系，世界坐标系为右手坐标系，z_near和z_far取正
template<typename _T = float>
inline matrix4x4<_T> frustum(_T z_near, _T z_far, _T fov, _T aspect_ratio)
{
	_T m_elem[] = {
			1 / (tan(fov / 2) * aspect_ratio), 0, 0, 0,
			0, 1 / tan(fov / 2), 0, 0,
			0, 0, -(z_far + z_near) / (z_far - z_near), -2 * z_near * z_far / (z_far - z_near),
			0, 0, -1, 0
	};
	return tplmp::cast<matrix4x4<_T> >(m_elem);
}

}

#endif//_MATH_TRANSFORM
