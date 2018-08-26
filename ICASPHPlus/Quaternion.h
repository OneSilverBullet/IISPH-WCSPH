#pragma once
#include <cmath>

class Quaternion
{
public:
	Quaternion();
	Quaternion(double _x, double _y, double _z, double _r = 0.0);
	Quaternion(const double* m);
	~Quaternion();

public:
	void setVertex(double _x, double _y, double _z, double _r = 0.0);
	// 乘法
	Quaternion dot(Quaternion left, Quaternion right);
	
	/* 计算旋转后的四元数。旋转轴正方向，顺时针。
	*	radian	 ------ 旋转角(弧度)
	*	_x,_y,_z ------ 旋转轴
	*/
	Quaternion rotational(double radian, double _x, double _y, double _z);
	
	/* 获取旋转四元数
	*	radian	   ------ 旋转角(弧度)
	*	_x, _y, _z ------ 旋转轴
	*/
	Quaternion getRotationalQuaternion(double radian, double _x, double _y, double _z);

public:
	Quaternion operator * (const Quaternion right);
	Quaternion normalize();

	float* GetMatrixLH();

public:
	double w; // 实部 
	double x;
	double y;
	double z;
};

