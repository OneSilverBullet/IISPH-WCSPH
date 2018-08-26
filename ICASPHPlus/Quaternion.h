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
	// �˷�
	Quaternion dot(Quaternion left, Quaternion right);
	
	/* ������ת�����Ԫ������ת��������˳ʱ�롣
	*	radian	 ------ ��ת��(����)
	*	_x,_y,_z ------ ��ת��
	*/
	Quaternion rotational(double radian, double _x, double _y, double _z);
	
	/* ��ȡ��ת��Ԫ��
	*	radian	   ------ ��ת��(����)
	*	_x, _y, _z ------ ��ת��
	*/
	Quaternion getRotationalQuaternion(double radian, double _x, double _y, double _z);

public:
	Quaternion operator * (const Quaternion right);
	Quaternion normalize();

	float* GetMatrixLH();

public:
	double w; // ʵ�� 
	double x;
	double y;
	double z;
};

