#pragma once
#include <vector>
#include <math.h>
#include "../ICASPHPlus/SPH/Vector3f.h"
#include <algorithm>

//////////////////////////////////////////////////////////////////////////
//Matrix3f类
//作用：用于刚体模拟当中的旋转等数据，在SPH计算当中，没有使用该类。
//////////////////////////////////////////////////////////////////////////
class Matrix3f
{
public:
	double value[3][3];

public:
	Matrix3f(){
		for (int i=0; i<3; i++)
		{
			for (int j=0; j<3; j++)
			{
				value[i][j] = 0.0f;
			}
		}
	}
	Matrix3f(Matrix3f& v)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				value[i][j] = v[i][j];
			}
		}
	}

	Vector3f& operator*(Vector3f& a);
	Vector3f& operator*(const Vector3f& a);
	Matrix3f& operator*(Matrix3f& a);
	Matrix3f& operator*(const Matrix3f& a);


	double* operator[](int x);
	const double* operator[](const int x)const;
	//重载赋值操作符
	Matrix3f& operator=(Matrix3f& a);
	Matrix3f& operator=(const Matrix3f& a);

	//关于旋转矩阵的操作
	//形参为角速度的x轴弧度
	Matrix3f RotateX(float omigaX);
	Matrix3f RotateY(float omigaY);
	Matrix3f RotateZ(float omigaZ);

	Matrix3f Matrix3Multiplication(Matrix3f m1, Matrix3f m2);
	Matrix3f Rotate(float omigaX, float omigaY, float omigaZ);

	//矩阵的求逆
	//计算逆矩阵
	Matrix3f Inverse();

	//根据Vec3返回对应的对角矩阵
	Matrix3f Diagonal(Vector3f& a);
	//返回当前矩阵的转置矩阵
	Matrix3f Tranpose();

	friend Vector3f& operator*(const Matrix3f& m, Vector3f& a);

};