#pragma once
#include <vector>
#include <math.h>
#include "../ICASPHPlus/SPH/Vector3f.h"
#include <algorithm>

//////////////////////////////////////////////////////////////////////////
//Matrix3f��
//���ã����ڸ���ģ�⵱�е���ת�����ݣ���SPH���㵱�У�û��ʹ�ø��ࡣ
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
	//���ظ�ֵ������
	Matrix3f& operator=(Matrix3f& a);
	Matrix3f& operator=(const Matrix3f& a);

	//������ת����Ĳ���
	//�β�Ϊ���ٶȵ�x�ỡ��
	Matrix3f RotateX(float omigaX);
	Matrix3f RotateY(float omigaY);
	Matrix3f RotateZ(float omigaZ);

	Matrix3f Matrix3Multiplication(Matrix3f m1, Matrix3f m2);
	Matrix3f Rotate(float omigaX, float omigaY, float omigaZ);

	//���������
	//���������
	Matrix3f Inverse();

	//����Vec3���ض�Ӧ�ĶԽǾ���
	Matrix3f Diagonal(Vector3f& a);
	//���ص�ǰ�����ת�þ���
	Matrix3f Tranpose();

	friend Vector3f& operator*(const Matrix3f& m, Vector3f& a);

};