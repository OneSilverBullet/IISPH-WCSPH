#include "Matrix3f.h"

Vector3f & Matrix3f::operator*(Vector3f & a)
{
	Vector3f result;
	result.SetValue(value[0][0] * a.x + value[0][1] * a.y + value[0][2] * a.z,
		value[1][0] * a.x + value[1][1] * a.y + value[1][2] * a.z,
		value[2][0] * a.x + value[2][1] * a.y + value[2][2] * a.z);
	return result;
}

Vector3f & Matrix3f::operator*(const Vector3f & a)
{
	Vector3f result;
	result.SetValue(value[0][0] * a.x + value[0][1] * a.y + value[0][2] * a.z,
		value[1][0] * a.x + value[1][1] * a.y + value[1][2] * a.z,
		value[2][0] * a.x + value[2][1] * a.y + value[2][2] * a.z);
	return result;
}

Matrix3f& Matrix3f::operator*(Matrix3f & a)
{
	Matrix3f result;
	int i, j;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result.value[i][j] = value[i][0] * a.value[0][j] +
				value[i][1] * a.value[1][j] +
				value[i][2] * a.value[2][j];
		}
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			value[i][j] = result.value[i][j];
		}
	}
	return *this;
}

Matrix3f & Matrix3f::operator*(const Matrix3f & a)
{
	Matrix3f result;
	int i, j;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			result.value[i][j] = value[i][0] * a.value[0][j] +
				value[i][1] * a.value[1][j] +
				value[i][2] * a.value[2][j];
		}
	}

	return result;
}

double * Matrix3f::operator[](int x)
{
	return value[x];
}

const double * Matrix3f::operator[](const int x) const
{
	return value[x];
}

//对数据进行拷贝
Matrix3f & Matrix3f::operator=(Matrix3f & a)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			value[i][j] = a[i][j];
		}
	}
	return *this;
}

Matrix3f & Matrix3f::operator=(const Matrix3f & a)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			value[i][j] = a[i][j];
		}
	}
	return *this;
}

Matrix3f Matrix3f::RotateX(float omigaX)
{
	Matrix3f result;
	result.value[0][0] = 1.0f;
	result.value[0][1] = 0.0f;
	result.value[0][2] = 0.0f;
	result.value[1][0] = 0.0f;
	result.value[1][1] = cos(omigaX);
	result.value[1][2] = -sin(omigaX);
	result.value[2][0] = 0.0f;
	result.value[2][1] = sin(omigaX);
	result.value[2][2] = cos(omigaX);
	return result;
}

Matrix3f Matrix3f::RotateY(float omigaY)
{
	Matrix3f result;
	result.value[0][0] = cos(omigaY);
	result.value[0][1] = 0.0f;
	result.value[0][2] = -sin(omigaY);
	result.value[1][0] = 0.0f;
	result.value[1][1] = 1.0f;
	result.value[1][2] = 0.0f;
	result.value[2][0] = sin(omigaY);
	result.value[2][1] = 0.0f;
	result.value[2][2] = cos(omigaY);
	return result;
}

Matrix3f Matrix3f::RotateZ(float omigaZ)
{
	Matrix3f result;
	result.value[0][0] = cos(omigaZ);
	result.value[0][1] = -sin(omigaZ);
	result.value[0][2] = 0.0f;
	result.value[1][0] = sin(omigaZ);
	result.value[1][1] = cos(omigaZ);
	result.value[1][2] = 0.0f;
	result.value[2][0] = 0.0f;
	result.value[2][1] = 0.0f;
	result.value[2][2] = 1.0f;
	return result;
}

Matrix3f Matrix3f::Matrix3Multiplication(Matrix3f m1, Matrix3f m2)
{
	Matrix3f result;
	int i, j;
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			result.value[i][j] = m1.value[i][0] * m2.value[0][j] +
				m1.value[i][1] * m2.value[1][j] +
				m1.value[i][2] * m2.value[2][j];
		}
	}
	return result;
}

Matrix3f Matrix3f::Rotate(float omigaX, float omigaY, float omigaZ)
{
	//依照x、y、z轴的顺序进行旋转
	Matrix3f result;
	result = Matrix3Multiplication(Matrix3Multiplication(RotateX(omigaX), RotateY(omigaY)), RotateZ(omigaZ));
	return result;
}

 
Matrix3f Matrix3f::Inverse()
{
	int i = 0;
	int j = 0;
	int k = 0;
	double data[3][6];
	//根据原函数构造初始求值函数
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 6; j++)
		{
			if (j < 3)
			{
				data[i][j] = value[i][j];
			}
			else if (j == i + 3)
				data[i][j] = 1.0;
			else
				data[i][j] = 0.0;
		}
	}

	int maxI = 0;
	for (i = 1; i < 3; i++)
	{
		if (fabs(data[maxI][0]) < fabs(data[i][0]))
			maxI = i;
	}
	if (maxI != 0)
	{
		double temp;
		for (j = 0; j < 6; j++)
		{
			temp = data[0][j];
			data[0][j] = data[maxI][j];
			data[maxI][j] = temp;
		}
	}
	double temp2;
	for (i = 0; i < 3; i++)
	{
		if (data[i][i] != 0)
			temp2 = 1.0 / data[i][i];
		else
		{
			//对无逆矩阵输出到后台上
			cout << "此矩阵无逆! " << endl;
			Matrix3f null;
			return null;
		}
		for (j = 0; j < 6; j++)
			data[i][j] *= temp2;
		for (j = 0; j < 3; j++)
		{
			if (j != i)
			{
				double temp3 = data[j][i];
				for (k = 0; k < 6; k++)
					data[j][k] -= temp3*data[i][k];
			}
		}
	}
	//求逆完毕，提取必要信息
	Matrix3f result;
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 6; j++)
		{
			if (j >= 3)
			{
				result.value[i][j] = data[i][j];
			}
		}
	}
	return result;
}

Matrix3f Matrix3f::Diagonal(Vector3f & a)
{
	Matrix3f res;
	res.value[0][0] = a.x;
	res.value[1][1] = a.y;
	res.value[2][2] = a.z;
	return res;
}

Matrix3f Matrix3f::Tranpose()
{
	Matrix3f res;
	for (int i=0;i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			res[i][j] = value[j][i];
		}
	}
	return res;
}

Vector3f & operator*(const Matrix3f& m, Vector3f & a)
{
	Vector3f result;
	result.SetValue(m.value[0][0] * a.x + m.value[0][1] * a.y + m.value[0][2] * a.z,
		m.value[1][0] * a.x + m.value[1][1] * a.y + m.value[1][2] * a.z,
		m.value[2][0] * a.x + m.value[2][1] * a.y + m.value[2][2] * a.z);
	return result;
}
