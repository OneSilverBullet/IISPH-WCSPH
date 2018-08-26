#include "Quaternion.h"


Quaternion::Quaternion()
{
	this->x = 0.0;
	this->y = 0.0;
	this->z = 0.0;
	this->w = 1.0;
}


Quaternion::~Quaternion()
{
}

Quaternion::Quaternion(double _x, double _y, double _z, double _r)
{
	this->x = _x;
	this->y = _y;
	this->z = _z;
	this->w = _r;
}

Quaternion Quaternion::normalize()
{
	double temp = sqrt(x*x + y*y + z*z + w*w);
	return Quaternion(x / temp, y / temp, z / temp, w / temp);
}

void Quaternion::setVertex(double _x, double _y, double _z, double _r)
{
	this->x = _x;
	this->y = _y;
	this->z = _z;
	this->w = _r;
}

Quaternion Quaternion::dot(Quaternion left, Quaternion right)
{
	Quaternion qua;
	double d1, d2, d3, d4;

	d1 = left.w * right.w;
	d2 = -left.x * right.x;
	d3 = -left.y * right.y;
	d4 = -left.z * right.z;
	qua.w = d1 + d2 + d3 + d4;

	d1 = left.w * right.x;
	d2 = right.w * left.x;
	d3 = left.y * right.z;
	d4 = -left.z * right.y;
	qua.x = d1 + d2 + d3 + d4;

	d1 = left.w * right.y;
	d2 = right.w * left.y;
	d3 = left.z * right.x;
	d4 = -left.x * right.z;
	qua.y = d1 + d2 + d3 + d4;

	d1 = left.w * right.z;
	d2 = right.w * left.z;
	d3 = left.x * right.y;
	d4 = -left.y * right.x;
	qua.z = d1 + d2 + d3 + d4;

	return qua;
}

Quaternion Quaternion::operator * (const Quaternion right)
{
	Quaternion qua;
	double d1, d2, d3, d4;

	d1 = this->w * right.w;
	d2 = -this->x * right.x;
	d3 = -this->y * right.y;
	d4 = -this->z * right.z;
	qua.w = d1 + d2 + d3 + d4;

	d1 = this->w * right.x;
	d2 = right.w * this->x;
	d3 = this->y * right.z;
	d4 = -this->z * right.y;
	qua.x = d1 + d2 + d3 + d4;

	d1 = this->w * right.y;
	d2 = right.w * this->y;
	d3 = this->z * right.x;
	d4 = -this->x * right.z;
	qua.y = d1 + d2 + d3 + d4;

	d1 = this->w * right.z;
	d2 = right.w * this->z;
	d3 = this->x * right.y;
	d4 = -this->y * right.x;
	qua.z = d1 + d2 + d3 + d4;

	return qua;
}

Quaternion Quaternion::getRotationalQuaternion(double radian, double _x, double _y, double _z)
{
	Quaternion qua;
	double temp;
	double cr, cv;

	qua.w = qua.x = qua.y = qua.z = 0.0;

	temp = _x *  _x + _y *  _y + _z *  _z;

	if (temp == 0.0)
		return qua;

	temp = 1.0 / sqrt(temp);
	_x *= temp;
	_y *= temp;
	_z *= temp;

	cr = cos(0.5 * radian);
	cv = sin(0.5 * radian);

	qua.w = cr;
	qua.x = cv * _x;
	qua.y = cv * _y;
	qua.z = cv * _z;

	return qua;
}

Quaternion Quaternion::rotational(double radian, double _x, double _y, double _z)
{
	if (radian == 0.0 || (_x == 0.0&&_y == 0.0&&_z == 0.0))
		return Quaternion(x, y, z, w);

	Quaternion p(x, y, z, w);
	Quaternion q = getRotationalQuaternion(radian, _x, _y, _z);
	Quaternion _q = getRotationalQuaternion(-radian, _x, _y, _z);

	p = q*p;
	p = p*_q;

	return p;
}

// 由旋转四元数计算旋转矩阵
float* Quaternion::GetMatrixLH()
{
	float matrix[16];
	float xx = x*x;
	float yy = y*y;
	float zz = z*z;
	float xy = x*y;
	float wz = w*z;
	float wy = w*y;
	float xz = x*z;
	float yz = y*z;
	float wx = w*x;

	matrix[0] = 1.0f - 2 * (yy + zz);
	matrix[1] = 2 * (xy - wz);
	matrix[2] = 2 * (wy + xz);
	matrix[3] = 0.0f;

	matrix[4] = 2 * (xy + wz);
	matrix[5] = 1.0f - 2 * (xx + zz);
	matrix[6] = 2 * (yz - wx);
	matrix[7] = 0.0f;

	matrix[8] = 2 * (xy - wy);
	matrix[9] = 2 * (yz + wx);
	matrix[10] = 1.0f - 2 * (xx + yy);
	matrix[11] = 0.0f;

	matrix[12] = 0.0f;
	matrix[13] = 0.0f;
	matrix[14] = 0.0f;
	matrix[15] = 1.0f;

	return matrix;
};

// 由旋转矩阵创建四元数
Quaternion::Quaternion(const double* m)
{
	double matrix[4][4];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix[i][j] = m[i * 4 + j];
		}
	}

	float tr, s, q[4];
	int i, j, k;

	int nxt[3] = { 1, 2, 0 };
	// 计算矩阵轨迹
	tr = matrix[0][0] + matrix[1][1] + matrix[2][2];

	// 检查矩阵轨迹是正还是负
	if (tr > 0.0f)
	{
		s = sqrt(tr + 1.0f);
		this->w = s / 2.0f;
		s = 0.5f / s;
		this->x = (matrix[1][2] - matrix[2][1]) * s;
		this->y = (matrix[2][0] - matrix[0][2]) * s;
		this->z = (matrix[0][1] - matrix[1][0]) * s;
	}
	else
	{
		// 轨迹是负
		// 寻找m11 m22 m33中的最大分量
		i = 0;
		if (matrix[1][1] > matrix[0][0]) i = 1;
		if (matrix[2][2] > matrix[i][i]) i = 2;
		j = nxt[i];
		k = nxt[j];

		s = sqrt((matrix[i][i] - (matrix[j][j] + matrix[k][k])) + 1.0f);
		q[i] = s * 0.5f;
		if (s != 0.0f) s = 0.5f / s;
		q[3] = (matrix[j][k] - matrix[k][j]) * s;
		q[j] = (matrix[i][j] - matrix[j][i]) * s;
		q[k] = (matrix[i][k] - matrix[k][i]) * s;
		this->x = q[0];
		this->y = q[1];
		this->z = q[2];
		this->w = q[3];
	}
};