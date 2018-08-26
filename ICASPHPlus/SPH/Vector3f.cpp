#include "Vector3f.h"

Vector3f & Vector3f::operator=(Vector3f* value)
{
	x = value->x;
	y = value->y;
	z = value->z;
	return *this;
}

double  Vector3f::operator*(Vector3f & value)
{
	double xTemp, yTemp, zTemp;
	xTemp = x*value.x;
	yTemp = y*value.y;
	zTemp = z*value.z;
	double result = xTemp + yTemp + zTemp;
	return result;
}

Vector3f Vector3f::operator*(double value)
{
	Vector3f result(0, 0, 0);
	result.x = x * value;
	result.y = y * value;
	result.z = z * value;
	return result;
}

Vector3f Vector3f::operator*(double value) const
{
	Vector3f result(0, 0, 0);
	result.x = x * value;
	result.y = y * value;
	result.z = z * value;
	return result;
}

Vector3f Vector3f::operator+(Vector3f & value)
{
	Vector3f result(0, 0, 0);
	result.x = x + value.x;
	result.y = y + value.y;
	result.z = z + value.z;
	return result;
}

Vector3f Vector3f::operator+(const Vector3f & value)
{
	Vector3f result(0, 0, 0);
	result.x = x + value.x;
	result.y = y + value.y;
	result.z = z + value.z;
	return result;
}

Vector3f Vector3f::operator-(Vector3f & value)
{
	Vector3f result(0, 0, 0);
	result.x = x - value.x;
	result.y = y - value.y;
	result.z = z - value.z;
	return result;
}

Vector3f& Vector3f::operator-()
{
	//返回倒置的vector
	this->x = -this->x;
	this->y = -this->y;
	this->z = -this->z;
	return *this;
}

Vector3f Vector3f::operator/(double value)
{
	Vector3f result(0, 0, 0);
	result.x = x / value;
	result.y = y / value;
	result.z = z / value;
	return result;
}

Vector3f & Vector3f::operator-=(Vector3f & a)
{
	x -= a.x;
	y -= a.y;
	z -= a.z;
	return *this;
}
Vector3f& Vector3f::operator-=(const Vector3f& a)
{
	x -= a.x;
	y -= a.y;
	z -= a.z;
	return *this;
}

Vector3f & Vector3f::operator+=(Vector3f & a)
{
	x += a.x;
	y += a.y;
	z += a.z;
	return *this;
}

Vector3f & Vector3f::operator+=(const Vector3f & a)
{
	x += a.x;
	y += a.y;
	z += a.z;
	return *this;
}

Vector3f & Vector3f::operator/=(int& a)
{
	x /= a;
	y /= a;
	z /= a;
	return *this;
}

Vector3f & Vector3f::operator/=(const int&  a)
{
	x /= a;
	y /= a;
	z /= a;
	return *this;
}

double Vector3f::operator[](const int a)
{
	switch (a)
	{
	case 0:
		return x;
		break;
	case 1:
		return y;
		break;
	case 2:
		return z;
		break;
	}
}

bool Vector3f::operator==(Vector3f & value)
{
	return (x == value.x&&y == value.y&&z == value.z);
}

double Vector3f::SquaredNorm() 
{
	return double(x*x + y*y + z*z);
}

double Vector3f::SquaredNorm() const
{
	return double(x*x + y*y + z*z);
}

double Vector3f::Norm()
{
	double result = double(x*x + y*y + z*z);
	return sqrt(result);
}

double Vector3f::Norm() const 
{
	double result = double(x*x + y*y + z*z);
	return sqrt(result);
}

Vector3f & Vector3f::Zero()
{
	x = 0;
	y = 0;
	z = 0;
	return *this;
}

double Vector3f::GetX()
{
	return x;
}

double Vector3f::GetX() const
{
	return x;
}

void Vector3f::SetX(float value)
{
	x = value;
}

double Vector3f::GetY()
{
	return y;
}

double Vector3f::GetY() const
{
	return y;
}

void Vector3f::SetY(float value)
{
	y = value;
}

double Vector3f::GetZ()
{
	return z;
}

double Vector3f::GetZ() const
{
	return z;
}

void Vector3f::SetZ(float value)
{
	z = value;
}

void Vector3f::SetValue(float xV, float yV, float zV)
{
	x = xV;
	y = yV;
	z = zV;
}

void Vector3f::SetValue(Vector3f & a)
{
	x = a.GetX();
	y = a.GetY();
	z = a.GetZ();
}

//返回向量的标准化
Vector3f Vector3f::Normalize()
{
	double length = sqrt(x*x + y*y + z*z);
	return Vector3f(x / length, y / length, z / length);
}

//返回向量的模
double Vector3f::Length()
{
	double result = sqrt(x*x + y*y + z*z);
	return result;
}

double Vector3f::Dot(Vector3f a)
{
	return double(a.x*x + a.y*y + a.z*z);
}

Vector3f Vector3f::Cross(Vector3f& a)
{
	Vector3f temp;
	temp.SetValue(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
	return temp;
}

Vector3f Vector3f::Cross(const Vector3f& a)
{
	Vector3f temp;
	temp.SetValue(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
	return temp;
}

Vector3f operator*(double x, Vector3f & a)
{
	Vector3f result;
	result.SetX(x*a.GetX());
	result.SetY(x*a.GetY());
	result.SetZ(x*a.GetZ());
	return result;
}

Vector3f operator*(double x, const Vector3f& a)
{
	Vector3f result;
	result.SetX(x*a.GetX());
	result.SetY(x*a.GetY());
	result.SetZ(x*a.GetZ());
	return result;
}

ostream & operator<<(ostream & out, const Vector3f & a)
{
	out << a.GetX() << " " << a.GetY() << " " << a.GetZ();
	return out;
}

Vector3f operator-(const Vector3f & a, const Vector3f & b)
{
	Vector3f result(0, 0, 0);
	result.x = a.x - b.x;
	result.y = a.y - b.y;
	result.z = a.z - b.z;
	return result;
}
