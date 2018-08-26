#include"Point3D.h"

/*��Ԫ�������������ά�����ļӷ�*/
Point3D operator + (const Point3D a, const Point3D b) {
	return Point3D(a.x + b.x, a.y + b.y, a.z + b.z);
}

/*��Ԫ�������������ά�����ļӷ�*/
Point3D operator - (const Point3D a, const Point3D b) {
	return Point3D(a.x - b.x, a.y - b.y, a.z - b.z);
}

/*��Ԫ�������������ά�����Ĳ��*/
Point3D operator * (const Point3D a, const Point3D b) {
	return Point3D(a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

/*��Ԫ�����������*/
Point3D operator * (const Point3D a, const double times)
{
	return Point3D(a.x * times, a.y * times, a.z * times);
}

Point3D operator / (const Point3D a, const double times)
{
	return Point3D(a.x / times, a.y / times, a.z / times);
}

bool operator == (const Point3D a, const Point3D b)
{
	if (a.x == b.x&&a.y == b.y&&a.z == b.z)
	{
		return true;
	}

	return false;
}

bool operator != (const Point3D a, const Point3D b)
{
	if (a.x == b.x&&a.y == b.y&&a.z == b.z)
	{
		return false;
	}

	return true;
}

/*������ά�����ĵ��*/
double dotProduct(const Point3D a, const Point3D b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*�����������ɵ������ε����*/
double getTriangleArea(const Point3D a, const Point3D b, const Point3D c) {
	Point3D crossProduct = (b - a) * (c - b);
	return 0.5 * crossProduct.length();
}

/*��������֮�乹�ɵĽǶȴ�С*/
double angleBetween(Point3D a, Point3D b) {
	double cosAngle = dotProduct(a, b) / (a.length() * b.length());
	if (cosAngle - 1 > -EPS)
	{
		cosAngle = 1;
	}
	else if (cosAngle + 1 < EPS)
	{
		cosAngle = -1;
	}
	return acos(cosAngle);
}

/*����oe������os���ɵĽǶȴ�С*/
double angleBetween(Point3D o, Point3D s, Point3D e) {
	return angleBetween(o - e, o - s);
}

/*����֮��ľ���*/
double distanc(Point3D a, Point3D b) {
	return (b - a).length();
}

/*�ж�a, b, c�����Ƿ��ߣ����߷���true�� ���򷵻�false*/
bool isCollineation(Point3D a, Point3D b, Point3D c) {
	return ((b - a)*(c - a)).length() < EPS;
}

Point3D operator - (const Point3D a)
{
	return Point3D(-a.x, -a.y, -a.z);
}