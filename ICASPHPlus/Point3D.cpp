#include"Point3D.h"

/*二元运算符：两个三维向量的加法*/
Point3D operator + (const Point3D a, const Point3D b) {
	return Point3D(a.x + b.x, a.y + b.y, a.z + b.z);
}

/*二元运算符：两个三维向量的加法*/
Point3D operator - (const Point3D a, const Point3D b) {
	return Point3D(a.x - b.x, a.y - b.y, a.z - b.z);
}

/*二元运算符：两个三维向量的叉乘*/
Point3D operator * (const Point3D a, const Point3D b) {
	return Point3D(a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x);
}

/*二元运算符：数乘*/
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

/*两个三维向量的点乘*/
double dotProduct(const Point3D a, const Point3D b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*三个点所构成的三角形的面积*/
double getTriangleArea(const Point3D a, const Point3D b, const Point3D c) {
	Point3D crossProduct = (b - a) * (c - b);
	return 0.5 * crossProduct.length();
}

/*两个向量之间构成的角度大小*/
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

/*向量oe和向量os构成的角度大小*/
double angleBetween(Point3D o, Point3D s, Point3D e) {
	return angleBetween(o - e, o - s);
}

/*两点之间的距离*/
double distanc(Point3D a, Point3D b) {
	return (b - a).length();
}

/*判断a, b, c三点是否共线，共线返回true， 否则返回false*/
bool isCollineation(Point3D a, Point3D b, Point3D c) {
	return ((b - a)*(c - a)).length() < EPS;
}

Point3D operator - (const Point3D a)
{
	return Point3D(-a.x, -a.y, -a.z);
}