#pragma once
#include <math.h>
#include "../ICASPHPlus/SPH/ShareData.h"

const double EPS = 1e-10;

/*该类既可以表示一个三维的点，也可以表示一个三维中的向量*/
class Point3D{

	/*点的三维坐标，也可以表示一个三维向量*/
public:
	double x, y, z;

	/*构造函数*/
public:

	/*无参数的构造函数*/
	Point3D(){}

	/*带参数的构造函数*/
	Point3D(double _x, double _y, double _z) {
		x = _x;
		y = _y;
		z = _z;
	}

	/*赋值*/
	void set(double _x, double _y, double _z) {
		x = _x;
		y = _y;
		z = _z;
	}

	/*重载单元运算符*/
public:

	/*单元运算符：三维向量的加等*/
	Point3D operator += (const Point3D others) {
		x += others.x;
		y += others.y;
		z += others.z;
		return Point3D(x, y, z);
	}

	/*单元运算符：三维向量的减等*/
	Point3D operator -= (const Point3D others) {
		x -= others.x;
		y -= others.y;
		z -= others.z;
		return Point3D(x, y, z);
	}

	/*单元运算符：三维向量的数乘*/
	Point3D operator *= (const double times) {
		x *= times;
		y *= times;
		z *= times;
		return Point3D(x, y, z);
	}

	/*单元运算符：三维向量的数除，除数不能为0*/
	Point3D operator /= (const double times) {
		x /= times;
		y /= times;
		z /= times;
		return Point3D(x, y, z);

	}

	/*三维点上的基本运算方法*/
public:

	/*向量的长度*/
	double length() {
		return sqrt(x * x + y * y + z * z);
	}
	/*向量单位化*/
	Point3D normalize() {
		double len = length();
		return Point3D(x / len, y / len, z / len);
	}
};

/*一元运算符：三维向量取反*/
Point3D operator - (const Point3D a);

/*二元运算符：两个三维向量的加法*/
Point3D operator + (const Point3D a, const Point3D b);

/*二元运算符：两个三维向量的减法*/
Point3D operator - (const Point3D a, const Point3D b);

/*二元运算符：两个三维向量的叉乘*/
Point3D operator * (const Point3D a, const Point3D b);

/*二元运算符：三维向量的数乘*/
Point3D operator * (const Point3D a, const double times);

/*二元运算符：三维向量的除法*/
Point3D operator / (const Point3D a, const double times);

/*二元运算符：两个三维向量的相等关系*/
bool operator == (const Point3D a, const Point3D b);

/*二元运算符：两个三维向量的相等关系*/
bool operator != (const Point3D a, const Point3D b);

/*两个三维向量的点乘*/
double dotProduct(const Point3D a, const Point3D b);

/*三个点所构成的三角形的面积*/
double getTriangleArea(const Point3D a, const Point3D b, const Point3D c);

/*两个向量之间构成的角度大小*/
double angleBetween(Point3D a, Point3D b);

/*向量oe和向量os构成的角度大小*/
double angleBetween(Point3D o, Point3D s, Point3D e);

/*两点之间的距离*/
double distanc(Point3D a, Point3D b);

/*判断a, b, c三点是否共线，共线返回true， 否则返回false*/
bool isCollineation(Point3D a, Point3D b, Point3D c);
