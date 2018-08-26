#pragma once
#include <math.h>
#include <math3d.h>
#include <iostream>
using namespace std;

//////////////////////////////////////////////////////////////////////////
//Vector3f类
//作用：数据结构，用于各种向量的测算
//////////////////////////////////////////////////////////////////////////
class Vector3f
{
public:
	double x, y, z;

public:
	//构造函数
	Vector3f() { x = 0;  y = 0; z = 0; }
	Vector3f(double Vx, double Vy, double Vz) {
		x = Vx;
		y = Vy;
		z = Vz;
	}
	~Vector3f(){}

	//基本运算
	//操作符重载
	Vector3f& operator=(Vector3f* value);
	double operator*(Vector3f& value);
	Vector3f operator*(double value);
	Vector3f operator*(double value)const;
	Vector3f operator+(Vector3f& value);
	Vector3f operator+(const Vector3f& value);
	Vector3f operator-(Vector3f& value);
	Vector3f& operator-();                    //关于前置符号的重载
	Vector3f operator/(double value);
	Vector3f& operator-=(Vector3f& a);
	Vector3f& operator-=(const Vector3f& a);
	Vector3f& operator+=(Vector3f& a);
	Vector3f& operator+=(const Vector3f& a);
	Vector3f& operator/=(int& a);
	Vector3f& operator/=(const int& a);
	double operator[](const int a);
	bool operator==(Vector3f& value);

	//特殊操作函数
	double SquaredNorm();
	double SquaredNorm()const;
	double Norm();
	double Norm() const;
	Vector3f& Zero();
	Vector3f Normalize();     	//返回当前向量的标准化
	double Length();            //返回当前向量的长度
	double Dot(Vector3f a);
	Vector3f Cross(Vector3f& a);
	Vector3f Cross(const Vector3f& a);

	//类的数据接口
	double GetX();
	double GetX() const;
	void SetX(float value);
	double GetY();
	double GetY() const;
	void SetY(float value);
	double GetZ();
	double GetZ() const;
	void SetZ(float value);
	void SetValue(float xV, float yV, float zV);
	void SetValue(Vector3f& a);

	//友元函数
	friend Vector3f operator*(double x, Vector3f& a);
	friend Vector3f operator*(double x, const Vector3f& a);
	friend ostream& operator<<(ostream & out, const Vector3f& a);
	friend Vector3f operator-(const Vector3f& a, const Vector3f& b);
};