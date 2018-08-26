#pragma once
#include <math.h>
#include <math3d.h>
#include <iostream>
using namespace std;

//////////////////////////////////////////////////////////////////////////
//Vector3f��
//���ã����ݽṹ�����ڸ��������Ĳ���
//////////////////////////////////////////////////////////////////////////
class Vector3f
{
public:
	double x, y, z;

public:
	//���캯��
	Vector3f() { x = 0;  y = 0; z = 0; }
	Vector3f(double Vx, double Vy, double Vz) {
		x = Vx;
		y = Vy;
		z = Vz;
	}
	~Vector3f(){}

	//��������
	//����������
	Vector3f& operator=(Vector3f* value);
	double operator*(Vector3f& value);
	Vector3f operator*(double value);
	Vector3f operator*(double value)const;
	Vector3f operator+(Vector3f& value);
	Vector3f operator+(const Vector3f& value);
	Vector3f operator-(Vector3f& value);
	Vector3f& operator-();                    //����ǰ�÷��ŵ�����
	Vector3f operator/(double value);
	Vector3f& operator-=(Vector3f& a);
	Vector3f& operator-=(const Vector3f& a);
	Vector3f& operator+=(Vector3f& a);
	Vector3f& operator+=(const Vector3f& a);
	Vector3f& operator/=(int& a);
	Vector3f& operator/=(const int& a);
	double operator[](const int a);
	bool operator==(Vector3f& value);

	//�����������
	double SquaredNorm();
	double SquaredNorm()const;
	double Norm();
	double Norm() const;
	Vector3f& Zero();
	Vector3f Normalize();     	//���ص�ǰ�����ı�׼��
	double Length();            //���ص�ǰ�����ĳ���
	double Dot(Vector3f a);
	Vector3f Cross(Vector3f& a);
	Vector3f Cross(const Vector3f& a);

	//������ݽӿ�
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

	//��Ԫ����
	friend Vector3f operator*(double x, Vector3f& a);
	friend Vector3f operator*(double x, const Vector3f& a);
	friend ostream& operator<<(ostream & out, const Vector3f& a);
	friend Vector3f operator-(const Vector3f& a, const Vector3f& b);
};