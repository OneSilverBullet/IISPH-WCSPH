#pragma once
#include "DistanceField.h"
#include "MCTriTable.h"

/***************
*�����࣬��������
*/
class Vec3D
{
public:
	Vec3D()
	{
	}
	Vec3D(double x, double y, double z)
	{
		this->x = x;
		this->y = y;
		this->z = z;
	}
	Vec3D operator+(Vec3D &v)
	{
		return Vec3D(this->x + v.x, this->y + v.y, this->z + v.z);
	}
	Vec3D operator-(Vec3D &v)
	{
		return Vec3D(this->x - v.x, this->y - v.y, this->z - v.z);
	}
	Vec3D operator*(double n)
	{
		return Vec3D(n*this->x, n*this->y, n*this->z);
	}
	Vec3D operator/(double n)
	{
		if (n == 0)
			return Vec3D(MAX_DIS, MAX_DIS, MAX_DIS);
		return Vec3D(this->x / n, this->y / n, this->z / n);
	}
public:
	double x;
	double y;
	double z;
};

/***************
*Marching Cubes�㷨ʵ����
*/
class MCDistField : public DistanceField
{
public:
	MCDistField();
	MCDistField(double len, uint hashSize, double xSize, double ySize, double zSize, double xOrigin = 0.0, double yOrigin = 0.0, double zOrigin = 0.0);
	~MCDistField();

public:
	vector<double> mesh;//mesh���ݣ�ÿ�����������һ�������Ρ�
	vector<double> meshVector;//�㷨������ÿ�������Ӧһ��������
	
	bool ClearGrid();
	//����Mesh������MC�㷨��
	void GenerateMesh();
	//�����ֵ��������ֵ���뵱ǰ��Ԫ��ߵĽ��㣬ͨ���ñ������˵㺯��ֵ�����Բ�ֵ�õ���
	Vec3D VertexInter(Vec3D &p1, Vec3D &p2, double v1, double v2, double isovalue = 0);

	void computeNormal(Vec3D &p1, Vec3D &p2, Vec3D &p3, Vec3D &n);
};