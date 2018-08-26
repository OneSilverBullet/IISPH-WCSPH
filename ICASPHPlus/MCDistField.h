#pragma once
#include "DistanceField.h"
#include "MCTriTable.h"

/***************
*顶点类，或向量类
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
*Marching Cubes算法实现类
*/
class MCDistField : public DistanceField
{
public:
	MCDistField();
	MCDistField(double len, uint hashSize, double xSize, double ySize, double zSize, double xOrigin = 0.0, double yOrigin = 0.0, double zOrigin = 0.0);
	~MCDistField();

public:
	vector<double> mesh;//mesh数据，每三个顶点组成一个三角形。
	vector<double> meshVector;//点法向量，每个顶点对应一个向量。
	
	bool ClearGrid();
	//构造Mesh，根据MC算法。
	void GenerateMesh();
	//顶点插值函数。等值面与当前单元体边的交点，通过该边两个端点函数值的线性插值得到。
	Vec3D VertexInter(Vec3D &p1, Vec3D &p2, double v1, double v2, double isovalue = 0);

	void computeNormal(Vec3D &p1, Vec3D &p2, Vec3D &p3, Vec3D &n);
};