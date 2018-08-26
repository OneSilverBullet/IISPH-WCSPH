#pragma once

#include "Vector3f.h"
#include "Vector3i.h"
#include "StaticRigidObject.h"
#include "ShareData.h"

//////////////////////////////////////////////////////////////////////////
//Boundary类
//作用：包含一些边界属性。根据边界属性采样得到当前边界的粒子位置数组
//////////////////////////////////////////////////////////////////////////
class Boundary
{
public:
	//边界的属性
	double boundaryOriginX = 0;
	double boundaryOriginY = 0;
	double boundaryOriginZ = 0; //用于边界测试的边界起始点
	double boundaryWidth = 1;
	double boundaryLength = 2;
	double boundaryHeight = 1;  //用于边界测试的边界长宽高
	//避免访问非法内存的解决方案
	double  offset =  0.3f;               //用于生成网格的部分的偏移，一般为两个网格的长度
	double gridBoundaryOriginX = -0.5f;
	double gridBoundaryOriginY = -0.5f;
	double gridBoundaryOriginZ = -0.5f; //用于生成网格的边界起始点
	double gridBoundaryWidth = 3.0f;
	double gridBoundaryLength = 5.0f;
	double gridBoundaryHeight = 3.0f;   //用于生成网格的长宽高
	
	double distanceFieldOffset = 0.2f;
	double distanceFieldOriginX = -0.2f;
	double distanceFieldOriginY = -0.2f;
	double distanceFieldOriginZ = -0.2f;
	double distanceFieldWidth = 1.4f;
	double distanceFieldLength = 2.4f;
	double distanceFieldHeight = 1.4f;

	double pRad = 0.025f;

	//是否绘制静态刚体
	bool staticRigidFlag = false;
	int boundaryNumber = 0;  //边界的粒子的数量

public:
	//构造析构函数
	Boundary();
	Boundary(Boundary& a);
	Boundary(double x, double y, double z, double w, double l, double h, double os);
	~Boundary();
	//依据特定的值进行新的网格生成
	void GenerateNewBoundary(double x, double y, double z, double w, double l, double h, double os); //生成新的网格
	void GenerateNewBoundary(Vector3f position, Vector3f fluidWLH, double offs);
	Boundary& operator=(Boundary& a);
	//依据采样得到的边界数据采样计算得到当前的边界粒子数组
	vector<Vector3f> initBoundaryData();
};