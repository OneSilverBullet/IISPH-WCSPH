#pragma once
#include "HashGridList.h"
#include "GridCell.h"
#include "Vector3f.h"
#include "Vector3i.h"
#include "ShareData.h"
#include "FunctionKit.h"
#include "SPHKernel.h"
#include "WCSPHFluidObject.h"
#include "Boundary.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
using namespace std;


//////////////////////////////////////////////////////////////////////////
//SPHComputation类
//作用：集合了FluidModel以及HashList两个类的对象，实现了SPH的计算。
//////////////////////////////////////////////////////////////////////////
class WCSPHComputation
{
public:
	HashGridList hashGridList; //哈希网格表
	WCSPHFluidObject fluidModel;    //流体表面模型
	Boundary  boundary;        //流体的边界模型
							   //计算流体当中的压强
	double stiffnesss = 50000; //流体当中的刚性
	int gamma = 7;             //流体当中的gamma值
	double viscosityConstant = 0.03; //计算黏度所需要的常量
	double gravity = -9.81;    //流体的重力

public:
	//初始化
	void Initialization(); //在SPH计算之前调用一次，对粒子的属性、位置进行初始化
	void MapParticleToGrid(); //根据粒子的位置进行网格分类
	void MapBoundaryParticleToGrid(); //将边界粒子映射到网格当中
									  //锁帧
	void Frame();   //调用这个函数对帧进行加速
					//计算部分的函数
	void Computation();  //SPH计算，该函数封装下列方法
						 //计算粒子周围邻居粒子
	vector<WCSPHParticle*> ComputeNeighborParticle(WCSPHParticle* origin);
	//计算粒子周围的边界粒子：用于补正粒子密度
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(WCSPHParticle* origin);
	//计算边界粒子周围的边界粒子：用于进行边界粒子PSI计算
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(StaticRigidParticle* origin);
	//根据粒子周围的邻居，计算粒子的密度，并且所有粒子的压强已经计算成功
	double ComputeDensity(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin);
	//根据粒子周围的邻居，计算粒子在压力下的加速度
	Vector3f ComputePressureAcceleration(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin);
	//根据粒子周围的邻居，计算粒子在粘性力下的加速度
	Vector3f ComputeViscostyAcceleration(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin);
	//计算粒子的加速度之和
	Vector3f ComputeAcceleration(Vector3f pressureA, Vector3f viscosityA);
	//将粒子的密度矫正至静止密度和当前密度的最大值
	double CorrectionParticleDensity(WCSPHParticle* a);
	//计算当前粒子的泰特方程得到粒子的压强
	double CalculateParticleTaite(double a);
	//计算流体的邻居表
	void UpdateFluidNeighborList();
	//计算流体的边界邻居表，并且需要不断更新
	void UpdateFluidBoundaryNeighborList();
	//计算边界的边界邻居表，只需要初始化调用一次
	void ComputeBoundaryNeighborList();
};
