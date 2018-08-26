#pragma once
#include "HashGridList.h"
#include "GridCell.h"
#include "Vector3f.h"
#include "Vector3i.h"
#include "ShareData.h"
#include "FunctionKit.h"
#include "SPHKernel.h"
#include "WCSPHFluidObject.h"
#include "IISPHParticle.h"
#include "IISPHFluidObject.h"
#include "Boundary.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
#include "MCDistField.h"
#include "DynamicRigidParticleObject.h"

using namespace std;

#define MAX_DIS 65536

class IISPHComputer
{
public:
	HashGridList hashGridList;       //哈希网格表
	IISPHFluidObject fluidModel;     //流体表面模型
	Boundary  boundary;              //流体的边界模型
	MCDistField* distanceField ;     //建立流体表面信息抽取
	double viscosityConstant = 0.03; //计算黏度所需要的常量
	double gravity = -9.81;          //流体的重力
	//关于迭代
	int m_iterations = 0;
	double m_maxError = 0.01;
	int m_maxIterations = 100;

	//////////////////////////////////////////////////////////////////////////
	// 无限自适应扩展
	//////////////////////////////////////////////////////////////////////////
	double baseDistance;       //最远的距离
	double fineDistance;//最近距离
	double baseMass = 0.2f;           //基础的质量
	double fineMass = 0.025f;           //最小的质量
	double adaptivityFactor = fineMass/baseMass;   //自适应因子

	
	//工程标志位
	bool fluidAtiFluid = false;
	bool fluidAtiRigid = false;
	bool fluidSurface = false;
	bool fluidAtiDynamicRigid = false;


public:
	//初始化
	void Initialization(); //在SPH计算之前调用一次，对粒子的属性、位置进行初始化

	//IISPH的主要步骤函数：以process函数开头的阶段计算函数
	//计算部分的阶段函数函数
	void Computation();  //SPH计算，该函数封装下列方法
	//1.计算临近粒子列表
	void ProcessNeighborList();
	//2.处理密度计算
	void ProcessDensity();
	//3.处理非压力的力
	void ProcessViscosity();
	//4.更新时间步骤长
	void ProcessUpdateTimeStep();
	//5.预测平流
	void ProcessPredictAdvection();
	//6.压力解算
	void ProcessPressureSolve();
	//7.积分计算
	void ProcessIntegration();

	//IISPH的计算细节函数，下面的函数是为了配合上面的阶段函数进行计算的
	//重置当前所有粒子的加速度
	void ClearAcceleration();
    //计算粒子周围邻居粒子
	vector<IISPHParticle*> ComputeNeighborParticle(IISPHParticle* origin);
	//计算粒子周围的边界粒子：用于补正粒子密度
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(IISPHParticle* origin);
	//计算流体粒子周围的动态刚体粒子：用于补正流体粒子的密度和压力
	vector<DynamicRigidParticle*> ComputeNeighborDynamicParticle(IISPHParticle* origin);
	//计算边界粒子周围的边界粒子：用于进行边界粒子PSI计算
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(StaticRigidParticle* origin);
	//计算动态刚体粒子周围的刚体粒子：用于进行刚体粒子PSI计算
	vector<DynamicRigidParticle*> ComputeNeighborDynamicParticle(DynamicRigidParticle* origin);
	//计算距离域节点周围的边界粒子，用于距离域的计算
	vector<IISPHParticle*> ComputeNeighborFluidParticleDistance2(DistanceNode* origin);
	//根据粒子周围的邻居，计算粒子的密度，并且所有粒子的压强已经计算成功
	double ComputeDensity(vector<IISPHParticle>& neighborList, vector<StaticRigidParticle>& boundaryNeighborList, IISPHParticle& origin);
    //根据粒子周围的邻居，计算粒子在粘性力下的加速度
	Vector3f ComputeViscostyAcceleration(vector<IISPHParticle>& neighborList, IISPHParticle& origin);
	//在计算好压力的情况下，计算粒子的压力加速度
	void ComputePressureAcceleration();
	//计算流体的邻居表
	void UpdateFluidNeighborList();
	//计算流体的边界邻居表，并且需要不断更新
	void UpdateFluidBoundaryNeighborList();
	//计算流体的动态刚体粒子邻居表，并且需要不断的更新
	void UpdateDynamicRigidNeighborsList();
	//计算刚体的刚体边界邻居表，只需要初始化一次（没有做刚体和刚体的碰撞，这里只生成唯一的一个动态刚体）
	void ComputeDynamicNeighborList();
	//计算边界的边界邻居表，只需要初始化调用一次
	void ComputeBoundaryNeighborList();
	//根据粒子的位置进行网格分类
	void MapParticleToGrid(); 
	//将边界粒子映射到网格当中
	void MapBoundaryParticleToGrid(); 
	//将动态刚体粒子映射到网格当中
	void MapDynamicRigidParticleToGrid();
	//基于压力计算得到每个粒子的压力加速度
	void ComputePressureAccels();

	//计算当前pos所对应的距离域值，并且产生液面
	bool ProcessDistanceField(double radius = 0.025);


	//////////////////////////////////////////////////////////////////////////
	// 无限自适应扩展
	//////////////////////////////////////////////////////////////////////////
	//论文中的方法
	void ComputePressureAccels2();
	//计算自适应因子
	void ComputeAdaptivityFactor();
	//计算所有粒子到表面的距离
	void ComputeDistance();
	//根据粒子的mark来对粒子进行分类
	void ClassifyParticle();
	//计算所有粒子的最优粒子质量
	void ComputeOptMass();
	//依据粒子的印记 计算所有粒子的O
	void ComputeO();
	//隐性融合密度
	void BlendDensity();
	//隐性融合速度
	void BlendVelocity();
	//更新所有粒子的融合因子
	void UpdateBlendFactor();

	//add
	//计算每个粒子的距离域DistanceField
	void ComputePtoSandMopt();
	void ProcessMergeandSplit();
};