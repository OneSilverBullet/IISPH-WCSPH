#pragma once
#include "WCSPHParticle.h"
#include "Vector3f.h"
#include "Vector3i.h"
#include "SPHKernel.h"
#include "Boundary.h"
#include "HashGridList.h"
#include "StaticRigidObject.h"
#include "ParticleObject.h"
#include <vector>

//////////////////////////////////////////////////////////////////////////
//FluidObject类
//作用：流体的模拟类。其中由流体粒子构成。包含一些计算的基本属性以及
//函数。更新每个粒子的属性与位置。
//////////////////////////////////////////////////////////////////////////
class WCSPHFluidObject
{
public:
	double timeStep = 0.001f;           //只有在流体模拟的时候，使用时间步骤来进行
										//FluidObject将粒子的固定属性提取出来，减少每个粒子的对象大小
										//粒子的固定属性
	double particleRad = 0.025f;        //粒子半径，固定大小
	double particleSupportRad = 0.1f;   //粒子的支撑半径，其为粒子半径的四倍
	double particleMass = 0.1f;         //每个粒子的质量：由当前的密度与粒子半径计算得到
	double restDensity = 1000.f;        //流体的静止密度，总会用到的参数
	int    particleNum;                 //粒子数量需要计算出来
										//流体的固定属性
	double fluidWitdth = 0.5f;
	double fluidLength = 0.5f;
	double fluidHeight = 0.5f;          //水体的大小
	double originX = 0.25f;
	double originY = 0.25f;
	double originZ = 0.f;               //水体的起始渲染位置
	Vector3i fluidWLH;                  //水体在X、Y、Z方向上的粒子个数
										//粒子的列表
	vector<WCSPHParticle*> particleList;          //当前水体的粒子列表
	vector<vector<WCSPHParticle*>> neighborList;  //水体的每个粒子的邻居表：用于流体内部进行属性计算
	Boundary boundary;                      //设置边界的属性

											//////////////////////////////////////////////////////////////////////////
											//***************粒子边界***********************************************//
											//////////////////////////////////////////////////////////////////////////
	RigidBodyParticleObject boundaryObj;                            //粒子组装成的边界对象
	vector<vector<StaticRigidParticle*>> boundaryNeighborList;       //边界粒子的邻居：用于边界PSI计算
	vector<vector<StaticRigidParticle*>> fluidBoundaryNeighborList;  //流体的边界粒子邻居：用于流体密度补正


public:
	//对当前水体进行初始化
	void Initialise();                //根据默认属性进行初始化
	void Initialise(Vector3f WLH, Vector3f fluidOrigin, double rD, double pR, double time, Boundary bound);
	void ComputeFluidWLH();
	void ComputeFluidWLH(double fluidW, double fluidL, double fluidH, double r);
	//计算水体的长宽高方向上到底有多少粒子
	void InitialiseParticlePosition();//初始化水体的粒子位置
	void InitialiseParticleDensity(); //初始化水体的密度
	void InitialiseParticleP();       //初始化水体的压强

	void UpdateParticlePosition();    //更新粒子的位置
	void UpdateParticleHashValue();   //更新每个粒子的哈希值
	void ParticleCorrection(WCSPHParticle* a);        //依据粒子的速度以及当前位置进行校正
	void ComputeParticleNum(int x, int y, int z);

	//更新粒子的半径，同时更新光滑核函数类
	void SetParticleRad();          //使用默认参数
	void SetParticleRad(double r);  //对粒子半径进行更改

									//关于邻居表的操作
									//其中更新邻居表的操作务必在Computer当中进行
	void ClearNeighborList();
	//清除流体粒子周围的边界粒子邻居
	void ClearFluidBoundaryNeighborList();

	//////////////////////////////////////////////////////////////////////////
	//***************粒子操作***********************************************//
	//////////////////////////////////////////////////////////////////////////

	//对所有边界粒子进行初始化
	void InitializeBoundary();
	//通过刚体计算从而完成所有的边界粒子的更新
	void AddRigidBodyObject(StaticRigidObject& rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles);
	void ComputeBoundaryPsi(RigidBodyParticleObject& bound);

};