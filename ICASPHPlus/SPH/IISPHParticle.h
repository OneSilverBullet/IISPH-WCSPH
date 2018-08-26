#pragma once
#include "Vector3f.h"
#include "Boundary.h"
#include "FunctionKit.h"
#include "StaticRigidObject.h"
#include "ParticleObject.h"
#include "ShareData.h"
#include "DynamicRigidParticle.h"
#include <vector>

class IISPHParticle
{
public:
	//通用粒子属性
	Vector3f position;        //粒子的位置
	Vector3f velocity;        //粒子的速度
	Vector3f acceleration;    //粒子的加速度
	double particleRad;        //粒子半径，固定大小
	double particleSupportRad; //粒子的支撑半径，其为粒子半径的四倍
	double particleMass;        //每个粒子的质量：由当前的密度与粒子半径计算得到
	double density;           //粒子的密度
	Vector3f viscosity;       //粒子的黏度力加速度
	Vector3f pressureAcceleration;        //粒子的压力加速度
	int hashIndex;            //当前粒子所属的网格哈希号

	//IISPH所用到的粒子属性
	double aii;
	Vector3f dii;
	Vector3f dij_pj;
	double density_adv;
	double pressure;
	double lastPressure;
	//依据粒子当前的大小动态的检索周围的网格
	int gridCellNum = 0;

	//add
	double distValue;
	double disToSurface;

	//Infinit Continous扩展用到的粒子属性
	double mopt;             //当前粒子距离的最优质量
	double distance;         //距离流体表面的距离
	double mark;             //每个粒子的标记
	double O;                //当前每个粒子都有O，用于适应不同支撑半径粒子之间的交互
	PARTICLE_TYPE pt;        //当前粒子的粒子类型
	double lastDensity;      //当前粒子的上一次的密度值
	Vector3f lastVelocity;     //当前粒子的上一次的速度值
	double β;               //当前粒子的密度混合、速度混合因子
	bool NEW_UPDATE = false; //当前粒子更新的状态，一旦粒子分裂或者合并，那么该值设置为true

	//关于属性的刷新修正
	vector<IISPHParticle*> fluidNeighbors;                //这里是指向流体邻居的指针
	vector<StaticRigidParticle*> boundaryNeighbors;       //这里是指向边界邻居的指针
	vector<DynamicRigidParticle*> dynamicRigidNeighbors;  //这里是指向动态刚体的指针

public:
	//构造函数直接调用初始化函数
	IISPHParticle() {
		//粒子正常属性初始化
		position.Zero();       
		velocity.Zero();
		acceleration.Zero();
		particleRad = 0.025f;
		particleSupportRad = 0.1f; 
		particleMass = 0.1f;       
		density = 0.0f;          
		viscosity.Zero();       
		pressureAcceleration.Zero();        
		hashIndex = 0;            

		//IISPH属性初始化
		aii = 0.0;
		dii.Zero();
		dij_pj.Zero();
		density_adv = 0.0;
		pressure = 0.0;
		lastPressure = 0.0;
		pressureAcceleration.Zero();

		//无限自适应的扩展属性初始化
		distance = 0.0f;
		mark = 0.0f;
		mopt = 0.0f;
		distance = 0.0f;
		mark = 0;
		O = 0;
		pt = PARTICLE_TYPE::S;        
		lastDensity = 0;      
		lastVelocity.Zero();
		β = 0;
		//临近表扩展
		fluidNeighbors.resize(0);
		boundaryNeighbors.resize(0);
		dynamicRigidNeighbors.resize(0);
	}
	~IISPHParticle() {}
	void Initialization(Vector3f position, int v, double radius); //对粒子进行初始化更新
	void InitializationPVS(Vector3f position, int v, double mass , double density);
	IISPHParticle& operator=(IISPHParticle& a);    //对粒子进行赋值操作重载
	void Intigration(Boundary bound, double timeStep);       //依据加速度对粒子的速度和位置进行积分求值
	void ParticlePositionCorrection(Boundary bound);               //对当前粒子的位置进行校正

};