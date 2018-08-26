#pragma once
#include "Vector3f.h"
#include "StaticRigidObject.h"
#include <vector>

//////////////////////////////////////////////////////////////////////////
//StaticRigidParticle类
//作用：组成静态刚体的粒子的类，其中包含了模拟刚体所使用的粒子的必要的属性
//////////////////////////////////////////////////////////////////////////
class StaticRigidParticle
{
public:
	Vector3f x0;          //粒子的初始位置
	Vector3f position;    //粒子的位置
	Vector3f velocity;    //粒子的速度
	double   boundaryPsi; //静态粒子对流体粒子的影响因子
	Vector3f m_f;         //粒子受力
	int      hashValue;   //这里的哈希值和其他流体一样，都需要经历对应粒子的搜索

public:
	StaticRigidParticle()
	{
		x0.Zero();
		position.Zero();
		velocity.Zero();
		boundaryPsi = 0;
		m_f.Zero();
		hashValue = 0;
	}

	StaticRigidParticle(StaticRigidParticle& a)
	{
		operator=(a);
	}

	~StaticRigidParticle()
	{

	}

	//必须的类，在包含该类的对象的物体的赋值操作中，该函数被隐式调用
	StaticRigidParticle& operator=(StaticRigidParticle& a)
	{
		x0 = a.x0;
		position = a.position;
		velocity = a.velocity;
		boundaryPsi = a.boundaryPsi;
		m_f = a.m_f;
		hashValue = a.hashValue;
		return *this;
	}
};


//////////////////////////////////////////////////////////////////////////
//RigidBodyParticleObject类
//作用：由粒子组成的静态刚体对象。用于模拟流体四面的玻璃墙。
//////////////////////////////////////////////////////////////////////////
class RigidBodyParticleObject
{
public:
	StaticRigidObject rd;                                 //刚体的各种属性
	vector<StaticRigidParticle*> staticRigidParticleList;  //静态粒子的列表    
	double supportRadius;                          //静态粒子的支撑半径，用于计算boundaryPSI,这里使用最小的粒子支撑半径

public:
	//返回当前由静态粒子所组成的静态物体一共有多少个
	int GetStaticParticleNum() { return staticRigidParticleList.size(); }
	RigidBodyParticleObject& operator=(RigidBodyParticleObject& a);
};