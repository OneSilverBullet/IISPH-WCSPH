#pragma once

#include "DynamicRigidObject.h";
#include "DynamicRigidParticle.h"
#include "../ICASPHPlus/SPH/Vector3i.h"
#include "../ICASPHPlus/SPH/FunctionKit.h"
#include <vector>

class DynamicRigidParticleObject
{
public:
	double m_h = 0.003f;
	vector<DynamicRigidParticle*> dynamicParticles;
	DynamicRigidObject m_rigidBody;

	//当前刚体的初始位置
	double ox = 0.5f;
	double oy = 0.5f;
	double oz = 0.6f;
	//当前生成刚体正方体的长度
	double length=0.2f;
	//当前生成的刚体正方体的粒子半径
	double radius = 0.025f;
	//当前用粒子生成的物体的粒子质量
	double particleMass = 0.1f;
	//当前生成刚体正方体的长宽高方向一共有多少个粒子
	Vector3i dynamicCubeWLH;

	//边界偏移值
	double boundaryOffset = 0.3f;


public:

	DynamicRigidParticleObject() {
		dynamicParticles.resize(0);
	}

	DynamicRigidParticleObject(DynamicRigidObject rigidBody): m_rigidBody(rigidBody), m_h(0.0) {}

	~DynamicRigidParticleObject() {}

	void updateTimeStepSize() { m_h = 0.001; }

	virtual bool isDynamic() const { return m_rigidBody.getMass() != 0.0; }

	virtual double const getMass() const { return m_rigidBody.getMass(); }
	virtual Vector3f const& getPosition() { return m_rigidBody.getPosition(); }
	virtual Vector3f const& getVelocity() const { return m_rigidBody.getVelocity(); }
	virtual Matrix3f const& getRotation() const { return m_rigidBody.getRotationMatrix(); }
	virtual Vector3f const& getAngularVelocity() const { return m_rigidBody.getAngularVelocity(); }
	virtual void addForce(const Vector3f &f) 
	{ 
		//包含重力加速度和粒子给当前刚体的加速度
		m_rigidBody.getAcceleration() += (1.0 / m_rigidBody.getMass()) * f;
		//根据当前加速度对当前的速度进行叠加
		m_rigidBody.getVelocity() += m_rigidBody.getAcceleration() * m_h;
	}
	virtual void addTorque(const Vector3f &t) { 
		m_rigidBody.getAngularVelocity() += m_rigidBody.getInertiaTensorInverseW() * t * m_h;
	}

	//在上面的基础上增添的函数
	//所有的基础上进行初始化
	void Initialize();
	//初始化粒子
	void InitParticle();
	//初始化当前刚体的属性
	void InitRigidBody();
	//计算当前刚体的惯性张量
	void CalculateInertiaTensor();
	//根据当前角速度求出当前刚体的方向矩阵
	void CalculateRotationMat();
	//设置粒子与刚体粒子的相对差距
	void CalculateDiffDistance();
	//对当前刚体的属性进行计算0
	void UpdateObjectPar();
	//对当前粒子的属性进行计算
	void UpdateParticlePar();
	//当前刚体刷新的一帧
	void Calculation();



};