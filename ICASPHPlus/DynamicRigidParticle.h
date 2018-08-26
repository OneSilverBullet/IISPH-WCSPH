#pragma once

#include "../ICASPHPlus/SPH/Vector3f.h"
#include "../ICASPHPlus/SPH/Matrix3f.h"

class DynamicRigidParticle
{
public:
	Vector3f m_x0;  //上一帧的粒子的位置
	Vector3f m_x;   //当前的粒子位置
	Vector3f m_v;   //当前的粒子速度
	double m_boundaryPsi; //当前粒子的boundaryPsi
	Vector3f m_f;   //当前粒子受力
	int hashValue;  //当前粒子的哈希值
	Vector3f diffDistance;

public:
	DynamicRigidParticle() {
		m_x0.Zero();
		m_x.Zero();
		m_v.Zero();
		m_boundaryPsi = 0;
		m_f.Zero();
		diffDistance.Zero();
	}

	DynamicRigidParticle(Vector3f pos, int hash)
	{
		m_x0 = pos;
		m_x = pos;
		m_v.Zero();
		m_boundaryPsi = 0;
		m_f.Zero();
		hashValue = hash;
		diffDistance.Zero();
	}

	~DynamicRigidParticle()
	{

	}

	DynamicRigidParticle& operator=(DynamicRigidParticle& a)
	{
		m_x0 = a.m_x0;
		m_x = a.m_x;
		m_v = a.m_v;
		m_boundaryPsi = a.m_boundaryPsi;
		m_f = a.m_f;
		hashValue = a.hashValue;
		diffDistance = a.diffDistance;
	}


};
