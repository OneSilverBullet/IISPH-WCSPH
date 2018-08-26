#pragma once

#include "../ICASPHPlus/SPH/Vector3f.h"
#include "../ICASPHPlus/SPH/Matrix3f.h"

class DynamicRigidParticle
{
public:
	Vector3f m_x0;  //��һ֡�����ӵ�λ��
	Vector3f m_x;   //��ǰ������λ��
	Vector3f m_v;   //��ǰ�������ٶ�
	double m_boundaryPsi; //��ǰ���ӵ�boundaryPsi
	Vector3f m_f;   //��ǰ��������
	int hashValue;  //��ǰ���ӵĹ�ϣֵ
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
