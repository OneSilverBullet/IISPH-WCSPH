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

	//��ǰ����ĳ�ʼλ��
	double ox = 0.5f;
	double oy = 0.5f;
	double oz = 0.6f;
	//��ǰ���ɸ���������ĳ���
	double length=0.2f;
	//��ǰ���ɵĸ�������������Ӱ뾶
	double radius = 0.025f;
	//��ǰ���������ɵ��������������
	double particleMass = 0.1f;
	//��ǰ���ɸ���������ĳ���߷���һ���ж��ٸ�����
	Vector3i dynamicCubeWLH;

	//�߽�ƫ��ֵ
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
		//�����������ٶȺ����Ӹ���ǰ����ļ��ٶ�
		m_rigidBody.getAcceleration() += (1.0 / m_rigidBody.getMass()) * f;
		//���ݵ�ǰ���ٶȶԵ�ǰ���ٶȽ��е���
		m_rigidBody.getVelocity() += m_rigidBody.getAcceleration() * m_h;
	}
	virtual void addTorque(const Vector3f &t) { 
		m_rigidBody.getAngularVelocity() += m_rigidBody.getInertiaTensorInverseW() * t * m_h;
	}

	//������Ļ���������ĺ���
	//���еĻ����Ͻ��г�ʼ��
	void Initialize();
	//��ʼ������
	void InitParticle();
	//��ʼ����ǰ���������
	void InitRigidBody();
	//���㵱ǰ����Ĺ�������
	void CalculateInertiaTensor();
	//���ݵ�ǰ���ٶ������ǰ����ķ������
	void CalculateRotationMat();
	//����������������ӵ���Բ��
	void CalculateDiffDistance();
	//�Ե�ǰ��������Խ��м���0
	void UpdateObjectPar();
	//�Ե�ǰ���ӵ����Խ��м���
	void UpdateParticlePar();
	//��ǰ����ˢ�µ�һ֡
	void Calculation();



};