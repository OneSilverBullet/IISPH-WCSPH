#pragma once
#include "Vector3f.h"
#include "StaticRigidObject.h"
#include <vector>

//////////////////////////////////////////////////////////////////////////
//StaticRigidParticle��
//���ã���ɾ�̬��������ӵ��࣬���а�����ģ�������ʹ�õ����ӵı�Ҫ������
//////////////////////////////////////////////////////////////////////////
class StaticRigidParticle
{
public:
	Vector3f x0;          //���ӵĳ�ʼλ��
	Vector3f position;    //���ӵ�λ��
	Vector3f velocity;    //���ӵ��ٶ�
	double   boundaryPsi; //��̬���Ӷ��������ӵ�Ӱ������
	Vector3f m_f;         //��������
	int      hashValue;   //����Ĺ�ϣֵ����������һ��������Ҫ������Ӧ���ӵ�����

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

	//������࣬�ڰ�������Ķ��������ĸ�ֵ�����У��ú�������ʽ����
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
//RigidBodyParticleObject��
//���ã���������ɵľ�̬�����������ģ����������Ĳ���ǽ��
//////////////////////////////////////////////////////////////////////////
class RigidBodyParticleObject
{
public:
	StaticRigidObject rd;                                 //����ĸ�������
	vector<StaticRigidParticle*> staticRigidParticleList;  //��̬���ӵ��б�    
	double supportRadius;                          //��̬���ӵ�֧�Ű뾶�����ڼ���boundaryPSI,����ʹ����С������֧�Ű뾶

public:
	//���ص�ǰ�ɾ�̬��������ɵľ�̬����һ���ж��ٸ�
	int GetStaticParticleNum() { return staticRigidParticleList.size(); }
	RigidBodyParticleObject& operator=(RigidBodyParticleObject& a);
};