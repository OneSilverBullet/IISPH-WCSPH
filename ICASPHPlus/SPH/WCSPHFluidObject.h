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
//FluidObject��
//���ã������ģ���ࡣ�������������ӹ��ɡ�����һЩ����Ļ��������Լ�
//����������ÿ�����ӵ�������λ�á�
//////////////////////////////////////////////////////////////////////////
class WCSPHFluidObject
{
public:
	double timeStep = 0.001f;           //ֻ��������ģ���ʱ��ʹ��ʱ�䲽��������
										//FluidObject�����ӵĹ̶�������ȡ����������ÿ�����ӵĶ����С
										//���ӵĹ̶�����
	double particleRad = 0.025f;        //���Ӱ뾶���̶���С
	double particleSupportRad = 0.1f;   //���ӵ�֧�Ű뾶����Ϊ���Ӱ뾶���ı�
	double particleMass = 0.1f;         //ÿ�����ӵ��������ɵ�ǰ���ܶ������Ӱ뾶����õ�
	double restDensity = 1000.f;        //����ľ�ֹ�ܶȣ��ܻ��õ��Ĳ���
	int    particleNum;                 //����������Ҫ�������
										//����Ĺ̶�����
	double fluidWitdth = 0.5f;
	double fluidLength = 0.5f;
	double fluidHeight = 0.5f;          //ˮ��Ĵ�С
	double originX = 0.25f;
	double originY = 0.25f;
	double originZ = 0.f;               //ˮ�����ʼ��Ⱦλ��
	Vector3i fluidWLH;                  //ˮ����X��Y��Z�����ϵ����Ӹ���
										//���ӵ��б�
	vector<WCSPHParticle*> particleList;          //��ǰˮ��������б�
	vector<vector<WCSPHParticle*>> neighborList;  //ˮ���ÿ�����ӵ��ھӱ����������ڲ��������Լ���
	Boundary boundary;                      //���ñ߽������

											//////////////////////////////////////////////////////////////////////////
											//***************���ӱ߽�***********************************************//
											//////////////////////////////////////////////////////////////////////////
	RigidBodyParticleObject boundaryObj;                            //������װ�ɵı߽����
	vector<vector<StaticRigidParticle*>> boundaryNeighborList;       //�߽����ӵ��ھӣ����ڱ߽�PSI����
	vector<vector<StaticRigidParticle*>> fluidBoundaryNeighborList;  //����ı߽������ھӣ����������ܶȲ���


public:
	//�Ե�ǰˮ����г�ʼ��
	void Initialise();                //����Ĭ�����Խ��г�ʼ��
	void Initialise(Vector3f WLH, Vector3f fluidOrigin, double rD, double pR, double time, Boundary bound);
	void ComputeFluidWLH();
	void ComputeFluidWLH(double fluidW, double fluidL, double fluidH, double r);
	//����ˮ��ĳ���߷����ϵ����ж�������
	void InitialiseParticlePosition();//��ʼ��ˮ�������λ��
	void InitialiseParticleDensity(); //��ʼ��ˮ����ܶ�
	void InitialiseParticleP();       //��ʼ��ˮ���ѹǿ

	void UpdateParticlePosition();    //�������ӵ�λ��
	void UpdateParticleHashValue();   //����ÿ�����ӵĹ�ϣֵ
	void ParticleCorrection(WCSPHParticle* a);        //�������ӵ��ٶ��Լ���ǰλ�ý���У��
	void ComputeParticleNum(int x, int y, int z);

	//�������ӵİ뾶��ͬʱ���¹⻬�˺�����
	void SetParticleRad();          //ʹ��Ĭ�ϲ���
	void SetParticleRad(double r);  //�����Ӱ뾶���и���

									//�����ھӱ�Ĳ���
									//���и����ھӱ�Ĳ��������Computer���н���
	void ClearNeighborList();
	//�������������Χ�ı߽������ھ�
	void ClearFluidBoundaryNeighborList();

	//////////////////////////////////////////////////////////////////////////
	//***************���Ӳ���***********************************************//
	//////////////////////////////////////////////////////////////////////////

	//�����б߽����ӽ��г�ʼ��
	void InitializeBoundary();
	//ͨ���������Ӷ�������еı߽����ӵĸ���
	void AddRigidBodyObject(StaticRigidObject& rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles);
	void ComputeBoundaryPsi(RigidBodyParticleObject& bound);

};