#pragma once
#include "Vector3f.h"
#include "IISPHParticle.h"
#include "StaticRigidObject.h"
#include "ParticleObject.h"
#include "Boundary.h"
#include "FunctionKit.h"
#include "SPHKernel.h"
#include "ShareData.h"
#include "DynamicRigidParticle.h"
#include "DynamicRigidParticleObject.h"
#include <vector>



class IISPHFluidObject
{
public:
	//ʱ�䲽��ĳ�ʼ��Ϣ
	double IISPH_TimeStep = 0.003f;           //ֻ��������ģ���ʱ��ʹ��ʱ�䲽��������
										      //FluidObject�����ӵĹ̶�������ȡ����������ÿ�����ӵĶ����С
						
	//���ӵĹ̶�����
	double IISPH_RestDensity = 1000.f;        //����ľ�ֹ�ܶȣ��ܻ��õ��Ĳ���
	int    IISPH_ParticleNum;                 //����������Ҫ�������

	//ˮ����Ⱦ�ĳ�ʼ��Ϣ
	double IISPH_FluidWitdth = 1.0f;
	double IISPH_FluidLength = 1.0f;
	double IISPH_FluidHeight = 1.0f;          //ˮ��Ĵ�С
	double IISPH_OriginX = 0.25f;
	double IISPH_OriginY = 0.25f;
	double IISPH_OriginZ = 0.025f;              //ˮ�����ʼ��Ⱦλ��									

	//˫����
	bool doubleDamBreak = false;                //�Ƿ���˫����ģ��
											   //��ҵ��Ƶ��е�˫��������ģ��
	double IISPH_FluidWitdth2 = 0.5f;
	double IISPH_FluidLength2 = 0.5f;
	double IISPH_FluidHeight2 = 0.5f;          //ˮ��Ĵ�С
	double IISPH_OriginX2 = 0.25f;
	double IISPH_OriginY2 = 1.25f;
	double IISPH_OriginZ2 = 0.025f;              //ˮ�����ʼ��Ⱦλ��		
	Vector3i FluidWLH2;

	//��̬����
	bool staticPilarFlag = false;
	int boundaryNumber = 0;   //��¼��ǰ�߽������

	//�߽������Ϣ
	Boundary boundary;                                             //���ñ߽�����ԣ����ڼ���߽����ӵ�λ�ó�ʼ��
	RigidBodyParticleObject boundaryObj;                           //������װ�ɵľ�̬����߽����
	DynamicRigidParticleObject dynamicCube;                       //��̬����

	//�ٽ���
	vector<IISPHParticle*> particleList;                             //��ǰˮ��������б�
	//vector<vector<IISPHParticle>> neighborList;                     //ˮ���ÿ�����ӵ��ھӱ����������ڲ��������Լ���
	vector<vector<StaticRigidParticle*>> boundaryNeighborList;       //�߽����ӵ��ھӣ����ڱ߽�PSI����
	//vector<vector<StaticRigidParticle>> fluidBoundaryNeighborList;  //����ı߽������ھӣ����������ܶȲ���
	vector<vector<DynamicRigidParticle*>> dynamicRigidNeighborList;   //����Ķ�̬�����ھӣ����ڶ�̬����PSI����


	//ˮ���ʼ������Ϣ
	Vector3i FluidWLH;          //�����ˮ���ʼ����ʱ�򣬳���߷������ж��ٸ�����,ֻ�ڳ�ʼ��ʱ������
	double particleInitialPad = 0.025f;  //ˮ�����ӵĳ�ʼ���뾶������֮�����ӵĴ�С����Ӧ��������ֵֻ�ڳ�ʼ����ʱ������


public:
	//�Ե�ǰˮ����г�ʼ��
	void Initialise();                //����Ĭ�����Խ��г�ʼ��

	void ComputeFluidWLH();
	void InitialiseParticlePosition();//��ʼ��ˮ�������λ��
	void InitialiseParticleDensity(); //��ʼ��ˮ����ܶ�

	void UpdateParticleHashValue();   //����ÿ�����ӵĹ�ϣֵ
	void ComputeParticleNum(int x, int y, int z);
	void ClearNeighborList();
	//�������������Χ�ı߽������ھ�
	void ClearFluidBoundaryNeighborList();
	//�������������Χ�Ķ�̬���������ھ�
	void ClearDynamicRigidNeighborList();

	//////////////////////////////////////////////////////////////////////////
	//***************���Ӳ���***********************************************//
	//////////////////////////////////////////////////////////////////////////
	//�����б߽����ӽ��г�ʼ��
	void InitializeBoundary();
	//ͨ���������Ӷ�������еı߽����ӵĸ���
	void AddRigidBodyObject(StaticRigidObject& rbo, const unsigned int numBoundaryParticles, vector<Vector3f> boundaryParticles);
	void ComputeBoundaryPsi(RigidBodyParticleObject& bound);

	//////////////////////////////////////////////////////////////////////////
	//****************��̬���Ӳ���*****************************************///
	//////////////////////////////////////////////////////////////////////////
	void ComputeBoundaryPsi(DynamicRigidParticleObject& cube);


};