#pragma once
#include "HashGridList.h"
#include "GridCell.h"
#include "Vector3f.h"
#include "Vector3i.h"
#include "ShareData.h"
#include "FunctionKit.h"
#include "SPHKernel.h"
#include "WCSPHFluidObject.h"
#include "IISPHParticle.h"
#include "IISPHFluidObject.h"
#include "Boundary.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
#include "MCDistField.h"
#include "DynamicRigidParticleObject.h"

using namespace std;

#define MAX_DIS 65536

class IISPHComputer
{
public:
	HashGridList hashGridList;       //��ϣ�����
	IISPHFluidObject fluidModel;     //�������ģ��
	Boundary  boundary;              //����ı߽�ģ��
	MCDistField* distanceField ;     //�������������Ϣ��ȡ
	double viscosityConstant = 0.03; //����������Ҫ�ĳ���
	double gravity = -9.81;          //���������
	//���ڵ���
	int m_iterations = 0;
	double m_maxError = 0.01;
	int m_maxIterations = 100;

	//////////////////////////////////////////////////////////////////////////
	// ��������Ӧ��չ
	//////////////////////////////////////////////////////////////////////////
	double baseDistance;       //��Զ�ľ���
	double fineDistance;//�������
	double baseMass = 0.2f;           //����������
	double fineMass = 0.025f;           //��С������
	double adaptivityFactor = fineMass/baseMass;   //����Ӧ����

	
	//���̱�־λ
	bool fluidAtiFluid = false;
	bool fluidAtiRigid = false;
	bool fluidSurface = false;
	bool fluidAtiDynamicRigid = false;


public:
	//��ʼ��
	void Initialization(); //��SPH����֮ǰ����һ�Σ������ӵ����ԡ�λ�ý��г�ʼ��

	//IISPH����Ҫ���躯������process������ͷ�Ľ׶μ��㺯��
	//���㲿�ֵĽ׶κ�������
	void Computation();  //SPH���㣬�ú�����װ���з���
	//1.�����ٽ������б�
	void ProcessNeighborList();
	//2.�����ܶȼ���
	void ProcessDensity();
	//3.�����ѹ������
	void ProcessViscosity();
	//4.����ʱ�䲽�賤
	void ProcessUpdateTimeStep();
	//5.Ԥ��ƽ��
	void ProcessPredictAdvection();
	//6.ѹ������
	void ProcessPressureSolve();
	//7.���ּ���
	void ProcessIntegration();

	//IISPH�ļ���ϸ�ں���������ĺ�����Ϊ���������Ľ׶κ������м����
	//���õ�ǰ�������ӵļ��ٶ�
	void ClearAcceleration();
    //����������Χ�ھ�����
	vector<IISPHParticle*> ComputeNeighborParticle(IISPHParticle* origin);
	//����������Χ�ı߽����ӣ����ڲ��������ܶ�
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(IISPHParticle* origin);
	//��������������Χ�Ķ�̬�������ӣ����ڲ����������ӵ��ܶȺ�ѹ��
	vector<DynamicRigidParticle*> ComputeNeighborDynamicParticle(IISPHParticle* origin);
	//����߽�������Χ�ı߽����ӣ����ڽ��б߽�����PSI����
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(StaticRigidParticle* origin);
	//���㶯̬����������Χ�ĸ������ӣ����ڽ��и�������PSI����
	vector<DynamicRigidParticle*> ComputeNeighborDynamicParticle(DynamicRigidParticle* origin);
	//���������ڵ���Χ�ı߽����ӣ����ھ�����ļ���
	vector<IISPHParticle*> ComputeNeighborFluidParticleDistance2(DistanceNode* origin);
	//����������Χ���ھӣ��������ӵ��ܶȣ������������ӵ�ѹǿ�Ѿ�����ɹ�
	double ComputeDensity(vector<IISPHParticle>& neighborList, vector<StaticRigidParticle>& boundaryNeighborList, IISPHParticle& origin);
    //����������Χ���ھӣ�����������ճ�����µļ��ٶ�
	Vector3f ComputeViscostyAcceleration(vector<IISPHParticle>& neighborList, IISPHParticle& origin);
	//�ڼ����ѹ��������£��������ӵ�ѹ�����ٶ�
	void ComputePressureAcceleration();
	//����������ھӱ�
	void UpdateFluidNeighborList();
	//��������ı߽��ھӱ�������Ҫ���ϸ���
	void UpdateFluidBoundaryNeighborList();
	//��������Ķ�̬���������ھӱ�������Ҫ���ϵĸ���
	void UpdateDynamicRigidNeighborsList();
	//�������ĸ���߽��ھӱ�ֻ��Ҫ��ʼ��һ�Σ�û��������͸������ײ������ֻ����Ψһ��һ����̬���壩
	void ComputeDynamicNeighborList();
	//����߽�ı߽��ھӱ�ֻ��Ҫ��ʼ������һ��
	void ComputeBoundaryNeighborList();
	//�������ӵ�λ�ý����������
	void MapParticleToGrid(); 
	//���߽�����ӳ�䵽������
	void MapBoundaryParticleToGrid(); 
	//����̬��������ӳ�䵽������
	void MapDynamicRigidParticleToGrid();
	//����ѹ������õ�ÿ�����ӵ�ѹ�����ٶ�
	void ComputePressureAccels();

	//���㵱ǰpos����Ӧ�ľ�����ֵ�����Ҳ���Һ��
	bool ProcessDistanceField(double radius = 0.025);


	//////////////////////////////////////////////////////////////////////////
	// ��������Ӧ��չ
	//////////////////////////////////////////////////////////////////////////
	//�����еķ���
	void ComputePressureAccels2();
	//��������Ӧ����
	void ComputeAdaptivityFactor();
	//�����������ӵ�����ľ���
	void ComputeDistance();
	//�������ӵ�mark�������ӽ��з���
	void ClassifyParticle();
	//�����������ӵ�������������
	void ComputeOptMass();
	//�������ӵ�ӡ�� �����������ӵ�O
	void ComputeO();
	//�����ں��ܶ�
	void BlendDensity();
	//�����ں��ٶ�
	void BlendVelocity();
	//�����������ӵ��ں�����
	void UpdateBlendFactor();

	//add
	//����ÿ�����ӵľ�����DistanceField
	void ComputePtoSandMopt();
	void ProcessMergeandSplit();
};