#pragma once
#include "HashGridList.h"
#include "GridCell.h"
#include "Vector3f.h"
#include "Vector3i.h"
#include "ShareData.h"
#include "FunctionKit.h"
#include "SPHKernel.h"
#include "WCSPHFluidObject.h"
#include "Boundary.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
using namespace std;


//////////////////////////////////////////////////////////////////////////
//SPHComputation��
//���ã�������FluidModel�Լ�HashList������Ķ���ʵ����SPH�ļ��㡣
//////////////////////////////////////////////////////////////////////////
class WCSPHComputation
{
public:
	HashGridList hashGridList; //��ϣ�����
	WCSPHFluidObject fluidModel;    //�������ģ��
	Boundary  boundary;        //����ı߽�ģ��
							   //�������嵱�е�ѹǿ
	double stiffnesss = 50000; //���嵱�еĸ���
	int gamma = 7;             //���嵱�е�gammaֵ
	double viscosityConstant = 0.03; //����������Ҫ�ĳ���
	double gravity = -9.81;    //���������

public:
	//��ʼ��
	void Initialization(); //��SPH����֮ǰ����һ�Σ������ӵ����ԡ�λ�ý��г�ʼ��
	void MapParticleToGrid(); //�������ӵ�λ�ý����������
	void MapBoundaryParticleToGrid(); //���߽�����ӳ�䵽������
									  //��֡
	void Frame();   //�������������֡���м���
					//���㲿�ֵĺ���
	void Computation();  //SPH���㣬�ú�����װ���з���
						 //����������Χ�ھ�����
	vector<WCSPHParticle*> ComputeNeighborParticle(WCSPHParticle* origin);
	//����������Χ�ı߽����ӣ����ڲ��������ܶ�
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(WCSPHParticle* origin);
	//����߽�������Χ�ı߽����ӣ����ڽ��б߽�����PSI����
	vector<StaticRigidParticle*> ComputeNeighborBoundaryParticle(StaticRigidParticle* origin);
	//����������Χ���ھӣ��������ӵ��ܶȣ������������ӵ�ѹǿ�Ѿ�����ɹ�
	double ComputeDensity(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin);
	//����������Χ���ھӣ�����������ѹ���µļ��ٶ�
	Vector3f ComputePressureAcceleration(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin);
	//����������Χ���ھӣ�����������ճ�����µļ��ٶ�
	Vector3f ComputeViscostyAcceleration(vector<WCSPHParticle*>& neighborList, vector<StaticRigidParticle*>& boundaryNeighborList, WCSPHParticle* origin);
	//�������ӵļ��ٶ�֮��
	Vector3f ComputeAcceleration(Vector3f pressureA, Vector3f viscosityA);
	//�����ӵ��ܶȽ�������ֹ�ܶȺ͵�ǰ�ܶȵ����ֵ
	double CorrectionParticleDensity(WCSPHParticle* a);
	//���㵱ǰ���ӵ�̩�ط��̵õ����ӵ�ѹǿ
	double CalculateParticleTaite(double a);
	//����������ھӱ�
	void UpdateFluidNeighborList();
	//��������ı߽��ھӱ�������Ҫ���ϸ���
	void UpdateFluidBoundaryNeighborList();
	//����߽�ı߽��ھӱ�ֻ��Ҫ��ʼ������һ��
	void ComputeBoundaryNeighborList();
};
