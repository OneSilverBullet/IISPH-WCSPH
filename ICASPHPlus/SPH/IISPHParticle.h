#pragma once
#include "Vector3f.h"
#include "Boundary.h"
#include "FunctionKit.h"
#include "StaticRigidObject.h"
#include "ParticleObject.h"
#include "ShareData.h"
#include "DynamicRigidParticle.h"
#include <vector>

class IISPHParticle
{
public:
	//ͨ����������
	Vector3f position;        //���ӵ�λ��
	Vector3f velocity;        //���ӵ��ٶ�
	Vector3f acceleration;    //���ӵļ��ٶ�
	double particleRad;        //���Ӱ뾶���̶���С
	double particleSupportRad; //���ӵ�֧�Ű뾶����Ϊ���Ӱ뾶���ı�
	double particleMass;        //ÿ�����ӵ��������ɵ�ǰ���ܶ������Ӱ뾶����õ�
	double density;           //���ӵ��ܶ�
	Vector3f viscosity;       //���ӵ��������ٶ�
	Vector3f pressureAcceleration;        //���ӵ�ѹ�����ٶ�
	int hashIndex;            //��ǰ���������������ϣ��

	//IISPH���õ�����������
	double aii;
	Vector3f dii;
	Vector3f dij_pj;
	double density_adv;
	double pressure;
	double lastPressure;
	//�������ӵ�ǰ�Ĵ�С��̬�ļ�����Χ������
	int gridCellNum = 0;

	//add
	double distValue;
	double disToSurface;

	//Infinit Continous��չ�õ�����������
	double mopt;             //��ǰ���Ӿ������������
	double distance;         //�����������ľ���
	double mark;             //ÿ�����ӵı��
	double O;                //��ǰÿ�����Ӷ���O��������Ӧ��֧ͬ�Ű뾶����֮��Ľ���
	PARTICLE_TYPE pt;        //��ǰ���ӵ���������
	double lastDensity;      //��ǰ���ӵ���һ�ε��ܶ�ֵ
	Vector3f lastVelocity;     //��ǰ���ӵ���һ�ε��ٶ�ֵ
	double ��;               //��ǰ���ӵ��ܶȻ�ϡ��ٶȻ������
	bool NEW_UPDATE = false; //��ǰ���Ӹ��µ�״̬��һ�����ӷ��ѻ��ߺϲ�����ô��ֵ����Ϊtrue

	//�������Ե�ˢ������
	vector<IISPHParticle*> fluidNeighbors;                //������ָ�������ھӵ�ָ��
	vector<StaticRigidParticle*> boundaryNeighbors;       //������ָ��߽��ھӵ�ָ��
	vector<DynamicRigidParticle*> dynamicRigidNeighbors;  //������ָ��̬�����ָ��

public:
	//���캯��ֱ�ӵ��ó�ʼ������
	IISPHParticle() {
		//�����������Գ�ʼ��
		position.Zero();       
		velocity.Zero();
		acceleration.Zero();
		particleRad = 0.025f;
		particleSupportRad = 0.1f; 
		particleMass = 0.1f;       
		density = 0.0f;          
		viscosity.Zero();       
		pressureAcceleration.Zero();        
		hashIndex = 0;            

		//IISPH���Գ�ʼ��
		aii = 0.0;
		dii.Zero();
		dij_pj.Zero();
		density_adv = 0.0;
		pressure = 0.0;
		lastPressure = 0.0;
		pressureAcceleration.Zero();

		//��������Ӧ����չ���Գ�ʼ��
		distance = 0.0f;
		mark = 0.0f;
		mopt = 0.0f;
		distance = 0.0f;
		mark = 0;
		O = 0;
		pt = PARTICLE_TYPE::S;        
		lastDensity = 0;      
		lastVelocity.Zero();
		�� = 0;
		//�ٽ�����չ
		fluidNeighbors.resize(0);
		boundaryNeighbors.resize(0);
		dynamicRigidNeighbors.resize(0);
	}
	~IISPHParticle() {}
	void Initialization(Vector3f position, int v, double radius); //�����ӽ��г�ʼ������
	void InitializationPVS(Vector3f position, int v, double mass , double density);
	IISPHParticle& operator=(IISPHParticle& a);    //�����ӽ��и�ֵ��������
	void Intigration(Boundary bound, double timeStep);       //���ݼ��ٶȶ����ӵ��ٶȺ�λ�ý��л�����ֵ
	void ParticlePositionCorrection(Boundary bound);               //�Ե�ǰ���ӵ�λ�ý���У��

};