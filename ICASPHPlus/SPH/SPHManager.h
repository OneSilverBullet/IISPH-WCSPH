#pragma once
#include "WCSPHComputer.h"
#include "IISPHComputer.h"

//����SPHģʽ
enum SPH_TYPE { IISPH, WCSPH };

class SPHManager
{
public:
	//���ֲ�ͬ���͵�SPH 
	WCSPHComputation wcsphComputer;
	IISPHComputer iisphComputer;

	//�ɵ��ڵĲ���
	double radius;  //�뾶
	double w;  
	double l;
	double h;       //ˮ��ĳ����
	double ox;
	double oy;
	double oz;      //ˮ�����ɵĳ�ʼλ��
	SPH_TYPE sphtype = SPH_TYPE::IISPH;  //Ĭ��ΪWCSPH
	bool runFlag = false;  //��ʼ���еı�־λ ����Ϊfalse

	//���̱�־λ
	bool fluidAtiFluid = false;
	bool fluidAtiRigid = false;
	bool fluidSurface = false;

public:
	//��ʼ��
	void Initialize();
	//���㵱ǰ��SPHģ��
	void Compute();

	//��־λ�ٿ�
	//��ʼ����
	void Run();
	//ֹͣ����,ģ��ͣ��
	void Stop();
	//�ı䵱ǰSPHģʽ
	void SetSPHType(SPH_TYPE T);

	//�������ƽӿ�
	void SetRadius(double rad);
	void SetFluidVolume(double x, double y, double z);
	void SetFluidOrigin(double ox, double oy, double oz);

	void SetRadius(double radius, double x, double y, double z, double ox, double oy, double oz);

	void ChangeToSigleDamBreak();
	void ChangeToDoubleDamBreak();
	void ChangeToStaticRigid();
	void ChangeToSurface();

};
