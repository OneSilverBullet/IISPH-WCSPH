#pragma once
#include "Vector3f.h"
#include "ShareData.h"
#include "FunctionKit.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////////
//Particle��
//���ã�����ģ������Ļ������ӣ������������ӵı�������
//////////////////////////////////////////////////////////////////////////
class WCSPHParticle
{
public:
	Vector3f position;        //���ӵ�λ��
	Vector3f velocity;        //���ӵ��ٶ�
	Vector3f acceleration;    //���ӵļ��ٶ�
	double P;                 //���ӵ�ѹǿ
	double density;           //���ӵ��ܶ�
	Vector3f viscosity;       //���ӵ�����
	Vector3f pressure;        //���ӵ�ѹ��
	int hashIndex;            //��ǰ���������������ϣ��

public:
	//���캯��ֱ�ӵ��ó�ʼ������
	WCSPHParticle() {}
	~WCSPHParticle() {}
	void Initialization(Vector3f position, int v); //�����ӽ��г�ʼ������

												   //�����ӽ��и�ֵ
	WCSPHParticle& operator=(WCSPHParticle& a);
};