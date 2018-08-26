#pragma once

#include "Vector3f.h"
#include "Vector3i.h"
#include "StaticRigidObject.h"
#include "ShareData.h"

//////////////////////////////////////////////////////////////////////////
//Boundary��
//���ã�����һЩ�߽����ԡ����ݱ߽����Բ����õ���ǰ�߽������λ������
//////////////////////////////////////////////////////////////////////////
class Boundary
{
public:
	//�߽������
	double boundaryOriginX = 0;
	double boundaryOriginY = 0;
	double boundaryOriginZ = 0; //���ڱ߽���Եı߽���ʼ��
	double boundaryWidth = 1;
	double boundaryLength = 2;
	double boundaryHeight = 1;  //���ڱ߽���Եı߽糤���
	//������ʷǷ��ڴ�Ľ������
	double  offset =  0.3f;               //������������Ĳ��ֵ�ƫ�ƣ�һ��Ϊ��������ĳ���
	double gridBoundaryOriginX = -0.5f;
	double gridBoundaryOriginY = -0.5f;
	double gridBoundaryOriginZ = -0.5f; //������������ı߽���ʼ��
	double gridBoundaryWidth = 3.0f;
	double gridBoundaryLength = 5.0f;
	double gridBoundaryHeight = 3.0f;   //������������ĳ����
	
	double distanceFieldOffset = 0.2f;
	double distanceFieldOriginX = -0.2f;
	double distanceFieldOriginY = -0.2f;
	double distanceFieldOriginZ = -0.2f;
	double distanceFieldWidth = 1.4f;
	double distanceFieldLength = 2.4f;
	double distanceFieldHeight = 1.4f;

	double pRad = 0.025f;

	//�Ƿ���ƾ�̬����
	bool staticRigidFlag = false;
	int boundaryNumber = 0;  //�߽�����ӵ�����

public:
	//������������
	Boundary();
	Boundary(Boundary& a);
	Boundary(double x, double y, double z, double w, double l, double h, double os);
	~Boundary();
	//�����ض���ֵ�����µ���������
	void GenerateNewBoundary(double x, double y, double z, double w, double l, double h, double os); //�����µ�����
	void GenerateNewBoundary(Vector3f position, Vector3f fluidWLH, double offs);
	Boundary& operator=(Boundary& a);
	//���ݲ����õ��ı߽����ݲ�������õ���ǰ�ı߽���������
	vector<Vector3f> initBoundaryData();
};