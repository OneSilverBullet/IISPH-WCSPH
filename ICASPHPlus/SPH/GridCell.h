#pragma once
#include "Vector3f.h"
#include "Vector3i.h"
#include "WCSPHParticle.h"
#include "IISPHParticle.h"
#include "FunctionKit.h"
#include "ParticleObject.h"
#include "DynamicRigidParticle.h"
#include <vector>
using namespace std;

//////////////////////////////////////////////////////////////////////////
//GridCell��
//���ã��������ɵ�ǰ�ռ������ĵ�Ԫ�ࡣ
//////////////////////////////////////////////////////////////////////////
class GridCell
{
private:
	Vector3i indexVec;                            //����ǰ��x�� y, z���������������
	int hashIndex;                                //��ǰ����Ĺ�ϣֵ�����ڴ������񼯺��õ���
	double gridCellLength;                        //��ǰ����ĳ���
	vector<WCSPHParticle*> particleList;                //��ǰ�����������������ӣ�ʹ��һ��vector������б���
	vector<IISPHParticle*> IISPHparticleList;      //��ǰ��������������IISPH���ӣ�ʹ��һ��vector������б���
	vector<StaticRigidParticle*> boundaryList;     //��ǰ�����а����ı߽羲̬���ӣ���һ��vector������б���
	vector<DynamicRigidParticle*> dynamicParticleList; //��ǰ��������Ķ�̬�������ӣ���һ��vector������б���

public:
	//��ʼ������
	GridCell(){}
	GridCell(int i, int j, int k, int hashV) {
		Initialization(i, j, k, hashV);
	}
     GridCell(const GridCell& a)
	{
		indexVec = a.indexVec;
		hashIndex = a.hashIndex;
		gridCellLength = a.gridCellLength;
		ResetParticleList();
		for (int i=0; i<a.particleList.size(); i++)
		{
			PushParticle(a.particleList[i]);
		}

		for (int j = 0; j < a.boundaryList.size(); j++)
		{
			PushBoundaryParticle(a.boundaryList[j]);
		}
		for (int k = 0; k < a.IISPHparticleList.size(); k++)
		{
			PushParticle(a.IISPHparticleList[k]);
		}
		for (int p=0; p<a.dynamicParticleList.size(); p++)
		{
			PushDynamicParticle(a.dynamicParticleList[p]);
		}

	}
	~GridCell() {}

	//�����������
	void Initialization(int i, int j, int k, int hashV);       //�Ըú������г�ʼ��
	void ComputeCellLength();                                  //�����������ĳ��ȣ���Ϊ���ӵ�֧�Ű뾶
	void PushParticle(WCSPHParticle* a);                             //������װ�뵽������
	void PushParticle(IISPHParticle* a);                        //����װ�뵽������
	void PushBoundaryParticle(StaticRigidParticle* a);          //���߽�����װ��������
	void PushDynamicParticle(DynamicRigidParticle* a);          //����̬����װ�뵽������
	void ResetParticleList();                                  //���������б��������
	void ResetBoundaryParticleList();                          //�������������б��������
	void ResetDynamicParticleList();
	vector<WCSPHParticle*> GetParticleList();
	vector<IISPHParticle*> GetIISPHParticleList();              //�õ���ǰIISPH�������б�
	vector<StaticRigidParticle*> GetBoundaryParticleList();
	vector<DynamicRigidParticle*> GetDynamicParticleList();
	bool CheckGridCell(Vector3f parPosition);                  //���ݵ�ǰ���ӵ�λ���жϵ�ǰ�����Ƿ����������ķ�Χ��
	bool CheckGridCell(Vector3i parIndex);                     //���ݵ�ǰ���ӵ�λ�ñ�ţ���������ǰ�����Ƿ��ڵ�ǰ����ķ�Χ��
	bool CheckGridCell(WCSPHParticle particle);
	bool CheckGridCell(IISPHParticle* particle);
	bool CheckGridCell(DynamicRigidParticle* particle);
	//��ֵ������
	GridCell& operator=(const GridCell& a);

	//���Եķ��ʽӿ�
	void SetGridCellLength(double length);
	double GetGridCellLength();
	void SetIndexVec(Vector3i value);
	void SetIndexVec(int i, int j, int k);
	Vector3i GetIndexVec();
	void SetHashIndex(int value);
	int GetHashIndex();
};