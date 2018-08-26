#pragma once
#include "GridCell.h"
#include "FunctionKit.h"
#include "Boundary.h"
#include <vector>
using namespace std;

//////////////////////////////////////////////////////////////////////////
//HashGridList��
//���ã��������ɵ�ǰ�ռ��������࣬��Ҫ�ڳ�ʼ������ʹ�á�
//�ں������㵱�У����ϸ���ÿ�������е��������ӡ�
//////////////////////////////////////////////////////////////////////////
class HashGridList
{
public:
	//ע����������񳤶ȣ���������е�����
	double cellLength = 0.1f;
	int hashCellNum = 0.0f;                //���ݵ�ǰ����ĸ�����ȷ����ѵĹ�ϣ�����
	int cellNum;                           //��ǰ����һ���ж��ٸ�
	double hashThred = 0.5;                //��ǰ�����װ������
	Vector3i gridWLH;                      //��ǰ�����ڳ�����ϵĸ���
	Boundary boundary;                     //���ݵ�ǰ�ı߽�����������������      
	vector<vector<GridCell>> gridCellList; //����Ĺ�ϣ���������񶼴�������


public:
	//Ĭ�Ͻ��г�ʼ��
	HashGridList() {}
	~HashGridList(){}
	void Intialization();                                               //������������ĳ�ʼ����������ʼ��һ��
	void Intialization(double cellL, double thred, Boundary bound);     //������������
	void ComputeGridWLH (double cellL, Boundary bound);                 //�����������ĳ��ȣ���Ϊ���ӵ�֧�Ű뾶
	void ComputeCellNum();                                              //�������������
	void InitializeHashCellNun();                                       //����ϣ��Ĵ�С����Ϊ��������֮��ĵ�һ���������г�ʼ��
	void EnsureTheGrid();                                               //���ݵ�ǰ���ȷ����ϣ��Ĵ�С���������ɶ�Ӧ�ı��
	void GenerateGridCell();                                            //�Կռ���л��֣�������������
	void Rehash();                                                      //����ѡ���ϣ��Ĵ�С
	bool CheckThred();                                                  //��⵱ǰ��ϣ���С�Ƿ���һ���ñ�
	void PushParticle(WCSPHParticle* a);                                      //��һ������ͨ������Listװ���Ӧ��������
	void PushParticle(StaticRigidParticle* a);                           //��һ���߽�����ͨ������Listװ���Ӧ��������
	void PushParticle(IISPHParticle* a);
	void PushParticle(DynamicRigidParticle* a);

	void ClearParticle();                                               //����ǰ�����е����ӽ������
	void ClearBoundaryParticle();                                       //����̬�����еľ�̬���ӽ������
    //����ӿڡ�
	double GetCellLength();
	int GetCellNum();
	int GetHashCellNum();

	//�����ѯ�������������ӡ�λ�õ����Բ���Ŀ������
	//���ݲ�ͬ����Ϣ�õ���ǰ�б��ж�Ӧ������
	vector<GridCell> GetGridCellVector(int hashIndex);
	//��������λ��index���в�ѯ���ڹ�ϣֵ������
	vector<GridCell> GetGridCellVector(Vector3i positionIndex);
	//��������λ�ò�ѯ���ڹ�ϣֵ������
	vector<GridCell> GetGridCellVector(Vector3f position);
	//���ݲ�ѯ�����������о��������ѯ
	GridCell GetGridCellWithSecVec(vector<GridCell>& gridVec, Vector3f parPosition);
	GridCell GetGridCellWithSecVec(vector<GridCell>& gridVec, WCSPHParticle a);
	//ֱ�Ӵӵ�ǰ�����������н��в�ѯ
	GridCell GetGridCell(WCSPHParticle a);
	GridCell GetGridCell(Vector3f parPosition);
};
