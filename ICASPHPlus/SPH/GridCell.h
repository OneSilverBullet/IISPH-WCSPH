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
//GridCell类
//作用：用于生成当前空间的网格的单元类。
//////////////////////////////////////////////////////////////////////////
class GridCell
{
private:
	Vector3i indexVec;                            //网格当前从x， y, z三个方向的索引。
	int hashIndex;                                //当前网格的哈希值，用于存入网格集合用到。
	double gridCellLength;                        //当前网格的长度
	vector<WCSPHParticle*> particleList;                //当前网格当中所包含的粒子，使用一个vector对其进行保存
	vector<IISPHParticle*> IISPHparticleList;      //当前网格当中所包含的IISPH粒子，使用一个vector对其进行保存
	vector<StaticRigidParticle*> boundaryList;     //当前网格当中包含的边界静态粒子，用一个vector对其进行保存
	vector<DynamicRigidParticle*> dynamicParticleList; //当前网络包含的动态物体粒子，用一个vector对其进行保存

public:
	//初始化函数
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

	//网格操作函数
	void Initialization(int i, int j, int k, int hashV);       //对该函数进行初始化
	void ComputeCellLength();                                  //这里计算网格的长度，即为粒子的支撑半径
	void PushParticle(WCSPHParticle* a);                             //将粒子装入到网格当中
	void PushParticle(IISPHParticle* a);                        //奖励装入到网格当中
	void PushBoundaryParticle(StaticRigidParticle* a);          //将边界粒子装入网格当中
	void PushDynamicParticle(DynamicRigidParticle* a);          //将动态粒子装入到网格当中
	void ResetParticleList();                                  //重置粒子列表，清空网格
	void ResetBoundaryParticleList();                          //重置网格粒子列表，清空网格
	void ResetDynamicParticleList();
	vector<WCSPHParticle*> GetParticleList();
	vector<IISPHParticle*> GetIISPHParticleList();              //得到当前IISPH的粒子列表
	vector<StaticRigidParticle*> GetBoundaryParticleList();
	vector<DynamicRigidParticle*> GetDynamicParticleList();
	bool CheckGridCell(Vector3f parPosition);                  //根据当前粒子的位置判断当前粒子是否处于这个网格的范围内
	bool CheckGridCell(Vector3i parIndex);                     //根据当前粒子的位置编号，来决定当前粒子是否处于当前网格的范围内
	bool CheckGridCell(WCSPHParticle particle);
	bool CheckGridCell(IISPHParticle* particle);
	bool CheckGridCell(DynamicRigidParticle* particle);
	//赋值的重载
	GridCell& operator=(const GridCell& a);

	//属性的访问接口
	void SetGridCellLength(double length);
	double GetGridCellLength();
	void SetIndexVec(Vector3i value);
	void SetIndexVec(int i, int j, int k);
	Vector3i GetIndexVec();
	void SetHashIndex(int value);
	int GetHashIndex();
};