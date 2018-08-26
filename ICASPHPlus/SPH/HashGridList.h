#pragma once
#include "GridCell.h"
#include "FunctionKit.h"
#include "Boundary.h"
#include <vector>
using namespace std;

//////////////////////////////////////////////////////////////////////////
//HashGridList类
//作用：用于生成当前空间的网格的类，主要在初始化当中使用。
//在后续计算当中，不断更新每个网格当中的流体粒子。
//////////////////////////////////////////////////////////////////////////
class HashGridList
{
public:
	//注：这里的网格长度，是针对所有的网格
	double cellLength = 0.1f;
	int hashCellNum = 0.0f;                //依据当前网格的个数，确定最佳的哈希表个数
	int cellNum;                           //当前网格一共有多少个
	double hashThred = 0.5;                //当前网格的装填因子
	Vector3i gridWLH;                      //当前网格在长宽高上的个数
	Boundary boundary;                     //依据当前的边界来对整个进行设置      
	vector<vector<GridCell>> gridCellList; //网格的哈希表，所有网格都存入其中


public:
	//默认进行初始化
	HashGridList() {}
	~HashGridList(){}
	void Intialization();                                               //对于所有网格的初始化，仅仅初始化一次
	void Intialization(double cellL, double thred, Boundary bound);     //带参数的网格
	void ComputeGridWLH (double cellL, Boundary bound);                 //这里计算网格的长度，即为粒子的支撑半径
	void ComputeCellNum();                                              //计算网格的数量
	void InitializeHashCellNun();                                       //将哈希表的大小设置为网格数量之后的第一个素数进行初始化
	void EnsureTheGrid();                                               //根据当前情况确定哈希表的大小，并且生成对应的表格
	void GenerateGridCell();                                            //对空间进行划分，进行网格生成
	void Rehash();                                                      //重新选择哈希表的大小
	bool CheckThred();                                                  //检测当前哈希表大小是否是一个好表
	void PushParticle(WCSPHParticle* a);                                      //将一个粒子通过网格List装入对应的网格中
	void PushParticle(StaticRigidParticle* a);                           //将一个边界粒子通过网格List装入对应的网格当中
	void PushParticle(IISPHParticle* a);
	void PushParticle(DynamicRigidParticle* a);

	void ClearParticle();                                               //将当前网格当中的粒子进行清除
	void ClearBoundaryParticle();                                       //将静态网格当中的静态粒子进行清除
    //对外接口。
	double GetCellLength();
	int GetCellNum();
	int GetHashCellNum();

	//网格查询函数，依据粒子、位置等属性查找目标网格
	//根据不同的信息得到当前列表当中对应的网格
	vector<GridCell> GetGridCellVector(int hashIndex);
	//依据粒子位置index进行查询所在哈希值的容器
	vector<GridCell> GetGridCellVector(Vector3i positionIndex);
	//依据粒子位置查询所在哈希值的容器
	vector<GridCell> GetGridCellVector(Vector3f position);
	//依据查询过的容器进行具体网格查询
	GridCell GetGridCellWithSecVec(vector<GridCell>& gridVec, Vector3f parPosition);
	GridCell GetGridCellWithSecVec(vector<GridCell>& gridVec, WCSPHParticle a);
	//直接从当前网格容器当中进行查询
	GridCell GetGridCell(WCSPHParticle a);
	GridCell GetGridCell(Vector3f parPosition);
};
