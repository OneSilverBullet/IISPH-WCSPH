#pragma once
#include <math.h>
#include "Vector3f.h"
#include "Vector3i.h"
#include "ShareData.h"


//////////////////////////////////////////////////////////////////////////
//FunctionKit类
//作用：包含一些散乱的工具函数，比如距离计算等等
//////////////////////////////////////////////////////////////////////////
class FunctionKit
{
public:
	//求两个点距离的函数
	static double DistanceComputation(Vector3f a, Vector3f b);
	//网格的哈希计算函数
	static int GridHash(int i, int j, int k, int cellNum); 
	//通过当前位置，计算映射当前网格索引
	static Vector3i PositionMapIndex(Vector3f parPosition, double cellLength, double offset);
	//通过粒子位置映射当前网格索引
	static int PositionMapHash(Vector3f parPosition, double cellLength, int cellNum, double offset); 
	//检测对应的数字是否为素数
	static bool CheckPrime(int n);
	//关于向量计算
	//返回两个向量的cos值
	static double GetVectorsCos(Vector3f a, Vector3f b);
	//依据发现向量返回对应的反射向量
	static Vector3f Reflection(Vector3f normal, Vector3f a);
	//这两个函数在计算状态方程的时候使用
	static double PressureConstantB(double density0, double soundVelocity, double r);
	//泰特方程计算
	static double TaitEquation(double density, double B,int gamma, double restDensity);
	//依据粒子的半径，计算得到当前粒子需要多少个方块
	static int CalculateGridCellNum(double radius);

};