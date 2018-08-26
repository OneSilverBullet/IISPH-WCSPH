#pragma once
#include "WCSPHComputer.h"
#include "IISPHComputer.h"

//两种SPH模式
enum SPH_TYPE { IISPH, WCSPH };

class SPHManager
{
public:
	//两种不同类型的SPH 
	WCSPHComputation wcsphComputer;
	IISPHComputer iisphComputer;

	//可调节的参数
	double radius;  //半径
	double w;  
	double l;
	double h;       //水体的长宽高
	double ox;
	double oy;
	double oz;      //水体生成的初始位置
	SPH_TYPE sphtype = SPH_TYPE::IISPH;  //默认为WCSPH
	bool runFlag = false;  //开始运行的标志位 设置为false

	//工程标志位
	bool fluidAtiFluid = false;
	bool fluidAtiRigid = false;
	bool fluidSurface = false;

public:
	//初始化
	void Initialize();
	//计算当前的SPH模型
	void Compute();

	//标志位操控
	//开始计算
	void Run();
	//停止计算,模拟停滞
	void Stop();
	//改变当前SPH模式
	void SetSPHType(SPH_TYPE T);

	//参数控制接口
	void SetRadius(double rad);
	void SetFluidVolume(double x, double y, double z);
	void SetFluidOrigin(double ox, double oy, double oz);

	void SetRadius(double radius, double x, double y, double z, double ox, double oy, double oz);

	void ChangeToSigleDamBreak();
	void ChangeToDoubleDamBreak();
	void ChangeToStaticRigid();
	void ChangeToSurface();

};
