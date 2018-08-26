#include "SPHManager.h"

void SPHManager::Initialize()
{
	SetRadius(0.025, 0.5, 0.5, 0.5, 0.1, 0.1, 0.025);
	//依据当前的WCSPH状态进行
	switch (sphtype)
	{
	case IISPH:
		cout << "2" << endl;
		iisphComputer.Initialization();
		break;
	case WCSPH:
		wcsphComputer.Initialization();
		break;
	default:
		break;
	}
}

void SPHManager::Compute()
{
	//如果当前runFlag设置为true，那么就进行计算
	if (runFlag == true)
	{
		//计算不同的SPH的情况
		switch (sphtype)
		{
		case IISPH:
			iisphComputer.Computation();
			break;
		case WCSPH:
			wcsphComputer.Computation();
			break;
		default:
			break;
		}
	}
}

void SPHManager::Run()
{
	runFlag = true;
}

void SPHManager::Stop()
{
	runFlag = false;
}

void SPHManager::SetSPHType(SPH_TYPE T)
{
	//修正当前SPH的类型
	sphtype = T;
	//停止当前的SPH计算
	Stop();
}

void SPHManager::SetRadius(double rad)
{
	radius = rad;
}

void SPHManager::SetFluidVolume(double x, double y, double z)
{
	w = x;
	l = y;
	h = z;
}

void SPHManager::SetFluidOrigin(double x, double y, double z)
{
	ox = x;
	oy = y;
	oz = z;
}

void SPHManager::SetRadius(double radius, double x, double y, double z, double ox, double oy, double oz)
{
	//设置IISPH的各种粒子半径
	//更新当前的半径
	iisphComputer.fluidModel.particleInitialPad = radius;
	//更新当前水体的长宽高
	iisphComputer.fluidModel.IISPH_FluidWitdth = x;
	iisphComputer.fluidModel.IISPH_FluidLength = y;
	iisphComputer.fluidModel.IISPH_FluidHeight = z;
	iisphComputer.fluidModel.IISPH_OriginX = ox;
	iisphComputer.fluidModel.IISPH_OriginY = oy;
	iisphComputer.fluidModel.IISPH_OriginZ = oz;
	//更新当前的boundary
	iisphComputer.boundary.pRad = radius;

	//设置WCSPH的各种粒子半径
	wcsphComputer.fluidModel.particleRad = radius;
	//更新当前水体的长宽高
	wcsphComputer.fluidModel.fluidWitdth = x;
	wcsphComputer.fluidModel.fluidLength = y;
	wcsphComputer.fluidModel.fluidHeight = z;
	wcsphComputer.fluidModel.originX = ox;
	wcsphComputer.fluidModel.originY = oy;
	wcsphComputer.fluidModel.originZ = oz;
	//更新当前的boundary
	wcsphComputer.boundary.pRad = radius;
}

void SPHManager::ChangeToSigleDamBreak()
{
    fluidAtiFluid = false;
	fluidAtiRigid = false;
	fluidSurface = false;
	iisphComputer.fluidAtiFluid = fluidAtiFluid;
	iisphComputer.fluidAtiRigid = fluidAtiRigid;
	iisphComputer.fluidSurface = fluidSurface;
}

void SPHManager::ChangeToDoubleDamBreak()
{
	fluidAtiFluid = true;
	fluidAtiRigid = false;
	fluidSurface = false;
	iisphComputer.fluidAtiFluid = fluidAtiFluid;
	iisphComputer.fluidAtiRigid = fluidAtiRigid;
	iisphComputer.fluidSurface = fluidSurface;
}

void SPHManager::ChangeToStaticRigid()
{
	fluidAtiFluid = false;
	fluidAtiRigid = true;
	fluidSurface = false;
	iisphComputer.fluidAtiFluid = fluidAtiFluid;
	iisphComputer.fluidAtiRigid = fluidAtiRigid;
	iisphComputer.fluidSurface = fluidSurface;
}

void SPHManager::ChangeToSurface()
{
	fluidAtiFluid = false;
	fluidAtiRigid = false;
	fluidSurface = true;
	iisphComputer.fluidAtiFluid = fluidAtiFluid;
	iisphComputer.fluidAtiRigid = fluidAtiRigid;
	iisphComputer.fluidSurface = fluidSurface;
}
