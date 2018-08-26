#include "SPHManager.h"

void SPHManager::Initialize()
{
	SetRadius(0.025, 0.5, 0.5, 0.5, 0.1, 0.1, 0.025);
	//���ݵ�ǰ��WCSPH״̬����
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
	//�����ǰrunFlag����Ϊtrue����ô�ͽ��м���
	if (runFlag == true)
	{
		//���㲻ͬ��SPH�����
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
	//������ǰSPH������
	sphtype = T;
	//ֹͣ��ǰ��SPH����
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
	//����IISPH�ĸ������Ӱ뾶
	//���µ�ǰ�İ뾶
	iisphComputer.fluidModel.particleInitialPad = radius;
	//���µ�ǰˮ��ĳ����
	iisphComputer.fluidModel.IISPH_FluidWitdth = x;
	iisphComputer.fluidModel.IISPH_FluidLength = y;
	iisphComputer.fluidModel.IISPH_FluidHeight = z;
	iisphComputer.fluidModel.IISPH_OriginX = ox;
	iisphComputer.fluidModel.IISPH_OriginY = oy;
	iisphComputer.fluidModel.IISPH_OriginZ = oz;
	//���µ�ǰ��boundary
	iisphComputer.boundary.pRad = radius;

	//����WCSPH�ĸ������Ӱ뾶
	wcsphComputer.fluidModel.particleRad = radius;
	//���µ�ǰˮ��ĳ����
	wcsphComputer.fluidModel.fluidWitdth = x;
	wcsphComputer.fluidModel.fluidLength = y;
	wcsphComputer.fluidModel.fluidHeight = z;
	wcsphComputer.fluidModel.originX = ox;
	wcsphComputer.fluidModel.originY = oy;
	wcsphComputer.fluidModel.originZ = oz;
	//���µ�ǰ��boundary
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
