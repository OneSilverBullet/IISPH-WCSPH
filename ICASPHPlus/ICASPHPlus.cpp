#include "ICASPHPlus.h"

ICASPHPlus::ICASPHPlus(QWidget *parent)
	: QMainWindow(parent)
{
	ui.setupUi(this);
	init();
}
ICASPHPlus::~ICASPHPlus()
{

}

void ICASPHPlus::init()
{
	glwidget = new GLWidget(&ui);
	ui.BaseLayout->addWidget(glwidget);

}
void ICASPHPlus::onStart()
{
	ui.actionStart->setChecked(true);
	ui.actionStop->setChecked(false);
	glwidget->StartRun();

}

void ICASPHPlus::onStop()
{
	ui.actionStop->setChecked(true);
	ui.actionStart->setChecked(false);
	glwidget->StopRun();
}

void ICASPHPlus::isSavePic()
{
	if (ui.actionSavePic->isChecked() == true)
	{
		glwidget->setOutputFramePicture(true);
	}
	else
	{
		glwidget->setOutputFramePicture(false);
	}
	
}

void ICASPHPlus::onDrawDistanceField()
{
	glwidget->isShowDisField = true;
}

void ICASPHPlus::offDrawDistanceField()
{
	glwidget->isShowDisField = false;
}

void ICASPHPlus::changeIISPH()
{
	//先停止
	glwidget->StopRun();
	//设置为IISPH的渲染方式
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//改为单DamBreak模式
	glwidget->sphManager.ChangeToSigleDamBreak();
	//初始化
	glwidget->sphManager.Initialize();
	//刷新界面
	glwidget->update();
	ui.SPHTypeName->setText("IISPH DamBreak");
}

void ICASPHPlus::changeSingleBreak()
{
	//先停止
	glwidget->StopRun();
	//设置为IISPH的渲染方式
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//改为单DamBreak模式
	glwidget->sphManager.ChangeToSigleDamBreak();
	//初始化
	glwidget->sphManager.Initialize();
	//刷新界面
	glwidget->update();
	ui.SPHTypeName->setText("IISPH DamBreak");
}

void ICASPHPlus::changeDoubleBreak()
{
	//先停止
	glwidget->StopRun();
	//设置为IISPH的渲染方式
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//改为单DamBreak模式
	glwidget->sphManager.ChangeToDoubleDamBreak();
	//初始化
	glwidget->sphManager.Initialize();
	//刷新界面
	glwidget->update();
	ui.SPHTypeName->setText("IISPH DoubleDamBreak");
}

void ICASPHPlus::changeStaticRigid()
{
	//先停止
	glwidget->StopRun();
	//设置为IISPH的渲染方式
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//改为单DamBreak模式
	glwidget->sphManager.ChangeToStaticRigid();
	//初始化
	glwidget->sphManager.Initialize();
	//刷新界面
	glwidget->update();
	ui.SPHTypeName->setText("IISPH StaticRigid");
}

void ICASPHPlus::changeSurface()
{
	//先停止
	glwidget->StopRun();
	//设置为IISPH的渲染方式
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//改为单DamBreak模式
	glwidget->sphManager.ChangeToSurface();
	//初始化
	glwidget->sphManager.Initialize();
	//刷新界面
	glwidget->update();
	ui.SPHTypeName->setText("IISPH Surface");
}




void ICASPHPlus::changeWCSPH()
{
	//先停止
	glwidget->StopRun();
	//设置为IISPH的渲染方式
	glwidget->sphManager.SetSPHType(SPH_TYPE::WCSPH);
	//初始化
	glwidget->sphManager.Initialize();
	//刷新界面
	glwidget->update();
	ui.SPHTypeName->setText("WCSPH DamBreak");
}

void ICASPHPlus::isShowBoundaryPartilces()
{
	if (ui.actionshowBoundaryPar->isChecked() == true)
	{
		glwidget->setShowBoundaryParticles(true);
	}
	else
	{
		glwidget->setShowBoundaryParticles(false);
	}
}