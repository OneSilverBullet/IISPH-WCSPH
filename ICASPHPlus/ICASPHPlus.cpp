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
	//��ֹͣ
	glwidget->StopRun();
	//����ΪIISPH����Ⱦ��ʽ
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//��Ϊ��DamBreakģʽ
	glwidget->sphManager.ChangeToSigleDamBreak();
	//��ʼ��
	glwidget->sphManager.Initialize();
	//ˢ�½���
	glwidget->update();
	ui.SPHTypeName->setText("IISPH DamBreak");
}

void ICASPHPlus::changeSingleBreak()
{
	//��ֹͣ
	glwidget->StopRun();
	//����ΪIISPH����Ⱦ��ʽ
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//��Ϊ��DamBreakģʽ
	glwidget->sphManager.ChangeToSigleDamBreak();
	//��ʼ��
	glwidget->sphManager.Initialize();
	//ˢ�½���
	glwidget->update();
	ui.SPHTypeName->setText("IISPH DamBreak");
}

void ICASPHPlus::changeDoubleBreak()
{
	//��ֹͣ
	glwidget->StopRun();
	//����ΪIISPH����Ⱦ��ʽ
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//��Ϊ��DamBreakģʽ
	glwidget->sphManager.ChangeToDoubleDamBreak();
	//��ʼ��
	glwidget->sphManager.Initialize();
	//ˢ�½���
	glwidget->update();
	ui.SPHTypeName->setText("IISPH DoubleDamBreak");
}

void ICASPHPlus::changeStaticRigid()
{
	//��ֹͣ
	glwidget->StopRun();
	//����ΪIISPH����Ⱦ��ʽ
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//��Ϊ��DamBreakģʽ
	glwidget->sphManager.ChangeToStaticRigid();
	//��ʼ��
	glwidget->sphManager.Initialize();
	//ˢ�½���
	glwidget->update();
	ui.SPHTypeName->setText("IISPH StaticRigid");
}

void ICASPHPlus::changeSurface()
{
	//��ֹͣ
	glwidget->StopRun();
	//����ΪIISPH����Ⱦ��ʽ
	glwidget->sphManager.SetSPHType(SPH_TYPE::IISPH);
	//��Ϊ��DamBreakģʽ
	glwidget->sphManager.ChangeToSurface();
	//��ʼ��
	glwidget->sphManager.Initialize();
	//ˢ�½���
	glwidget->update();
	ui.SPHTypeName->setText("IISPH Surface");
}




void ICASPHPlus::changeWCSPH()
{
	//��ֹͣ
	glwidget->StopRun();
	//����ΪIISPH����Ⱦ��ʽ
	glwidget->sphManager.SetSPHType(SPH_TYPE::WCSPH);
	//��ʼ��
	glwidget->sphManager.Initialize();
	//ˢ�½���
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