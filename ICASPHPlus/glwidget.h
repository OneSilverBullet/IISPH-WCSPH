#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QOpenGLWidget>
#include <QtOpenGL/QtOpenGL>
#include <QWheelEvent>
#include <qimage.h>
#include <qtimer.h>
#include <fstream>
#include <direct.h> 
#include <cstdlib>
#include "ui_ICASPHPlus.h"
#include <GL\freeglut.h>
#include <time.h>
#include <iostream>
#include "../ICASPHPlus/SPH/IISPHComputer.h"
#include "../ICASPHPlus/SPH/WCSPHFluidObject.h"
#include "../ICASPHPlus/SPH/SPHManager.h"

using namespace std;

class GLWidget : public QOpenGLWidget
{
	Q_OBJECT

public:
	GLWidget(Ui::ICASPHPlusClass *ui = NULL, int id = 0, QWidget* parent = 0);
	~GLWidget();
	bool isLock;//锁定视图
	bool isBackground;//显示背景
	float ballRadius; //粒子半径
	float timeStep;	  //每帧时间步长
	bool isShowParticle;
	bool isShowDisField;
	bool isShowSurface;
	bool isShowHashGrid;
	bool isShowBoundaryParticles;

	//计算模型
	SPHManager sphManager;                 //SPH的管理者

	//设置是否输出每一帧的图片
	void setOutputFramePicture(bool isOutput);
	void OutputFramePicture();
	void setShowBoundaryParticles(bool isShow);
	//开始
	void StartRun();
	//停止
	void StopRun();
	
private slots:
	//动画帧
	void timeFrame();
	//void updateFrame();

private:
	Ui::ICASPHPlusClass *ui;
	GLint iWinWidth;			//窗口尺寸
	GLint iWinHeight;

	QTimer updateTimer;//计数器
	QTimer timer;//计数器
	QTime time;	//时钟
	int frame;	//当前帧数,用于记录每秒的帧率
	int sumFrame;//总帧
	int frameSec;
	double ratio;//缩放比例
	/* 环境变量 */
	int mouseX, mouseY;			//当前鼠标位置
	bool isMouseDown;			//鼠标是否按下
	float rotX = 0, rotY = 0;	//旋转角度
	double zoom;
	double posX, posY;			//视点位置
	int id;
	GLfloat lightDirection[3];		//光照方向
	GLfloat lightPosition[4];		//光源位置,xyz,dis
	bool isOutputFramePicture = true;


	//IISPHComputer       iisphcomputer;   //IISPH的数学模型
	/* 背景 */
	vector<float> bottomBackgroud;
	vector<float> colorBackgroud;

	/* OpenGL */
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);
	void setLight();
	/* 绘制模型 */
	//绘制背景，地板
	void drawBoundary(GLenum mode = GL_RENDER);
	void drawBoundaryParticles();
	void drawParticles();
	void drawDynamicRigidParticles();
	void drawFluid();
	void drawDistFieldGrid();
	void initWidget();

protected:
	/* 响应函数，openg */
	void keyPressEvent(QKeyEvent *e);
	void wheelEvent(QWheelEvent *e);
	void mouseMoveEvent(QMouseEvent *e);
	void mousePressEvent(QMouseEvent *e);
	void mouseReleaseEvent(QMouseEvent *e);
	void mouseDoubleClickEvent(QMouseEvent *e);


};

#endif // GLWIDGET_H
