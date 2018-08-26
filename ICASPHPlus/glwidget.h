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
	bool isLock;//������ͼ
	bool isBackground;//��ʾ����
	float ballRadius; //���Ӱ뾶
	float timeStep;	  //ÿ֡ʱ�䲽��
	bool isShowParticle;
	bool isShowDisField;
	bool isShowSurface;
	bool isShowHashGrid;
	bool isShowBoundaryParticles;

	//����ģ��
	SPHManager sphManager;                 //SPH�Ĺ�����

	//�����Ƿ����ÿһ֡��ͼƬ
	void setOutputFramePicture(bool isOutput);
	void OutputFramePicture();
	void setShowBoundaryParticles(bool isShow);
	//��ʼ
	void StartRun();
	//ֹͣ
	void StopRun();
	
private slots:
	//����֡
	void timeFrame();
	//void updateFrame();

private:
	Ui::ICASPHPlusClass *ui;
	GLint iWinWidth;			//���ڳߴ�
	GLint iWinHeight;

	QTimer updateTimer;//������
	QTimer timer;//������
	QTime time;	//ʱ��
	int frame;	//��ǰ֡��,���ڼ�¼ÿ���֡��
	int sumFrame;//��֡
	int frameSec;
	double ratio;//���ű���
	/* �������� */
	int mouseX, mouseY;			//��ǰ���λ��
	bool isMouseDown;			//����Ƿ���
	float rotX = 0, rotY = 0;	//��ת�Ƕ�
	double zoom;
	double posX, posY;			//�ӵ�λ��
	int id;
	GLfloat lightDirection[3];		//���շ���
	GLfloat lightPosition[4];		//��Դλ��,xyz,dis
	bool isOutputFramePicture = true;


	//IISPHComputer       iisphcomputer;   //IISPH����ѧģ��
	/* ���� */
	vector<float> bottomBackgroud;
	vector<float> colorBackgroud;

	/* OpenGL */
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);
	void setLight();
	/* ����ģ�� */
	//���Ʊ������ذ�
	void drawBoundary(GLenum mode = GL_RENDER);
	void drawBoundaryParticles();
	void drawParticles();
	void drawDynamicRigidParticles();
	void drawFluid();
	void drawDistFieldGrid();
	void initWidget();

protected:
	/* ��Ӧ������openg */
	void keyPressEvent(QKeyEvent *e);
	void wheelEvent(QWheelEvent *e);
	void mouseMoveEvent(QMouseEvent *e);
	void mousePressEvent(QMouseEvent *e);
	void mouseReleaseEvent(QMouseEvent *e);
	void mouseDoubleClickEvent(QMouseEvent *e);


};

#endif // GLWIDGET_H
