#include "glwidget.h"

GLWidget::GLWidget(Ui::ICASPHPlusClass *ui, int id, QWidget* parent)
	: QOpenGLWidget(parent)
{
	this->id = id;
	this->ui = ui;
	initializeGL();
	initWidget();

}

GLWidget::~GLWidget()
{

}

//**********************************初始化*********************************************
void GLWidget::initWidget()
{
	//调用当前初始化
	cout << "1" << endl;
	sphManager.Initialize();

	isMouseDown = false;			//鼠标是否按下
	frame = 0;
	sumFrame = 0;
	frameSec = 0;
	rotX = -90;//初始化旋转角度
	rotY = 0;
	zoom = -3.0;
	posX = 0;
	posY = 0;
	//isLock = false;
	//isBackground = false;

	//isShowParticle = Param::ctrlParam[IS_PARTICLE];
	//isShowDisField = Param::ctrlParam[IS_DIS_FIELD];
	//isShowSurface = Param::ctrlParam[IS_SURFACE];
	//isShowHashGrid = Param::ctrlParam[IS_HASH_GRID];
	//isOutputFrame = Param::ctrlParam[IS_OUTPUT_FRAME];

	//设置参数
	//ratio = 3.6 / pcisph.wBottle;//缩放比例
	ratio = 0.5;
	//ballRadius = 0.01;	//粒子半径
	timeStep = 10;		//帧间隔,ms

	lightPosition[0] = 0.0f;	//光源位置,xyz,dis
	lightPosition[1] = 1.5f;
	lightPosition[2] = 0.0f;
	lightPosition[3] = 0.0f;
	lightDirection[0] = 0.0f;	//光照方向
	lightDirection[0] = -1.0f;
	lightDirection[0] = 0.0f;

	//设置并获取背景
	//double bgy = pcisph.yBottle * ratio - 0.1;//背景y坐标
	//setBackground(30, 30, bgy, 1.0);
	setFocusPolicy(Qt::StrongFocus);

	/*****开始模拟*****/
	
	//timer.start(timeStep);
	timer.stop();

	//开始时钟，以及播放一帧
	time.start();
	//timeFrame();

	QObject::connect(&timer, SIGNAL(timeout()), this, SLOT(timeFrame()));
	//qDebug() << QString::fromLocal8Bit("当前线程id：") << QThread::currentThread();

}

void GLWidget::initializeGL()
{
	glClearColor(0, 0, 0, 1);//背景色

	glEnable(GL_DEPTH_TEST);//深度测试，深度值有效
	glEnable(GL_BLEND);		//颜色混合
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);//启用颜色追踪
	glShadeModel(GL_SMOOTH);
	glClearDepth(1.0f);		//设置深度范围为最大

							//设置光照参数
	setLight();
}

void GLWidget::setLight()
{
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	//GLfloat lightPosition[] = { 0.0, 1.2, -1.0, 1.0 };		//光源位置,xyz,dis
	GLfloat lightAmbient[] = { 0.6, 0.6, 0.6, 1.0 };		//全局光属性,光强度
	GLfloat lightDiffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	//GLfloat lightDirection[] = { 0.0f, -1.0f, 0.0f };	//光照方向
	GLfloat lightExponent[] = { 1.0f };					//聚光程度
	GLfloat lightCutoff[] = { 60.0f };					//光源夹角

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);//指定需要颜色追踪的材料属性
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);//指定第0号光源的位置 
													 //反射
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse); //漫反射后
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);//镜面反射后
													 //聚光
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, lightDirection);
	//glLightfv(GL_LIGHT0, GL_SPOT_EXPONENT, lightExponent);
	//glLightfv(GL_LIGHT0, GL_SPOT_CUTOFF, lightCutoff);

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lightAmbient);
}

void GLWidget::resizeGL(int width, int height)
{
	iWinWidth = (GLint)width;
	iWinHeight = (GLint)height;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	//投影方式：透视投影
	gluPerspective(60, (GLfloat)iWinWidth / (GLfloat)iWinHeight, 0.1, 100);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
}

//**********************************事件*********************************************
void GLWidget::mouseMoveEvent(QMouseEvent *e)
{
	if (isMouseDown)
	{
		int deltX, deltY;
		// 计算鼠标移动距离，旋转图像
		deltX = e->x() - mouseX;
		deltY = e->y() - mouseY;

		//if (false == isLock)
		{
			// 旋转角  
			rotX += deltX / 2;
			rotY += deltY / 2;
			// 旋转角不超过360度 
			rotX = fmodf(rotX, 360.0);
			rotY = fmodf(rotY, 360.0);
		}
		//else
		{
			//pcisph.xBottle += double(deltX) / this->width() / ratio * 4;
			//pcisph.yBottle -= double(deltY) / this->height() / ratio * 4;
		}

		//更新当前鼠标位置，使图像实时旋转
		mouseX = e->x();
		mouseY = e->y();
		update();
	}
}
void GLWidget::mousePressEvent(QMouseEvent *e)
{
	mouseX = e->x();
	mouseY = e->y();

	isMouseDown = true;


}
void GLWidget::mouseReleaseEvent(QMouseEvent *e)
{
	isMouseDown = false;
}
void GLWidget::mouseDoubleClickEvent(QMouseEvent *e)
{

}
void GLWidget::wheelEvent(QWheelEvent *e)
{
	int numDegrees = e->delta() / 8;//滚动的角度，*8就是鼠标滚动的距离
	double numSteps = numDegrees / 120.0;//滚动的步数，*15就是鼠标滚动的角度
	zoom += numSteps / 2;
	e->accept();      //接收该事件
	update();
}
void GLWidget::keyPressEvent(QKeyEvent *e)
{
	double dstep = 0.1;
	switch (e->key())
	{
	case Qt::Key_W:
		posY -= dstep;
		//zoom += (GLfloat)0.05;
		break;
	case Qt::Key_S:
		posY += dstep;
		//zoom -= (GLfloat)0.05;
		break;
	case Qt::Key_A:
		posX += dstep;
		//rotX -= (GLfloat)5;
		break;
	case Qt::Key_D:
		posX -= dstep;
		//rotX += (GLfloat)5;
		break;
	}
	update();
}
void GLWidget::timeFrame()
{
	//cout << "timeFrame" << endl;
	static int s = 0, e = 0;
	s = clock();
	sphManager.Compute();
	//SPHComputer.Computation();
	//刷新当前的图像
	this->update();

	//绘制当前帧率
	frame++;
	sumFrame++;
	frameSec = clock() - s;
	if (1000 <= time.elapsed())
	{
		if (sphManager.sphtype==SPH_TYPE::IISPH)
		{
			ui->statusBar->showMessage(QString("FPS : %1 / s\t  Frame: %2 ms\t  ParticleNum: %3\t  Frames:%4").arg(frame).arg(frameSec)
				.arg(sphManager.iisphComputer.fluidModel.particleList.size()).arg(sumFrame));
			time.restart();
			frame = 0;
		}
		else
		{
			ui->statusBar->showMessage(QString("FPS : %1 / s\t  Frame: %2 ms\t  ParticleNum: %3\t  Frames:%4").arg(frame).arg(frameSec)
				.arg(sphManager.wcsphComputer.fluidModel.particleList.size()).arg(sumFrame));
			time.restart();
			frame = 0;
		}
	}
	OutputFramePicture();
	/*if (isOutputFramePicture)
	{
		OutputFramePicture();
	}*/

	e = clock();
}

//**********************************绘制*********************************************
void GLWidget::paintGL()
{
	static int sd = 0, ed = 0;
	//sd = clock();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();

	glTranslatef(posX, posY, zoom);
	glRotatef(rotX, 0.0, 1.0, 0.0);
	glRotatef(rotY, 1.0, 0.0, 0.0);
	//glRotatef(90, 1.0, 0.0, 0.0);

	//移动光源同步移动
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);//指定第0号光源的位置 
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, lightDirection);

	//绘制原点
	GLUquadricObj *quadricObj = gluNewQuadric();
	gluQuadricDrawStyle(quadricObj, GLU_SMOOTH);
	glColor4f(0.5, 0.5, 0.5, 0.5);
	gluSphere(quadricObj, 0.001, 5, 5);

	if (true == isBackground) {
		//drawBackground();//绘制背景
	}

	if (isShowDisField == true) {
		drawDistFieldGrid();//距离场
	}

	//if (isShowParticle){
	drawParticles();//绘制粒子
					//}

	if (sphManager.fluidSurface== true) {
		drawFluid();//液面
	}

	if (isShowHashGrid) {
		//drawHashGrid();
	}

	drawBoundary();//绘制边界
	
	//drawDistFieldGrid();
	if (isShowBoundaryParticles) {
		drawBoundaryParticles();
	}

	glPopMatrix();

}

//绘制玻璃板
void GLWidget::drawBoundary(GLenum mode)
{
	double woffset = sphManager.wcsphComputer.boundary.boundaryWidth/2.0;
	double hoffset = sphManager.wcsphComputer.boundary.boundaryHeight/2.0;
	double loffset = sphManager.wcsphComputer.boundary.boundaryLength/2.0;

	glPushMatrix();
	glScaled(ratio, ratio, ratio);
	//边界位置和大小
	//double xb = pcisph.xBottle - ballRadius;
	//double yb = pcisph.yBottle - ballRadius;
	//double zb = pcisph.zBottle - ballRadius;
	//double wb = pcisph.wBottle + ballRadius;
	//double hb = pcisph.hBottle + ballRadius;
	//double lb = pcisph.lBottle + ballRadius;

	double xb, yb, zb, wb, hb, lb;

	if (sphManager.sphtype==SPH_TYPE::IISPH)
	{
		 xb = sphManager.iisphComputer.boundary.boundaryOriginX - woffset;
		 yb = sphManager.iisphComputer.boundary.boundaryOriginY - hoffset;
		 zb = sphManager.iisphComputer.boundary.boundaryOriginZ - loffset;
		 wb = sphManager.iisphComputer.boundary.boundaryWidth;
		 hb = sphManager.iisphComputer.boundary.boundaryHeight;
		 lb = sphManager.iisphComputer.boundary.boundaryLength;
	}
	else
	{
		 xb = sphManager.wcsphComputer.boundary.boundaryOriginX - woffset;
		 yb = sphManager.wcsphComputer.boundary.boundaryOriginY - hoffset;
		 zb = sphManager.wcsphComputer.boundary.boundaryOriginZ - loffset;
		 wb = sphManager.wcsphComputer.boundary.boundaryWidth;
		 hb = sphManager.wcsphComputer.boundary.boundaryHeight;
		 lb = sphManager.wcsphComputer.boundary.boundaryLength;
	}
	//绘制边界
	glLineWidth(2);
	glColor4f(0.7, 0.7, 0.9, 0.5);//边界颜色
								  //边线
	glBegin(GL_LINES);
	//后面
	glVertex3d(xb, yb + hb, zb);	//左上
	glVertex3d(xb, yb, zb);
	glVertex3d(xb, yb, zb);			//左下
	glVertex3d(xb + wb, yb, zb);
	glVertex3d(xb + wb, yb, zb);	//右下
	glVertex3d(xb + wb, yb + hb, zb);
	glVertex3d(xb + wb, yb + hb, zb);//右上
	glVertex3d(xb, yb + hb, zb);
	//前面
	glVertex3d(xb, yb + hb, zb + lb);
	glVertex3d(xb, yb, zb + lb);
	glVertex3d(xb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb + lb);
	glVertex3d(xb + wb, yb + hb, zb + lb);
	glVertex3d(xb + wb, yb + hb, zb + lb);
	glVertex3d(xb, yb + hb, zb + lb);
	//左侧
	glVertex3d(xb, yb + hb, zb);	//左上
	glVertex3d(xb, yb + hb, zb + lb);
	glVertex3d(xb, yb, zb);			//左下
	glVertex3d(xb, yb, zb + lb);
	//右侧
	glVertex3d(xb + wb, yb + hb, zb);//右上
	glVertex3d(xb + wb, yb + hb, zb + lb);
	glVertex3d(xb + wb, yb, zb);	//右下
	glVertex3d(xb + wb, yb, zb + lb);
	glEnd();

	//玻璃面
	glColor4f(0.9, 0.9, 1, 0.10);//边界颜色
	glBegin(GL_QUADS);
	//后玻璃
	glVertex3d(xb, yb + hb, zb);
	glVertex3d(xb, yb, zb);
	glVertex3d(xb + wb, yb, zb);
	glVertex3d(xb + wb, yb + hb, zb);
	//下玻璃
	glVertex3d(xb, yb, zb);			//后左下
	glVertex3d(xb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb);	//后右下
									//左玻璃
	glVertex3d(xb, yb + hb, zb);	//后左上
	glVertex3d(xb, yb + hb, zb + lb);//前左上
	glVertex3d(xb, yb, zb + lb);
	glVertex3d(xb, yb, zb);			//后左下
									//右玻璃
	glVertex3d(xb + wb, yb + hb, zb);//后右上
	glVertex3d(xb + wb, yb + hb, zb + lb);
	glVertex3d(xb + wb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb);	//后右下
	glEnd();

	glPopMatrix();
	
}
void GLWidget::drawBoundaryParticles()
{

	double woffset = sphManager.wcsphComputer.boundary.boundaryWidth / 2.0;
	double hoffset = sphManager.wcsphComputer.boundary.boundaryHeight / 2.0;
	double loffset = sphManager.wcsphComputer.boundary.boundaryLength / 2.0f;

	GLUquadricObj *quadricObj = gluNewQuadric();
	gluQuadricDrawStyle(quadricObj, GLU_SMOOTH);
	//粒子颜色
	glColor4f(0.5, 0.55, 0.7, 0.3);
	glPushMatrix();
	if (sphManager.sphtype == SPH_TYPE::IISPH)
	{
		int particleNumber = sphManager.iisphComputer.fluidModel.boundaryObj.staticRigidParticleList.size();
		int boundaryNumber = sphManager.iisphComputer.fluidModel.boundaryNumber;
		for (unsigned int i = 0; i < boundaryNumber; i++)
		{
			float x = (float)(sphManager.iisphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.x - woffset);
			float y = (float)(sphManager.iisphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.z - hoffset);
			float z = (float)(sphManager.iisphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.y - loffset);
			glPushMatrix();
			glTranslatef(x * ratio, y * ratio, z * ratio);
			gluSphere(quadricObj, sphManager.iisphComputer.boundary.pRad*ratio, 16, 16);
			glPopMatrix();
		}
	}
	else
	{
		int particleNumber = sphManager.wcsphComputer.fluidModel.boundaryObj.staticRigidParticleList.size();
		for (unsigned int i = 0; i < particleNumber; i++)
		{
			float x = (float)(sphManager.wcsphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.x - woffset);
			float y = (float)(sphManager.wcsphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.z - hoffset);
			float z = (float)(sphManager.wcsphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.y - loffset);
			glPushMatrix();
			glTranslatef(x * ratio, y * ratio, z * ratio);
			gluSphere(quadricObj, sphManager.wcsphComputer.boundary.pRad*ratio, 16, 16);

			glPopMatrix();
		}
	}
	glPopMatrix();
}
void GLWidget::drawParticles()
{
	double woffset = sphManager.wcsphComputer.boundary.boundaryWidth / 2.0;
	double hoffset = sphManager.wcsphComputer.boundary.boundaryHeight / 2.0;
	double loffset = sphManager.wcsphComputer.boundary.boundaryLength / 2.0f;

	GLUquadricObj *quadricObj = gluNewQuadric();
	gluQuadricDrawStyle(quadricObj, GLU_SMOOTH);

	//粒子颜色
	glColor4f(0.5, 0.55, 0.7, 0.9);
	glPushMatrix();
	if (sphManager.sphtype == SPH_TYPE::IISPH)
	{
		//绘制IISPH
		int particleNumber = sphManager.iisphComputer.fluidModel.particleList.size();
		for (unsigned int i = 0; i < particleNumber / 2; i++)
		{
			float x = (float)(sphManager.iisphComputer.fluidModel.particleList[i]->position.x - woffset);
			float y = (float)(sphManager.iisphComputer.fluidModel.particleList[i]->position.z - hoffset);
			float z = (float)(sphManager.iisphComputer.fluidModel.particleList[i]->position.y - loffset);
			glPushMatrix();
			glTranslatef(x * ratio, y * ratio, z * ratio);
			//double rrr = pow(sphManager.iisphComputer.fluidModel.particleList[i]->mopt*1.25 / sphManager.iisphComputer.fluidModel.particleList[i]->density, double(1.0 / 3.0))*0.5;
			//cout << sphManager.iisphComputer.fluidModel.particleList[i]->mopt << endl;
			gluSphere(quadricObj, sphManager.iisphComputer.fluidModel.particleList[i]->particleRad*ratio, 16, 16);
			//gluSphere(quadricObj, rrr*ratio, 16, 16);
			glPopMatrix();
		}
		if (sphManager.iisphComputer.fluidModel.doubleDamBreak==true)
		{
			glColor4f(0.7, 0.7, 0.0, 0.9);
		}
		for (unsigned int i = particleNumber / 2; i < particleNumber ; i++)
		{
			float x = (float)(sphManager.iisphComputer.fluidModel.particleList[i]->position.x - woffset);
			float y = (float)(sphManager.iisphComputer.fluidModel.particleList[i]->position.z - hoffset);
			float z = (float)(sphManager.iisphComputer.fluidModel.particleList[i]->position.y - loffset);
			glPushMatrix();
			glTranslatef(x * ratio, y * ratio, z * ratio);
			//double rrr = pow(sphManager.iisphComputer.fluidModel.particleList[i]->mopt*1.25 / sphManager.iisphComputer.fluidModel.particleList[i]->density, double(1.0 / 3.0))*0.5;
			//cout << sphManager.iisphComputer.fluidModel.particleList[i]->mopt << endl;
			gluSphere(quadricObj, sphManager.iisphComputer.fluidModel.particleList[i]->particleRad*ratio, 16, 16);
			//gluSphere(quadricObj, rrr*ratio, 16, 16);
			glPopMatrix();
		}
		if (sphManager.iisphComputer.fluidModel.staticPilarFlag == true)
		{
			glColor4f(0.9, 0.55, 0.7, 0.9);
			int particleNumber = sphManager.iisphComputer.fluidModel.boundaryObj.staticRigidParticleList.size();
			int boundaryNumber = sphManager.iisphComputer.fluidModel.boundaryNumber;
			for (unsigned int i = boundaryNumber; i < particleNumber; i++)
			{
				float x = (float)(sphManager.iisphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.x - woffset);
				float y = (float)(sphManager.iisphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.z - hoffset);
				float z = (float)(sphManager.iisphComputer.fluidModel.boundaryObj.staticRigidParticleList[i]->position.y - loffset);
				glPushMatrix();
				glTranslatef(x * ratio, y * ratio, z * ratio);
				gluSphere(quadricObj, sphManager.iisphComputer.boundary.pRad*ratio, 16, 16);
				glPopMatrix();
			}
		}
	}
	else
	{
		//绘制WCSPH
		int particleNumber = sphManager.wcsphComputer.fluidModel.particleList.size();
		for (unsigned int i = 0; i < particleNumber; i++)
		{
			float x = (float)(sphManager.wcsphComputer.fluidModel.particleList[i]->position.x - woffset);
			float y = (float)(sphManager.wcsphComputer.fluidModel.particleList[i]->position.z - hoffset);
			float z = (float)(sphManager.wcsphComputer.fluidModel.particleList[i]->position.y - loffset);
			glPushMatrix();
			glTranslatef(x * ratio, y * ratio, z * ratio);
			//WCSPH无法通过粒子半径变化进行交互，因此半径固定
			gluSphere(quadricObj, sphManager.wcsphComputer.fluidModel.particleRad*ratio, 16, 16);
			glPopMatrix();
		}
	}
	glPopMatrix();
}

void GLWidget::drawDynamicRigidParticles()
{

	double woffset = sphManager.wcsphComputer.boundary.boundaryWidth / 2.0;
	double hoffset = sphManager.wcsphComputer.boundary.boundaryHeight / 2.0;
	double loffset = sphManager.wcsphComputer.boundary.boundaryLength / 2.0f;

	GLUquadricObj *quadricObj = gluNewQuadric();
	gluQuadricDrawStyle(quadricObj, GLU_SMOOTH);
	//粒子颜色
	glColor4f(0.9, 0.55, 0.1, 0.9);
	glPushMatrix();
	if (sphManager.sphtype == SPH_TYPE::IISPH)
	{
		int particleNumber = sphManager.iisphComputer.fluidModel.dynamicCube.dynamicParticles.size();
		for (unsigned int i = 0; i < particleNumber; i++)
		{
			float x = (float)(sphManager.iisphComputer.fluidModel.dynamicCube.dynamicParticles[i]->m_x.x - woffset);
			float y = (float)(sphManager.iisphComputer.fluidModel.dynamicCube.dynamicParticles[i]->m_x.z - hoffset);
			float z = (float)(sphManager.iisphComputer.fluidModel.dynamicCube.dynamicParticles[i]->m_x.y - loffset);
			glPushMatrix();
			glTranslatef(x * ratio, y * ratio, z * ratio);
			gluSphere(quadricObj, sphManager.iisphComputer.boundary.pRad*ratio, 16, 16);
			glPopMatrix();
		}
	}
	glPopMatrix();
}

void GLWidget::drawDistFieldGrid()
{
	double woffset = 0.5f;
	double loffset = 1.0f;
	double hoffset = 0.5f;

	//颜色
	glColor4f(0.8, 0.1, 0.1, 1.0);
	glPointSize(2);

	glPushMatrix();
	glScaled(ratio, ratio, ratio);
	glBegin(GL_POINTS);

	DistanceNode* gridNode = NULL;
	for (int i = 0; i < sphManager.iisphComputer.distanceField->hash.size(); i++)
	{
		gridNode = sphManager.iisphComputer.distanceField->hash[i];

		while (gridNode != NULL)
		{

			if (gridNode->data < 0)
			{
				//液面内
				glColor4f(0.9, 0.9, 0.1, 1.0);
			}
			else if (gridNode->data < 0.07f)
			{
				//液面内
				glColor4f(0.1, 0.9, 0.9, 1.0);
			}
			else
			{
				//液面外
				glColor4f(0.1, 0.1, 0.5, 0.0);
			}

			glPushMatrix();
			glVertex3d((gridNode->x - woffset) , (gridNode->z - hoffset), (gridNode->y - loffset));
			glPopMatrix();

			gridNode = gridNode->next;
		}

	}

	glEnd();
	glPopMatrix();
}

//绘制页面
void GLWidget::drawFluid()
{
	glPushMatrix();

	//液面颜色
	glColor4f(0.4, 0.5, 0.7, 1);
	//glColor4f(0, 0, 0.8, 1);

	glScaled(ratio, ratio, ratio);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	//glEnableClientState(GL_COLOR_ARRAY);

	//定义顶点数组  
	glVertexPointer(3, GL_DOUBLE, 0, sphManager.iisphComputer.distanceField->mesh.data());
	glNormalPointer(GL_DOUBLE, 0, sphManager.iisphComputer.distanceField->meshVector.data());
	//glColorPointer(4, GL_FLOAT, 0, colorBackgroud.data());

	//绘制
	glDrawArrays(GL_TRIANGLES, 0, sphManager.iisphComputer.distanceField->mesh.size() / 3);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	//glDisableClientState(GL_COLOR_ARRAY);

	glPopMatrix();
}

void GLWidget::StartRun()
{
	//启动动画帧
	timer.start(timeStep);
	sphManager.Run();
}

void GLWidget::StopRun()
{
	//停止动画帧
	timer.stop();
	sphManager.Stop();
}

void GLWidget::setOutputFramePicture(bool isOutput)
{
	isOutputFramePicture = isOutput;
}

void GLWidget::OutputFramePicture()
{
	char filePath[32] = "PCISPHImage\\frameImage";//当前帧picture的路径

	char fName[32] = "";
	char file_out_path[32] = "d:\\";
	strcat(fName, file_out_path);
	strcat(fName, "PCISPHImage");
	mkdir(fName);

	QString fullPath = QString("%1%2%3.png").arg(file_out_path).arg(filePath).arg(sumFrame);

	QImage image = grabFramebuffer();

	image.save(fullPath, "PNG");
}
void GLWidget::setShowBoundaryParticles(bool isShow)
{
	isShowBoundaryParticles = isShow;
}
