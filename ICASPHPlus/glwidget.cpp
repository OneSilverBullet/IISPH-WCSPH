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

//**********************************��ʼ��*********************************************
void GLWidget::initWidget()
{
	//���õ�ǰ��ʼ��
	cout << "1" << endl;
	sphManager.Initialize();

	isMouseDown = false;			//����Ƿ���
	frame = 0;
	sumFrame = 0;
	frameSec = 0;
	rotX = -90;//��ʼ����ת�Ƕ�
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

	//���ò���
	//ratio = 3.6 / pcisph.wBottle;//���ű���
	ratio = 0.5;
	//ballRadius = 0.01;	//���Ӱ뾶
	timeStep = 10;		//֡���,ms

	lightPosition[0] = 0.0f;	//��Դλ��,xyz,dis
	lightPosition[1] = 1.5f;
	lightPosition[2] = 0.0f;
	lightPosition[3] = 0.0f;
	lightDirection[0] = 0.0f;	//���շ���
	lightDirection[0] = -1.0f;
	lightDirection[0] = 0.0f;

	//���ò���ȡ����
	//double bgy = pcisph.yBottle * ratio - 0.1;//����y����
	//setBackground(30, 30, bgy, 1.0);
	setFocusPolicy(Qt::StrongFocus);

	/*****��ʼģ��*****/
	
	//timer.start(timeStep);
	timer.stop();

	//��ʼʱ�ӣ��Լ�����һ֡
	time.start();
	//timeFrame();

	QObject::connect(&timer, SIGNAL(timeout()), this, SLOT(timeFrame()));
	//qDebug() << QString::fromLocal8Bit("��ǰ�߳�id��") << QThread::currentThread();

}

void GLWidget::initializeGL()
{
	glClearColor(0, 0, 0, 1);//����ɫ

	glEnable(GL_DEPTH_TEST);//��Ȳ��ԣ����ֵ��Ч
	glEnable(GL_BLEND);		//��ɫ���
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);//������ɫ׷��
	glShadeModel(GL_SMOOTH);
	glClearDepth(1.0f);		//������ȷ�ΧΪ���

							//���ù��ղ���
	setLight();
}

void GLWidget::setLight()
{
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	//GLfloat lightPosition[] = { 0.0, 1.2, -1.0, 1.0 };		//��Դλ��,xyz,dis
	GLfloat lightAmbient[] = { 0.6, 0.6, 0.6, 1.0 };		//ȫ�ֹ�����,��ǿ��
	GLfloat lightDiffuse[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	//GLfloat lightDirection[] = { 0.0f, -1.0f, 0.0f };	//���շ���
	GLfloat lightExponent[] = { 1.0f };					//�۹�̶�
	GLfloat lightCutoff[] = { 60.0f };					//��Դ�н�

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);//ָ����Ҫ��ɫ׷�ٵĲ�������
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);//ָ����0�Ź�Դ��λ�� 
													 //����
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse); //�������
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);//���淴���
													 //�۹�
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

	//ͶӰ��ʽ��͸��ͶӰ
	gluPerspective(60, (GLfloat)iWinWidth / (GLfloat)iWinHeight, 0.1, 100);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
}

//**********************************�¼�*********************************************
void GLWidget::mouseMoveEvent(QMouseEvent *e)
{
	if (isMouseDown)
	{
		int deltX, deltY;
		// ��������ƶ����룬��תͼ��
		deltX = e->x() - mouseX;
		deltY = e->y() - mouseY;

		//if (false == isLock)
		{
			// ��ת��  
			rotX += deltX / 2;
			rotY += deltY / 2;
			// ��ת�ǲ�����360�� 
			rotX = fmodf(rotX, 360.0);
			rotY = fmodf(rotY, 360.0);
		}
		//else
		{
			//pcisph.xBottle += double(deltX) / this->width() / ratio * 4;
			//pcisph.yBottle -= double(deltY) / this->height() / ratio * 4;
		}

		//���µ�ǰ���λ�ã�ʹͼ��ʵʱ��ת
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
	int numDegrees = e->delta() / 8;//�����ĽǶȣ�*8�����������ľ���
	double numSteps = numDegrees / 120.0;//�����Ĳ�����*15�����������ĽǶ�
	zoom += numSteps / 2;
	e->accept();      //���ո��¼�
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
	//ˢ�µ�ǰ��ͼ��
	this->update();

	//���Ƶ�ǰ֡��
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

//**********************************����*********************************************
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

	//�ƶ���Դͬ���ƶ�
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);//ָ����0�Ź�Դ��λ�� 
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, lightDirection);

	//����ԭ��
	GLUquadricObj *quadricObj = gluNewQuadric();
	gluQuadricDrawStyle(quadricObj, GLU_SMOOTH);
	glColor4f(0.5, 0.5, 0.5, 0.5);
	gluSphere(quadricObj, 0.001, 5, 5);

	if (true == isBackground) {
		//drawBackground();//���Ʊ���
	}

	if (isShowDisField == true) {
		drawDistFieldGrid();//���볡
	}

	//if (isShowParticle){
	drawParticles();//��������
					//}

	if (sphManager.fluidSurface== true) {
		drawFluid();//Һ��
	}

	if (isShowHashGrid) {
		//drawHashGrid();
	}

	drawBoundary();//���Ʊ߽�
	
	//drawDistFieldGrid();
	if (isShowBoundaryParticles) {
		drawBoundaryParticles();
	}

	glPopMatrix();

}

//���Ʋ�����
void GLWidget::drawBoundary(GLenum mode)
{
	double woffset = sphManager.wcsphComputer.boundary.boundaryWidth/2.0;
	double hoffset = sphManager.wcsphComputer.boundary.boundaryHeight/2.0;
	double loffset = sphManager.wcsphComputer.boundary.boundaryLength/2.0;

	glPushMatrix();
	glScaled(ratio, ratio, ratio);
	//�߽�λ�úʹ�С
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
	//���Ʊ߽�
	glLineWidth(2);
	glColor4f(0.7, 0.7, 0.9, 0.5);//�߽���ɫ
								  //����
	glBegin(GL_LINES);
	//����
	glVertex3d(xb, yb + hb, zb);	//����
	glVertex3d(xb, yb, zb);
	glVertex3d(xb, yb, zb);			//����
	glVertex3d(xb + wb, yb, zb);
	glVertex3d(xb + wb, yb, zb);	//����
	glVertex3d(xb + wb, yb + hb, zb);
	glVertex3d(xb + wb, yb + hb, zb);//����
	glVertex3d(xb, yb + hb, zb);
	//ǰ��
	glVertex3d(xb, yb + hb, zb + lb);
	glVertex3d(xb, yb, zb + lb);
	glVertex3d(xb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb + lb);
	glVertex3d(xb + wb, yb + hb, zb + lb);
	glVertex3d(xb + wb, yb + hb, zb + lb);
	glVertex3d(xb, yb + hb, zb + lb);
	//���
	glVertex3d(xb, yb + hb, zb);	//����
	glVertex3d(xb, yb + hb, zb + lb);
	glVertex3d(xb, yb, zb);			//����
	glVertex3d(xb, yb, zb + lb);
	//�Ҳ�
	glVertex3d(xb + wb, yb + hb, zb);//����
	glVertex3d(xb + wb, yb + hb, zb + lb);
	glVertex3d(xb + wb, yb, zb);	//����
	glVertex3d(xb + wb, yb, zb + lb);
	glEnd();

	//������
	glColor4f(0.9, 0.9, 1, 0.10);//�߽���ɫ
	glBegin(GL_QUADS);
	//����
	glVertex3d(xb, yb + hb, zb);
	glVertex3d(xb, yb, zb);
	glVertex3d(xb + wb, yb, zb);
	glVertex3d(xb + wb, yb + hb, zb);
	//�²���
	glVertex3d(xb, yb, zb);			//������
	glVertex3d(xb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb);	//������
									//����
	glVertex3d(xb, yb + hb, zb);	//������
	glVertex3d(xb, yb + hb, zb + lb);//ǰ����
	glVertex3d(xb, yb, zb + lb);
	glVertex3d(xb, yb, zb);			//������
									//�Ҳ���
	glVertex3d(xb + wb, yb + hb, zb);//������
	glVertex3d(xb + wb, yb + hb, zb + lb);
	glVertex3d(xb + wb, yb, zb + lb);
	glVertex3d(xb + wb, yb, zb);	//������
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
	//������ɫ
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

	//������ɫ
	glColor4f(0.5, 0.55, 0.7, 0.9);
	glPushMatrix();
	if (sphManager.sphtype == SPH_TYPE::IISPH)
	{
		//����IISPH
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
		//����WCSPH
		int particleNumber = sphManager.wcsphComputer.fluidModel.particleList.size();
		for (unsigned int i = 0; i < particleNumber; i++)
		{
			float x = (float)(sphManager.wcsphComputer.fluidModel.particleList[i]->position.x - woffset);
			float y = (float)(sphManager.wcsphComputer.fluidModel.particleList[i]->position.z - hoffset);
			float z = (float)(sphManager.wcsphComputer.fluidModel.particleList[i]->position.y - loffset);
			glPushMatrix();
			glTranslatef(x * ratio, y * ratio, z * ratio);
			//WCSPH�޷�ͨ�����Ӱ뾶�仯���н�������˰뾶�̶�
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
	//������ɫ
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

	//��ɫ
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
				//Һ����
				glColor4f(0.9, 0.9, 0.1, 1.0);
			}
			else if (gridNode->data < 0.07f)
			{
				//Һ����
				glColor4f(0.1, 0.9, 0.9, 1.0);
			}
			else
			{
				//Һ����
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

//����ҳ��
void GLWidget::drawFluid()
{
	glPushMatrix();

	//Һ����ɫ
	glColor4f(0.4, 0.5, 0.7, 1);
	//glColor4f(0, 0, 0.8, 1);

	glScaled(ratio, ratio, ratio);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	//glEnableClientState(GL_COLOR_ARRAY);

	//���嶥������  
	glVertexPointer(3, GL_DOUBLE, 0, sphManager.iisphComputer.distanceField->mesh.data());
	glNormalPointer(GL_DOUBLE, 0, sphManager.iisphComputer.distanceField->meshVector.data());
	//glColorPointer(4, GL_FLOAT, 0, colorBackgroud.data());

	//����
	glDrawArrays(GL_TRIANGLES, 0, sphManager.iisphComputer.distanceField->mesh.size() / 3);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	//glDisableClientState(GL_COLOR_ARRAY);

	glPopMatrix();
}

void GLWidget::StartRun()
{
	//��������֡
	timer.start(timeStep);
	sphManager.Run();
}

void GLWidget::StopRun()
{
	//ֹͣ����֡
	timer.stop();
	sphManager.Stop();
}

void GLWidget::setOutputFramePicture(bool isOutput)
{
	isOutputFramePicture = isOutput;
}

void GLWidget::OutputFramePicture()
{
	char filePath[32] = "PCISPHImage\\frameImage";//��ǰ֡picture��·��

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
