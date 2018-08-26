#include "Boundary.h"

Boundary::Boundary()
{
	GenerateNewBoundary(boundaryOriginX, boundaryOriginY, boundaryOriginZ, boundaryWidth, boundaryLength, boundaryHeight, offset);
}

Boundary::Boundary(Boundary & a)
{
	operator=(a);
}

Boundary::Boundary(double x, double y, double z, double w, double l, double h, double os)
{
	GenerateNewBoundary(x, y, z, w, l, h, os);

}

Boundary::~Boundary()
{
}

void Boundary::GenerateNewBoundary(double x, double y, double z, double w, double l, double h, double os)
{
	boundaryOriginX = x;
	boundaryOriginY = y;
	boundaryOriginZ = z;
	boundaryWidth = w;
	boundaryLength = l;
	boundaryHeight = h;
	offset = os;
	gridBoundaryOriginX = boundaryOriginX - offset;
	gridBoundaryOriginY = boundaryOriginY - offset;
	gridBoundaryOriginZ = boundaryOriginZ - offset;
	gridBoundaryWidth = boundaryWidth + 2 * offset;
	gridBoundaryLength = boundaryLength + 2 * offset;
	gridBoundaryHeight = boundaryHeight + 2 * offset;


	distanceFieldOffset = double(2.0 * offset / 3.0);
	distanceFieldOriginX = boundaryOriginX - distanceFieldOffset;
	distanceFieldOriginY = boundaryOriginX - distanceFieldOffset;
	distanceFieldOriginZ = boundaryOriginX - distanceFieldOffset;
	distanceFieldWidth = boundaryWidth + 2 * distanceFieldOffset;
	distanceFieldLength = boundaryLength + 2 * distanceFieldOffset;
	distanceFieldHeight = boundaryHeight + 2 * distanceFieldOffset;


	//设置新值
	Vector3f origin(x, y, z);
	Vector3f wlh(w, l, h);
	//当初始化一个boundary时候，立即对SharedData进行设置
	SharedData::SetBoundaryOrigin(origin);
	SharedData::SetBoundaryWLH(wlh);
	SharedData::SetOffset(os);
}

void Boundary::GenerateNewBoundary(Vector3f position, Vector3f fluidWLH, double offs)
{
	boundaryOriginX = position.GetX();
	boundaryOriginY = position.GetX();
	boundaryOriginZ = position.GetZ();
	boundaryWidth = fluidWLH.GetX();
	boundaryLength = fluidWLH.GetY();
	boundaryHeight = fluidWLH.GetZ();
	offset = offs;
	gridBoundaryOriginX = boundaryOriginX - offset;
	gridBoundaryOriginY = boundaryOriginY - offset;
	gridBoundaryOriginZ = boundaryOriginZ - offset;
	gridBoundaryWidth = boundaryWidth + 2 * offset;
	gridBoundaryLength = boundaryLength + 2 * offset;
	gridBoundaryHeight = boundaryHeight + 2 * offset;

	distanceFieldOffset = double(2.0 * offset / 3.0);
	distanceFieldOriginX = distanceFieldOriginX - distanceFieldOffset;
	distanceFieldOriginY = distanceFieldOriginY - distanceFieldOffset;
	distanceFieldOriginZ = distanceFieldOriginZ - distanceFieldOffset;
	distanceFieldWidth = boundaryWidth + 2 * distanceFieldOffset;
	distanceFieldLength = boundaryLength + 2 * distanceFieldOffset;
	distanceFieldHeight = boundaryHeight + 2 * distanceFieldOffset;
}

Boundary & Boundary::operator=(Boundary & a)
{
	boundaryOriginX = a.boundaryOriginX;
	boundaryOriginY = a.boundaryOriginY;
	boundaryOriginZ = a.boundaryOriginZ;
	boundaryWidth = a.boundaryWidth;
	boundaryLength = a.boundaryLength;
	boundaryHeight = a.boundaryHeight;
	offset = a.offset;
	gridBoundaryOriginX = a.gridBoundaryOriginX;
	gridBoundaryOriginY = a.gridBoundaryOriginY;
	gridBoundaryOriginZ = a.gridBoundaryOriginZ;
	gridBoundaryWidth = a.gridBoundaryWidth;
	gridBoundaryLength = a.gridBoundaryLength;
	gridBoundaryHeight = a.gridBoundaryHeight;

	distanceFieldOffset = a.distanceFieldOffset;
	distanceFieldOriginX = a.distanceFieldOriginX;
	distanceFieldOriginY = a.distanceFieldOriginY;
	distanceFieldOriginZ = a.distanceFieldOriginZ;
	distanceFieldWidth = a.distanceFieldWidth;
	distanceFieldLength = a.distanceFieldLength;
	distanceFieldHeight = a.distanceFieldHeight;

	pRad = a.pRad;
	return *this;
}

//依据边界的具体信息，建立粒子墙，计算好每个粒子的位置
vector<Vector3f> Boundary::initBoundaryData()
{
	std::cout << "初始化边界数据" << std::endl;
	std::vector<Vector3f> boundaryParticles;
	//StaticRigidObject *rb = new StaticRigidObject();
	const double particleRad = pRad;   //硬性规定为最小半径
	const double diam = particleRad * 2; //得到当前粒子的直径

	Vector3f boundaryMinP(boundaryOriginX-particleRad, boundaryOriginY- particleRad, boundaryOriginZ- particleRad);
	Vector3f boundaryMaxP(boundaryOriginX + boundaryWidth + particleRad, boundaryOriginY + boundaryLength + particleRad, boundaryOriginZ + boundaryHeight + particleRad);

	//使用直径来作为每一步向前的变化量
	double xshift = diam;
	double yshift = diam;
	double zshift = diam;

	const int stepsX = (int)round((boundaryWidth + 2 * diam) / xshift);
	const int stepsY = (int)round((boundaryLength + 2 * diam) / yshift);
	const int stepsZ = (int)round((boundaryHeight + 2 * diam) / zshift);

	Vector3f start = boundaryMinP;

	//六面墙，使用粒子进行墙的创建
	for (int i = 0; i < stepsX; i++)
	{
		for (int j = 0; j < stepsY; j++)
		{
			Vector3f currPos = Vector3f(i*xshift, j*yshift, 0) + boundaryMinP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = Vector3f(i*xshift, 0, k*zshift) + boundaryMinP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int j = 0; j < stepsY; j++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = Vector3f(0, j*yshift, k*zshift) + boundaryMinP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int j = 0; j < stepsY; j++)
		{
			Vector3f currPos = -Vector3f(i*xshift, j*yshift, 0) + boundaryMaxP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int i = 0; i < stepsX; i++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = -Vector3f(i*xshift, 0, k*zshift) + boundaryMaxP;
			boundaryParticles.push_back(currPos);
		}
	}
	for (int j = 0; j < stepsY; j++)
	{
		for (int k = 0; k < stepsZ; k++)
		{
			Vector3f currPos = -Vector3f(0, j*yshift, k*zshift) + boundaryMaxP;
			boundaryParticles.push_back(currPos);
		}
	}
	
	boundaryNumber = boundaryParticles.size();
	cout << boundaryNumber << endl;
	//绘制了一个柱子
	if (staticRigidFlag==true)
	{
		//得到当前装入容器当中边界粒子的个数
		int boundaryParticleNum = boundaryParticles.size();
		//下面开始绘制静态柱子
		//起始点
		double staticPilar_OriginX = boundaryOriginX + 0.4*boundaryWidth;
		double staticPilar_OriginY = boundaryOriginY + 0.4*boundaryLength;
		double staticPilar_OriginZ = boundaryOriginZ + 0.0*boundaryHeight;
		//柱子的长宽高
		double staticPilar_Width = 0.2f*boundaryWidth;
		double staticPilar_Length = 0.2f*boundaryLength;
		double staticPilar_Height = 0.8f*boundaryHeight;

		const int stepX1 = staticPilar_Width / xshift;
		const int stepY1 = staticPilar_Length / yshift;
		const int stepsZ1 = staticPilar_Height / zshift;
		cout << stepX1 << " " << stepY1 << " " << stepsZ1 << endl;

		for (int m = 0; m < stepX1; m++)
		{
			for (int n = 0; n < stepY1; n++)
			{
				for (int p = 0; p < stepsZ1; p++)
				{
					Vector3f currPos = Vector3f(m*xshift + staticPilar_OriginX, n*yshift + staticPilar_OriginY, p*zshift + staticPilar_OriginZ);
					boundaryParticles.push_back(currPos);
				}
			}
		}
	}
	cout << boundaryParticles.size() << endl;
	//std::cout << "boundaryParticles.size():" << boundaryParticles.size() << std::endl;
	//model.addRigidBodyObject(rb, static_cast<unsigned int>(boundaryParticles.size()), &boundaryParticles[0]);
	std::cout << "初始化边界数据完成" << std::endl;
	return boundaryParticles;
}
