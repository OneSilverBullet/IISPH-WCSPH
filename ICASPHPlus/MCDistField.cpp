#include "MCDistField.h"


MCDistField::MCDistField()
{
}

MCDistField::~MCDistField()
{
}

MCDistField::MCDistField(double len, uint hashSize, double xSize, double ySize, double zSize, double xOrigin,
	double yOrigin, double zOrigin) : DistanceField(len, hashSize, xSize, ySize, zSize, xOrigin, yOrigin, zOrigin)
{
	mesh.clear();
}

void MCDistField::GenerateMesh()
{

	double woffset = -0.7f;
	double loffset = 1.0f;
	double hoffset = 0.5f;

	mesh.clear();
	meshVector.clear();

	Vec3D vertexCell[8];//��Ԫ�嶥��
	Vec3D edgeVertexInter[12];//��Ԫ��߲�ֵ��ĵ�
	double vertexDisValue[8];//�������ֵ
	int vertexFlag = 0;//�����־λ

	//��������Ԫ������
	int xSum = xLen / lenght + 1;
	int ySum = yLen / lenght + 1;
	int zSum = zLen / lenght + 1;

	//��ǰ�����嵥Ԫ����㣬(OpenGL�����º�)
	double curx = ox;
	double cury = oy;
	double curz = oz;

	DistanceNode* gridNode[8];
	int triCell[16];//��ǰ��Ԫ���mesh
	//����ÿ���ڵ�����꣬����ÿ���ڵ�
    #pragma omp parallel for private(gridNode,triCell,vertexCell,edgeVertexInter,vertexDisValue,vertexFlag,curx,cury,curz)
	for (int i = 0; i < xSum*ySum*zSum; i++)
	{
		curx = idxGridCatalog[i].x;
		cury = idxGridCatalog[i].y;
		curz = idxGridCatalog[i].z;
		vertexFlag = 0;

		//��Ԫ���0���㣬��ǰ��
		/*
				1 - - - -5
				/|	    /|
				/ |      / |
				0 -|- - -4  |
				|	2 - - |- 6
				| /	  | /
				|/       |/
				3 - - - -7
				*/
		double xv0 = curx;
		double yv0 = cury + lenght;
		double zv0 = curz + lenght;
		int cellVertexIdx = 0;

		//������Ԫ��ÿ������
		vertexCell[0] = Vec3D(xv0, yv0, zv0);
		vertexCell[1] = Vec3D(xv0, yv0, zv0 - lenght);
		vertexCell[2] = Vec3D(xv0, yv0 - lenght, zv0 - lenght);
		vertexCell[3] = Vec3D(xv0, yv0 - lenght, zv0);
		vertexCell[4] = Vec3D(xv0 + lenght, yv0, zv0);
		vertexCell[5] = Vec3D(xv0 + lenght, yv0, zv0 - lenght);
		vertexCell[6] = Vec3D(xv0 + lenght, yv0 - lenght, zv0 - lenght);
		vertexCell[7] = Vec3D(xv0 + lenght, yv0 - lenght, zv0);
		//��ȡ��Ԫ�嶥��
		gridNode[0] = GetGridNode(xv0, yv0, zv0);
		gridNode[1] = GetGridNode(xv0, yv0, zv0 - lenght);
		gridNode[2] = GetGridNode(xv0, yv0 - lenght, zv0 - lenght);
		gridNode[3] = GetGridNode(xv0, yv0 - lenght, zv0);
		gridNode[4] = GetGridNode(xv0 + lenght, yv0, zv0);
		gridNode[5] = GetGridNode(xv0 + lenght, yv0, zv0 - lenght);
		gridNode[6] = GetGridNode(xv0 + lenght, yv0 - lenght, zv0 - lenght);
		gridNode[7] = GetGridNode(xv0 + lenght, yv0 - lenght, zv0);
		int iv = 0;//��������
		int positive = 0;
		int negative = 0;
		for (; iv < 8; iv++)
		{
			if (gridNode[iv] == NULL){
				break;
			}
			else if (gridNode[iv]->data>0)
			{
				positive++;
				vertexDisValue[iv] = gridNode[iv]->data;
			}
			else if (gridNode[iv]->data < 0)
			{
				negative++;
				vertexDisValue[iv] = gridNode[iv]->data;
			}
		}
		if (iv >= 8 && positive != 8 && negative != 8)
		{
			double isov = 0.0;//��ֵ��
			edgeVertexInter[0] = VertexInter(vertexCell[0], vertexCell[1], vertexDisValue[0], vertexDisValue[1], isov);
			edgeVertexInter[1] = VertexInter(vertexCell[1], vertexCell[2], vertexDisValue[1], vertexDisValue[2], isov);
			edgeVertexInter[2] = VertexInter(vertexCell[2], vertexCell[3], vertexDisValue[2], vertexDisValue[3], isov);
			edgeVertexInter[3] = VertexInter(vertexCell[3], vertexCell[0], vertexDisValue[3], vertexDisValue[0], isov);
			edgeVertexInter[4] = VertexInter(vertexCell[4], vertexCell[5], vertexDisValue[4], vertexDisValue[5], isov);
			edgeVertexInter[5] = VertexInter(vertexCell[5], vertexCell[6], vertexDisValue[5], vertexDisValue[6], isov);
			edgeVertexInter[6] = VertexInter(vertexCell[6], vertexCell[7], vertexDisValue[6], vertexDisValue[7], isov);
			edgeVertexInter[7] = VertexInter(vertexCell[7], vertexCell[4], vertexDisValue[7], vertexDisValue[4], isov);
			edgeVertexInter[8] = VertexInter(vertexCell[0], vertexCell[4], vertexDisValue[0], vertexDisValue[4], isov);
			edgeVertexInter[9] = VertexInter(vertexCell[1], vertexCell[5], vertexDisValue[1], vertexDisValue[5], isov);
			edgeVertexInter[10] = VertexInter(vertexCell[2], vertexCell[6], vertexDisValue[2], vertexDisValue[6], isov);
			edgeVertexInter[11] = VertexInter(vertexCell[3], vertexCell[7], vertexDisValue[3], vertexDisValue[7], isov);


			//���¶����Ӧ�ı�־λ
			for (int ivd = 0; ivd < 8; ivd++)
			{
				//�����ڵĵ���Ӧλ�ã���ʽ����С�ڵ�ֵ��ĵ�Ϊ�����ڲ��ĵ�
				if (vertexDisValue[ivd] <= isov)
				{
					vertexFlag = vertexFlag | (1 << ivd);
				}
			}

			if (vertexFlag >= 256)
				continue;

			//��ȡ��Ԫ��mesh�Ķ�������
			for (int itt = 0; itt < 16; itt++)
			{
				triCell[itt] = triTable[vertexFlag][itt];
			}

			//���ݶ�����������ȡ�����浱ǰ��Ԫ��������ÿ������
            #pragma omp critical(max_arx)
			{
			for (int itc = 0; itc < 16 && triCell[itc] != -1; itc++)
			{
				mesh.push_back(edgeVertexInter[triCell[itc]].x - woffset);
				mesh.push_back(edgeVertexInter[triCell[itc]].z - hoffset);
				mesh.push_back(edgeVertexInter[triCell[itc]].y - loffset);
			}

			////��������
			Vec3D normal;
			for (int itcn = 0; itcn < 16 && triCell[itcn] != -1; itcn += 3)
			{
				computeNormal(edgeVertexInter[triCell[itcn + 0]], edgeVertexInter[triCell[itcn + 1]],
					edgeVertexInter[triCell[itcn + 2]], normal);
				for (int itcnv = 0; itcnv < 3; itcnv++)
				{
					meshVector.push_back(normal.x - woffset);
					meshVector.push_back(normal.z - hoffset);
					meshVector.push_back(normal.y - loffset);
				}
			}
			}
		}
	}
}

Vec3D MCDistField::VertexInter(Vec3D &p1, Vec3D &p2, double v1, double v2, double isovalue)
{
	if (v1 == v2)
		return Vec3D(MAX_DIS, MAX_DIS, MAX_DIS);

	if (v1 == MAX_DIS || v2 == MAX_DIS)
	{
		v1 = -1;
		v2 = 1;
	}
	return p1 + ((p2 - p1)*(isovalue - v1)) / (v2 - v1);
}

//���㷨���������������㣬������ķ�������ʹ�ò�����㣬n=a��b��
void MCDistField::computeNormal(Vec3D &p1, Vec3D &p2, Vec3D &p3, Vec3D &normal)
{
	double v1[3] = { p1.x, p1.y, p1.z };
	double v2[3] = { p2.x, p2.y, p2.z };
	double v3[3] = { p3.x, p3.y, p3.z };
	double vc1[3], vc2[3];//��������������ƽ�е�����
	double a, b, c;
	double r;
	vc1[0] = v2[0] - v1[0]; vc1[1] = v2[1] - v1[1]; vc1[2] = v2[2] - v1[2];
	vc2[0] = v3[0] - v1[0]; vc2[1] = v3[1] - v1[1]; vc2[2] = v3[2] - v1[2];
	a = vc1[1] * vc2[2] - vc2[1] * vc1[2];
	b = vc2[0] * vc1[2] - vc1[0] * vc2[2];
	c = vc1[0] * vc2[1] - vc2[0] * vc1[1];
	r = sqrt(a * a + b * b + c * c);
	normal.x = a / r;	//��һ��
	normal.y = b / r;
	normal.z = c / r;
}

bool MCDistField::ClearGrid()
{
	DistanceField::ClearGrid();
	mesh.clear();
	meshVector.clear();

	return true;
}