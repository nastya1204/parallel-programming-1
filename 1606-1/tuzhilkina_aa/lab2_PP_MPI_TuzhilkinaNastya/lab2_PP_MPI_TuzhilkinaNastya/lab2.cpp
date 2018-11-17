#include "mpi.h"
#include <math.h>  
#include <stdio.h> 
#include <stdlib.h>
#include <iostream>
#include <random>
#include <ctime>
#include <cstdlib>
#include <conio.h>

using namespace std;

void main(int argc, char* argv[])
{
	double* LinearResult;
	double* pMatrix;                                 // �������� �������
	double* pVector;                                 // �������� ������
	double* pResult;                                 // ��������� ��������� ������� �� ������     
	int numbereleminportions;
	int ProcNum, ProcRank;
	double* PortionBuffer;                           // ����� ������ ������ ��� MPI_Scatter
	double* ProcElem;                                // ��������, ���������� ������ ��������
	double* ProcResult;
	int* NumofElem;                                  //������ ��� MPI_Gatherv
	int* Offsetofelem;                               //������ ��� MPI_Gatherv
	int Size = 10;
	double tstart, tfinish, tstart1, tfinish1;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	numbereleminportions = (Size*Size) / ProcNum;    //����� ��������� � ������

	//��������� ������
	pMatrix = new double[Size*Size];
	pVector = new double[Size];
	pResult = new double[Size];
	NumofElem = new int[ProcNum];
	Offsetofelem = new int[ProcNum];

	//������������� ������
	if (ProcRank == 0)
	{
		for (int i = 0; i < Size*Size; i++)
		{
			pMatrix[i] = rand() % 20;
		}
		for (int j = 0; j < Size; j++)
		{
			pVector[j] = rand() % 20;
			pResult[j] = 0;
		}
		for (int y = 0; y < ProcNum; y++)
		{
			NumofElem[y] = numbereleminportions / Size;
			Offsetofelem[y] = y*(numbereleminportions / Size);
		}
	}
	
	// ���������� ������� �� �����
	ProcResult = new double[numbereleminportions];			//��������� �� ����������
	
	// �������� ������ � ������� ��������� 
	ProcElem = new double[numbereleminportions / Size];     //������ ��� �������� ������ ���������� ��� ������� ��������
	PortionBuffer = new double[numbereleminportions];       //������ ������ ������ ������ ����� ����� ��������� � ������
	tstart = MPI_Wtime();
	MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ProcElem, numbereleminportions / Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ProcResult, numbereleminportions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(NumofElem, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Offsetofelem, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(pMatrix, numbereleminportions, MPI_DOUBLE, PortionBuffer, numbereleminportions, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// ���������� ��������� ��������������� ������� � ����������� �� root �������
	int k = 0;
	int j = 0;
	double sumeleminrow = 0.0;
	for (int i = 0; i < numbereleminportions; i++)
	{
		ProcResult[i] = PortionBuffer[i] * pVector[j];
		j++;
		sumeleminrow = sumeleminrow + ProcResult[i];

		if (j == Size)
		{
			ProcElem[k] = sumeleminrow;
			sumeleminrow = 0.0;
			k++;
			j = 0;	  // �������� ��������� �������� �� �������
		}
	}

	// ������ ����������
	int n = numbereleminportions / Size;
	MPI_Gatherv(ProcElem, n, MPI_DOUBLE, pResult, NumofElem, Offsetofelem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	tfinish = MPI_Wtime();
	// ������ ����������
	if (ProcRank == 0)
	{
		cout << "Source Matrix:" << endl;
		for (int q = 0; q < Size*Size; q++)
		{
			if (q%Size == 0 && q != 0) cout << endl;
			cout << pMatrix[q] << " ";
		}

		cout << endl << endl << "Source Vector:" << endl;
		for (int z = 0; z < Size; z++)
		{
			cout << pVector[z] << " ";
		}

		cout << endl << endl << "Parallel Result:" << endl;
		for (int p = 0; p < Size; p++)
		{
			cout << pResult[p] << " ";
		}

		cout << endl << "ParallelTime:" << tfinish - tstart << endl;
	}

	LinearResult = new double[Size];
	if (ProcRank == 0)
	{
		double res = 0.0;
		int jj = 0;
		int ii = 0;
		tstart1 = MPI_Wtime();
		for (int i = 0; i < Size*Size; i++)
		{
			res = res + pMatrix[i] * pVector[jj];
			jj++;
			if (jj == Size)
			{
				LinearResult[ii] = res;
				ii++;
				jj = 0;
				res = 0.0;
			}
		}
		tfinish1 = MPI_Wtime();

		cout << endl << "Linear Result:" << endl;
		for (int p = 0; p < Size; p++)
		{
			cout << LinearResult[p] << " ";
		}

		cout << endl << "LinearTime:" << tfinish1 - tstart1 << endl << endl;

		system("pause");
	}

	//�������
	delete[] pMatrix;
	delete[] pVector;
	delete[] pResult;
	delete[] PortionBuffer;
	delete[] ProcElem;
	delete[] ProcResult;
	delete[] NumofElem;
	delete[] Offsetofelem;
	delete[] LinearResult;

	MPI_Finalize();
}