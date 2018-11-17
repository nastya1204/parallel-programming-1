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
	double* pMatrix;                                 // Исходная матрица
	double* pVector;                                 // Исходный вектор
	double* pResult;                                 // Результат умножения матрицы на вектор     
	int numbereleminportions;
	int ProcNum, ProcRank;
	double* PortionBuffer;                           // Буфер приема порций для MPI_Scatter
	double* ProcElem;                                // Элементы, полученные внутри процесса
	double* ProcResult;
	int* NumofElem;                                  //Массив для MPI_Gatherv
	int* Offsetofelem;                               //Массив для MPI_Gatherv
	int Size = 10;
	double tstart, tfinish, tstart1, tfinish1;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	numbereleminportions = (Size*Size) / ProcNum;    //Число элементов в порции

	//Выделение памяти
	pMatrix = new double[Size*Size];
	pVector = new double[Size];
	pResult = new double[Size];
	NumofElem = new int[ProcNum];
	Offsetofelem = new int[ProcNum];

	//Инициализация данных
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
	
	// Разделение матрицы на части
	ProcResult = new double[numbereleminportions];			//Результат на процессоре
	
	// Передача порций и вектора процессам 
	ProcElem = new double[numbereleminportions / Size];     //Массив для хранения частей результата для каждого процесса
	PortionBuffer = new double[numbereleminportions];       //Размер буфера приема порций равен числу элементов в порции
	tstart = MPI_Wtime();
	MPI_Bcast(pVector, Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ProcElem, numbereleminportions / Size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ProcResult, numbereleminportions, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(NumofElem, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Offsetofelem, ProcNum, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(pMatrix, numbereleminportions, MPI_DOUBLE, PortionBuffer, numbereleminportions, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Вычисление элементов результирующего вектора и отправление на root процесс
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
			j = 0;	  // Обнуляем индикатор движения по вектору
		}
	}

	// Сборка результата
	int n = numbereleminportions / Size;
	MPI_Gatherv(ProcElem, n, MPI_DOUBLE, pResult, NumofElem, Offsetofelem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	tfinish = MPI_Wtime();
	// Печать результата
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

	//Очистка
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