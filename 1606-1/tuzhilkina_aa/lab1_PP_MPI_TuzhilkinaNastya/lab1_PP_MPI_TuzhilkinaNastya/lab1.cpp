#include <mpi.h>
#include <iostream>
#include <cstdlib>

using namespace std;

int* Creation_of_the_matrix(int row, int column)
{
	int *matrix;
		matrix = new int[row*column];
	return matrix;
}

int* Creation_of_the_vector(int rows)
{
	int *vector;
	vector = new int[rows];
	return vector;
}

void Filling_of_the_matrix(int* matrix, int row, int column) 
{
	for (int i = 0; i < row; i++)
		for (int j = 0; j < column; j++)
			matrix[i * row + j] = rand() % 500;
}

void Show_the_matrix(int* matrix, int row, int column)
{
	cout << "Matrix:" << endl;
	if (row <= 20 && column <= 20) 
	{
		for (int i = 0; i < row; i++) 
		{
			for (int j = 0; j < column; j++)
				cout << matrix[i * row + j] << " ";
			cout << endl;
		}
	}
}



int main(int argc, char **argv) 
{
	int rows, columns;
	int rank, size;
	int *matrix = NULL, *vector = NULL, *resultPar = NULL, *mas = NULL;
	double endSeq = 0;
	double startSeq = 0;
	double endPar = 0;
	double startPar = 0;
	double timeSeq = 0;
	double timePar = 0;

	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) 
	{
		cout << "Enter the number of rows=";
		cin >> rows;
			while (rows == 0)
			{
				cout << "The number of rows cannot be 0. Enter the number of rows=";
				cin >> rows;
			}

		cout << "Enter the number of columns=";
		cin >> columns;
		while (columns == 0)
		{
			cout << "The number of columns cannot be 0. Enter the number of columns=";
			cin >> columns;
		}
		cout << endl;

		matrix = Creation_of_the_matrix(rows, columns);
		vector = Creation_of_the_vector(rows);
		resultPar = new int[rows];
		for (int i = 0; i < rows; i++)
			resultPar[i] = INT_MIN;

		Filling_of_the_matrix(matrix, rows, columns);

		if (rows <= 20 && columns <= 20)
			Show_the_matrix(matrix, rows, columns);

		/*Последовательный алгоритм*/
		startSeq = MPI_Wtime();

		cout << endl << "  ---SEQUENCE VERSION---  " << endl;

		for (int i = 0; i < rows; i++) {
			vector[i] = INT_MIN;
			for (int j = 0; j < columns; j++)
			{
				if (matrix[i * columns + j] > vector[i])
					vector[i] = matrix[i * columns + j];
			}
		}

		endSeq = MPI_Wtime();
		timeSeq = (endSeq - startSeq) * 1000;
		cout << "Time of sequence version: " << timeSeq << " ms" << endl;

		cout << "Result: ";
		for (int i = 0; i < rows; i++)
			cout << endl << " Maximum value in row " << i + 1 << " is " << vector[i] << " ";
		cout << endl;

		/*Паралелльный алгоритм*/
		startPar = MPI_Wtime();

		cout << endl << "  ---PARALLEL VERSION---  " << endl;

		for (int i = 1; i < size; i++) {
			MPI_Send(&columns, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&rows, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(matrix, columns*rows, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		for (int i = 1; i < size; i++) {
			MPI_Recv(vector, rows, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
			for (int j = 0; j < rows; j++) {
				if (vector[j] > resultPar[j])
					resultPar[j] = vector[j];
			}
		}

		endPar = MPI_Wtime();

		timePar = (endPar - startPar);
		cout << "Time of parallel version: " << timePar << " ms" << endl;

		cout << "Result: ";
		for (int i = 0; i < rows; i++)
			cout << endl << " Maximum value in row " << i + 1 << " is " << resultPar[i] << " ";
		cout << endl;

		// Сравнение времени работы алгоритмов
		if (timePar <= timeSeq)
			cout << endl << "Parallel version faster, then sequence" << endl;
		else
			cout << "Sequence version faster, then parallel" << endl;

		delete matrix;
		delete vector;
		delete resultPar;
	}
	
	else 
	{
		MPI_Recv(&columns, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
		MPI_Recv(&rows, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

		matrix = Creation_of_the_matrix(rows, columns);
		vector = Creation_of_the_vector(rows);

		MPI_Recv(matrix, columns*rows, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);

		for (int i = 0; i < rows; i++) 
		{
			vector[i] = INT_MIN;
			for (int j = (rank - 1); j < columns; j += (size - 1))
				if (matrix[i * columns + j] > vector[i])
					vector[i] = matrix[i * columns + j];
		}

		MPI_Send(vector, rows, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	system("pause");
	return 0;
}