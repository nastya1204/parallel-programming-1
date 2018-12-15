#define _USE_MATH_DEFINES
#define STB_IMAGE_IMPLEMENTATION
#define SIGMA 2

#include <stdint.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

typedef unsigned char uint8_t;

void generateKernel(double* kernel, int sigma) 
{
	for (int i = -1; i <= 1; i++) 
	{
		for (int j = -1; j <= 1; j++) 
		{
			kernel[i + j + 2 * (i + 2)] = (double)(exp(-(i*i + j*j) / (sigma*sigma))) / (sigma*sqrt(2 * M_PI));
		}
	}

}

void countElemProc(int* elem_proc, int* displs, int height, int width, int psize) 
{
	int port = height / psize;
	port = (port < 1) ? 1 : port;
	int rem = height % psize;
	int elem = port*width;

	int sent_rows_num = 0;
	for (int i = 0; i < psize; ++i) 
	{
		elem_proc[i] = elem;
		if (rem > 0) 
		{
			elem_proc[i] += width;
			rem--;
		}

		displs[i] = sent_rows_num;
		sent_rows_num += elem_proc[i];
	}
}

uint8_t* imageProcess(uint8_t* recvbuf, double* kernel, int width, int elem_proc) 
{
	uint8_t* res = new uint8_t[elem_proc];
	for (int i = width; i < elem_proc - width; i++) 
	{
		if (((i % width) != (width - 1)) && ((i % width) != 0)) 
		{
			int k = 0;
			for (int j = -width; j <= width; j += width) 
			{
				res[i] += static_cast<uint8_t>(recvbuf[(i - 1) + j] * kernel[k] + recvbuf[i + j] * kernel[k + 1] + recvbuf[(i + 1) + j] * kernel[k + 2]);
				k += 3;
			}
		}
		else 
		{
			res[i] = recvbuf[i];
		}
	}
	return res;
}

int checkEquality(uint8_t* seq_image, uint8_t* par_image, int height, int width) 
{
	for (int i = 0; i < height * width; ++i) 
	{
		if (seq_image[i] != par_image[i]) 
		{
			printf("%d\n", i);
			return 0;
		}
	}
	return 1;
}

void main(int argc, char* argv[]) 
{
	int ProcRank;
	int psize = 0;
	int width = 0, height = 0;
	int port = 0;
	int offset = 0;
	double rem = 0;
	double* kernel = new double[9];
	double start_par_time = 0, end_par_time = 0;
	double start_seq_time = 0, end_seq_time = 0;
	uint8_t* par_image = nullptr;
	uint8_t* seq_image = nullptr;

	if (argc < 3) 
	{
		printf("Image size not set\n");
		return;
	}
	else 
	{
		if (atoi(argv[1]) < 100 && atoi(argv[2]) < 100)
		{
			printf("Invalid arguments\n");
			return;
		}
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &psize);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);


	double* buf = new double[11];
	if (ProcRank == 0) 
	{
		width = atoi(argv[1]);
		height = atoi(argv[2]);

		generateKernel(kernel, SIGMA);
		par_image = new uint8_t[width*height]; 
		for (int i = 0; i < width*height; i++) 
		{
			par_image[i] = rand() % 255;
		}

		buf[0] = width;
		buf[1] = height;
		for (int i = 2; i < 11; i++) 
		{
			buf[i] = kernel[i - 2];
		}
	}
	 
	//Sequential algorithm
	 

	start_seq_time = MPI_Wtime();
	if (ProcRank == 0)
		seq_image = imageProcess(par_image, kernel, width, height*width);
	end_seq_time = MPI_Wtime();

	 
	//Prepare additional date
	 

	start_par_time = MPI_Wtime();

	MPI_Bcast(buf, 11, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (ProcRank != 0) 
	{
		width = buf[0];
		height = buf[1];
		for (int i = 2; i < 11; i++) 
		{
			kernel[i - 2] = buf[i];
		}
		offset = width;
	}

	port = height / psize;
	port = (port < 1) ? 1 : port;
	rem = height - port * psize;      //remainder

	int* elem_proc = new int[psize];  // count elements for one proc
	int* displs = new int[psize];     // displ for each portion
	int* sizerec = new int[psize];    // size recv buf


	countElemProc(elem_proc, displs, height, width, psize);

	 
	//Image division
	 
	uint8_t* recvbuf = nullptr;
	if ((ProcRank == 0) || (ProcRank == psize - 1)) 
	{
		sizerec[ProcRank] = elem_proc[ProcRank] + width;
		recvbuf = new uint8_t[sizerec[ProcRank]];
	}
	else 
	{
		sizerec[ProcRank] = elem_proc[ProcRank] + 2 * width;
		recvbuf = new uint8_t[sizerec[ProcRank]];
	}
	for (int i = 0; i < elem_proc[ProcRank]; i++)
	{
		recvbuf[i] = 0;
	}

	MPI_Scatterv(par_image, elem_proc, displs, MPI_UINT8_T, recvbuf + offset, elem_proc[ProcRank], MPI_UINT8_T, 0, MPI_COMM_WORLD);

	 
	//Block edges transfer
	 
	MPI_Status status;

	if (ProcRank != psize - 1) 
	{
		MPI_Recv(recvbuf + elem_proc[ProcRank] + offset, width, MPI_UINT8_T, ProcRank + 1, 0, MPI_COMM_WORLD, &status);
	}
	if (ProcRank != 0) 
	{
		MPI_Send(recvbuf + offset, width, MPI_UINT8_T, ProcRank - 1, 0, MPI_COMM_WORLD);
	}

	if (ProcRank != 0) 
	{
		MPI_Recv(recvbuf, width, MPI_UINT8_T, ProcRank - 1, 0, MPI_COMM_WORLD, &status);
	}
	if (ProcRank != psize - 1) 
	{
		MPI_Send(recvbuf + elem_proc[ProcRank] - width + offset, width, MPI_UINT8_T, ProcRank + 1, 0, MPI_COMM_WORLD);
	}

	 
	//Image processing and build new image
	 

	recvbuf = imageProcess(recvbuf, kernel, width, sizerec[ProcRank]);

	MPI_Gatherv((void*)(recvbuf + offset), elem_proc[ProcRank], MPI_UINT8_T, par_image, elem_proc, displs, MPI_UINT8_T, 0, MPI_COMM_WORLD);

	end_par_time = MPI_Wtime();

	 
	//Correctness
	 

	int res = 0;
	if (ProcRank == 0) 
	{
		res = checkEquality(seq_image, par_image, height, width);
	}

	if (ProcRank == 0) 
	{
		delete[] par_image;
		delete[] seq_image;
		printf("time_par = %f\n", end_par_time - start_par_time);
		printf("time_seq = %f\n", end_seq_time - start_seq_time);
		if (res)
			printf("The algorithm is correct.\n");
	}
	delete[] kernel;
	delete[] buf;

	MPI_Finalize();
}