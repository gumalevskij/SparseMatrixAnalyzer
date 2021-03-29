#pragma once
#include <chrono>
#include <map>
#include <sstream>
#include <string>
#include <algorithm>
#include "stdafx.h"
#include "lodepng.h"

class MatrixCSR
{
private:
	//Matrix info
	int sz_row;
	int sz_col;
	int sz_elem;
	string path_to_matrix;

	//COO
	int nnz;
	int* rowA;
	int* colA;
	double* valA;

	//CSR
	int * csr_ia;
	int * csr_ja;
	double * csr_aa;

	//CSR_T
	bool csr_t_exist = false;
	int* csr_t_ia;
	int* csr_t_ja;
	double* csr_t_aa;
	int* nn;
	int* reverse_row_index;

	//properties
	int Nonzeros = 0;
	bool is_init_data;
	bool is_symmetrical = true;
	bool is_structural_symmetrical = true;
	int is_symmetrical_count = 0;
	int is_structural_symmetrical_count = 0;
	double MinElement;
	double MaxElement;
	double MinModElement;
	double MaxModElement;
	double MaxDiag;
	double MinDiag;
	double MaxModDiag;
	double MinModDiag;
	double *MinRows;
	double *MaxRows;
	double *MinModRows;
	double *MaxModRows;
	bool DiagonallyDominantMatrix = true;
	bool DiagonalSignDefinite = true;
	bool ExistenceOfIsolatedSubmatrices = false;
	int countComponents;
	int countUpElements;
	int countDownElements;
	int countDiag;

	//Graph
	string path_to_picture1;
	bool *used;
	int *stack;

public:
	//Logger
	double TimeRead;
	double TimeConvertation;
	double TimeTransportation;
	double TimeWrite;
	double TimeCalculation;
	void Create_log();
	//
	int ReadMtx();
	int ReadSortMtx(string input);
	int CSRtoCSR_t();
	int Write_csr();
	int Write_csr_t();
	int ConvertMatrixMtxToCSR();
	int CheckExistenceOfIsolatedSubmatrices();
	int CalculateParameters();
	void Create_out_html();

	//portrait
	void plot(int* ptr, int* y, double* data, int n, double max, double min);
	void encodeTwoSteps(const char* filename, const unsigned char* image, unsigned width, unsigned height);

	void Clear();
};
