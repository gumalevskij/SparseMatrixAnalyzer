#pragma once

#include "stdafx.h"

struct MatrixCSR
{
	//Matrix info
	int sz_row;
	int sz_col;
	int sz_elem;
	string path_to_matrix;

	//COO
	int nnz;
	int* colA;
	int* rowA;
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

	//properties
	int Nonzeros = 0;
	bool is_init_data;
	bool is_symmetrical = true;
	bool is_structural_symmetrical = true;
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


	//Graph
	string path_to_picture1;
	bool *used;
	int *stack;

	int ReadMtx();
	int ReadSortMtx();
	int CSRtoCSR_t();
	int Write_csr();
	int Write_csr_t();
	int ConvertMatrixMtxToCSR();
	int CheckExistenceOfIsolatedSubmatrices();
	int CalculateParameters();
	void Create_out_html();
	void Clear();
};
