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

	//CCS
	int* ccs_ia;
	int* ccs_ja;
	double* ccs_aa;

	//properties
	bool is_init_data;
	bool is_symmetrical;
	double MinElement;
	double MaxElement;
	double MinModElement;
	double MaxModElement;
	double MaxDiag;
	double MinDiag;
	double MaxModDiag;
	double MinModDiag;


	//Graph
	vector<bool> used;
	vector<int> comp;
	string path_to_picture1;

	int ReadMtx();
	int ReadSortMtx();
	int CCStoCRS();
	int CRStoCCS();
	int Write_crs();
	int Write_ccs();
	int ConvertMatrixMtxToCCS();
	int ConvertMatrixMtxToCSR();
	int CalculateProperties();
	void Clear();
};
