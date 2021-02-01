#pragma once

#include "stdafx.h"

struct MatrixCSR
{
	//Matrix info
	int sz_col;
	int sz_row;
	int sz_elem;

	//COO
	int nnz;
	int* colA;
	int* rowA;
	double* valA;

	//CSR
	int * crs_ia;
	int * crs_ja;
	double * crs_aa;

	//CCS
	int* ccs_ia;
	int* ccs_ja;
	double* ccs_aa;

	//properties
	bool is_init_data;
	bool is_symmetrical;

	//Graph
	vector<bool> used;
	vector<int> comp;

	int ReadMtx();
	int CCStoCRS();
	int Write_crs();
	int Write_ccs();
	int ConvertMatrixMtxToCCS();
	void Clear();
};