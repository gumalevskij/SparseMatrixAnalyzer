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
	int* ia;
	int* ja;
	double* aa;

	//properties
	bool is_init_data;
	bool is_symmetrical;

	//Graph
	vector<bool> used;
	vector<int> comp;

	int ReadMtx();
	int Write();
	int ConvertMatrixMtxToCSR();
	void Clear();
};