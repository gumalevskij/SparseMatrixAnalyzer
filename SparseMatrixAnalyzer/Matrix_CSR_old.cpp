#include "Matrix_CSR_old.h"

int MatrixCSR::ReadMtx()
{
	ifstream ofp;

	int nthreads;
	ofp.open("input.mtx");
	ofp >> sz_col;
	ofp >> sz_row;
	ofp >> sz_elem;
	nnz = sz_elem;
	rowA = new int[sz_elem];
	colA = new int[sz_elem];
	valA = new double[sz_elem];

	for (int i = 0; i < sz_elem; i++)
	{
		ofp >> colA[i];
		ofp >> rowA[i];
		ofp >> valA[i];
		rowA[i]--;
		colA[i]--;
	}

	ofp.clear();
	ofp.close();

	return 0;
}

int MatrixCSR::Write()
{
	ofstream ofp;
	ofp.open("ia.txt");
	ofp << scientific << setprecision(15);
	for (int i = 0; i < sz_row + 1; i++)
	{
		ofp << ia[i] << endl;
	}
	ofp.clear();
	ofp.close();

	ofp.open("ja.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << ja[i] << endl;
	}

	ofp.clear();
	ofp.close();
	ofp.open("aa.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << aa[i] << endl;
	}

	ofp.clear();
	ofp.close();

	return 0;
}

int MatrixCSR::ConvertMatrixMtxToCSR()
{
	ia = new int[sz_row + 1];
	ja = new int[sz_elem];
	aa = new double[sz_elem];

	for (int i = 0; i < sz_row + 1; i++)
	{
		ia[i] = 0;
	}
	ia[0] = 0;
	int currentRow = 0;
	for (int i = 0; i < sz_elem; i++)
	{
		aa[i] = valA[i];
		ja[i] = colA[i];
		ia[rowA[i] + 1] += 1;
	}

	for (int i = 0; i < sz_row; i++)
	{
		ia[i + 1] += ia[i];
	}

	return 0;
}

void MatrixCSR::Clear()
{
	delete[] rowA;
	delete[] colA;
	delete[] valA;
	delete[] ia;
	delete[] ja;
	delete[] aa;
	is_init_data = false;
}
