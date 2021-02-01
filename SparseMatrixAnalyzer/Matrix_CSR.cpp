#include "Matrix_CSR.h"
#include <map>

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
		ofp >> rowA[i];
		ofp >> colA[i];
		ofp >> valA[i];
		rowA[i]--;
		colA[i]--;
	}

	ofp.clear();
	ofp.close();

	return 0;
}

int MatrixCSR::CCStoCRS()
{
	crs_aa = new double[sz_elem];
	crs_ja = new int[sz_elem];
	crs_ia = new int[sz_row + 1];
	int* nn = new int[sz_row + 1];
	for (int i = 0; i <= sz_row; i++)
	{
		crs_ia[i] = 0;
	}
	int k = 0;
	int *reverse_col_index = new int[sz_elem]; // each value belongs to which col
	for (int i = 0; i < sz_col; i++)
	{
		for (int j = 0; j < ccs_ja[i + 1] - ccs_ja[i]; j++)
		{
			reverse_col_index[k] = i;
			k++;
		}
	}

	for (int i = 0; i < sz_elem; i++)
	{
		crs_ia[ccs_ia[i] + 1]++;
	}
	for (int i = 1; i <= sz_row; i++)
	{
		crs_ia[i] += crs_ia[i - 1];
	}
	std::copy(crs_ia, crs_ia + sz_row+1, nn);

	for (int i = 0; i < sz_elem; i++) 
	{
		int x = nn[ccs_ia[i]];
		nn[ccs_ia[i]] += 1;
		crs_aa[x] = ccs_aa[i];
		crs_ja[x] = reverse_col_index[i];
	}
	return 1;
}

int MatrixCSR::Write_crs()
{
	ofstream ofp;
	ofp.open("ia.txt");
	ofp << scientific << setprecision(15);
	for (int i = 0; i <= sz_row; i++)
	{
		ofp << crs_ia[i] << endl;
	}
	ofp.clear();
	ofp.close();

	ofp.open("ja.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << crs_ja[i] << endl;
	}
	
	ofp.clear();
	ofp.close();
	ofp.open("aa.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << crs_aa[i] << endl;
	}

	ofp.clear();
	ofp.close();

	return 0;
}

int MatrixCSR::Write_ccs()
{
	ofstream ofp;
	ofp.open("ccs_ia.txt");
	ofp << scientific << setprecision(15);
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << ccs_ia[i] << endl;
	}
	ofp.clear();
	ofp.close();

	ofp.open("ccs_ja.txt");
	for (int i = 0; i <= sz_col; i++)
	{
		ofp << ccs_ja[i] << endl;
	}

	ofp.clear();
	ofp.close();
	ofp.open("ccs_aa.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << ccs_aa[i] << endl;
	}

	ofp.clear();
	ofp.close();

	return 0;
}

int MatrixCSR::ConvertMatrixMtxToCCS()
{
	ccs_ja = new int[sz_col];
	ccs_ia = new int[sz_elem];
	ccs_aa = new double[sz_elem];

	for (int i = 0; i < sz_col + 1; i++)
	{
		ccs_ja[i] = 0;
	}

	ccs_ja[0] = 0;
	int currentRow = 0;
	for (int i = 0; i < sz_elem; i++)
	{
		ccs_aa[i] = valA[i];
		ccs_ia[i] = rowA[i];
		ccs_ja[colA[i]+1] += 1;
	}

	for (int i = 0; i < sz_row; i++)
	{
		ccs_ja[i+1] += ccs_ja[i];
	}

	return 0;
}

void MatrixCSR::Clear()
{
	delete [] rowA;
	delete [] colA;
	delete [] valA;
	delete [] crs_ia;
	delete [] crs_ja;
	delete [] crs_aa;
	delete [] ccs_ia;
	delete [] ccs_ja;
	delete [] ccs_aa;
	is_init_data = false;
}