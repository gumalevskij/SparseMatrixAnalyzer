#include "Matrix_CSR.h"
#include <map>
#include <sstream>

using namespace std;

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

static std::vector<std::string> split(const std::string& s, char delim)
{
	std::vector<std::string> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

int MatrixCSR::ReadSortMtx()
{
	ifstream ofp;
	string line, line2;
	int nthreads;
	ofp.open("input.mtx");
	getline(ofp, line);
	auto parameters = split(line, ' ');
	while(getline(ofp, line2))
	{
		if (line2[0] != '%')
		{
			break;
		}
	}
	auto sizes = split(line2, ' ');
	sz_col = atoi(sizes[0].c_str());
	sz_row = atoi(sizes[1].c_str());
	sz_elem = atoi(sizes[2].c_str());

	nnz = sz_elem;
	rowA = new int[sz_elem];
	colA = new int[sz_elem];
	valA = new double[sz_elem];
	map<int, map<int, double>> elements;
	map<int, map<int, double>>::iterator it;

	for (int i = 0; i < sz_elem; i++)
	{
		ofp >> rowA[i];
		ofp >> colA[i];
		ofp >> valA[i];
		rowA[i]--;
		colA[i]--;

		it = elements.find(rowA[i]);
		if (it != elements.end())
		{
			elements[rowA[i]][colA[i]] = valA[i];
		}
		else
		{
			map<int, double> ColandVal;
			ColandVal[colA[i]] = valA[i];
			elements[rowA[i]] = ColandVal;
		}
	}

	int i = 0;
	for(auto row : elements)
	{
		for(auto col : row.second)
		{
			rowA[i] = row.first;
			colA[i] = col.first;
			valA[i] = col.second;
			i++;
		}
	}

	ofp.clear();
	ofp.close();

	return 0;
}

int MatrixCSR::CCStoCRS()
{
	csr_aa = new double[sz_elem];
	csr_ja = new int[sz_elem];
	csr_ia = new int[sz_row + 1];
	int* nn = new int[sz_row + 1];
	for (int i = 0; i <= sz_row; i++)
	{
		csr_ia[i] = 0;
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
		csr_ia[ccs_ia[i] + 1]++;
	}
	for (int i = 1; i <= sz_row; i++)
	{
		csr_ia[i] += csr_ia[i - 1];
	}
	std::copy(csr_ia, csr_ia + sz_row+1, nn);

	for (int i = 0; i < sz_elem; i++) 
	{
		int x = nn[ccs_ia[i]];
		nn[ccs_ia[i]] += 1;
		csr_aa[x] = ccs_aa[i];
		csr_ja[x] = reverse_col_index[i];
	}
	return 1;
}

int MatrixCSR::CRStoCCS()
{
	ccs_aa = new double[sz_elem];
	ccs_ja = new int[sz_col + 1];
	ccs_ia = new int[sz_elem];
	int* nn = new int[sz_col + 1];
	for (int i = 0; i <= sz_col; i++)
	{
		ccs_ja[i] = 0;
	}
	int k = 0;
	int* reverse_row_index = new int[sz_elem]; // each value belongs to which col
	for (int i = 0; i < sz_row; i++)
	{
		for (int j = 0; j < csr_ia[i + 1] - csr_ia[i]; j++)
		{
			reverse_row_index[k] = i;
			k++;
		}
	}

	for (int i = 0; i < sz_elem; i++)
	{
		ccs_ja[csr_ja[i] + 1]++;
	}
	for (int i = 1; i <= sz_col; i++)
	{
		ccs_ja[i] += ccs_ja[i - 1];
	}
	std::copy(ccs_ja, ccs_ja + sz_col + 1, nn);

	for (int i = 0; i < sz_elem; i++)
	{
		int x = nn[csr_ja[i]];
		nn[csr_ja[i]] += 1;
		ccs_aa[x] = csr_aa[i];
		ccs_ia[x] = reverse_row_index[i];
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
		ofp << csr_ia[i] << endl;
	}
	ofp.clear();
	ofp.close();

	ofp.open("ja.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << csr_ja[i] << endl;
	}
	
	ofp.clear();
	ofp.close();
	ofp.open("aa.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << csr_aa[i] << endl;
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

int MatrixCSR::ConvertMatrixMtxToCSR()
{
	csr_ia = new int[sz_row + 1];
	csr_ja = new int[sz_elem];
	csr_aa = new double[sz_elem];

	for (int i = 0; i < sz_row + 1; i++)
	{
		csr_ia[i] = 0;
	}
	csr_ia[0] = 0;
	int currentRow = 0;
	for (int i = 0; i < sz_elem; i++)
	{
		csr_aa[i] = valA[i];
		csr_ja[i] = colA[i];
		csr_ia[rowA[i] + 1] += 1;
	}

	for (int i = 0; i < sz_row; i++)
	{
		csr_ia[i + 1] += csr_ia[i];
	}

	return 0;
}

int MatrixCSR::CalculateProperties()
{
	for (int i = 0; i < sz_row; i++)
	{

	}

	return 0;
}

void MatrixCSR::Clear()
{
	delete [] rowA;
	delete [] colA;
	delete [] valA;
	delete [] csr_ia;
	delete [] csr_ja;
	delete [] csr_aa;
	is_init_data = false;
}