#include "Matrix_CSR.h"
#include <map>
#include <sstream>
#include <string>

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
	//getline(ofp, line);
	//auto parameters = split(line, ' ');
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

int MatrixCSR::CSRtoCSR_t()
{
	csr_t_aa = new double[sz_elem];
	csr_t_ia = new int[sz_col + 1];
	csr_t_ja = new int[sz_elem];
	int* nn = new int[sz_col + 1];
	for (int i = 0; i <= sz_col; i++)
	{
		csr_t_ia[i] = 0;
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
		csr_t_ia[csr_ja[i] + 1]++;
	}
	for (int i = 1; i <= sz_col; i++)
	{
		csr_t_ia[i] += csr_t_ia[i - 1];
	}
	std::copy(csr_t_ia, csr_t_ia + sz_col + 1, nn);

	for (int i = 0; i < sz_elem; i++)
	{
		int x = nn[csr_ja[i]];
		nn[csr_ja[i]] += 1;
		csr_t_aa[x] = csr_aa[i];
		csr_t_ja[x] = reverse_row_index[i];
	}
	return 1;
}

int MatrixCSR::Write_csr()
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

int MatrixCSR::Write_csr_t()
{
	ofstream ofp;
	ofp.open("ia_t.txt");
	ofp << scientific << setprecision(15);
	for (int i = 0; i <= sz_col; i++)
	{
		ofp << csr_t_ia[i] << endl;
	}
	ofp.clear();
	ofp.close();

	ofp.open("ja_t.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << csr_t_ja[i] << endl;
	}

	ofp.clear();
	ofp.close();
	ofp.open("aa_t.txt");
	for (int i = 0; i < sz_elem; i++)
	{
		ofp << csr_t_aa[i] << endl;
	}

	ofp.clear();
	ofp.close();

	return 0;
}

int MatrixCSR::WriteProperties()
{
	ofstream ofp;
	ofp.open("prop.txt");
	ofp << scientific << setprecision(15);
	ofp << ExistenceOfIsolatedSubmatrices << endl;
	ofp << currentComp << endl;
	ofp.clear();
	ofp.close();

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

int MatrixCSR::CheckExistenceOfIsolatedSubmatrices()
{
	used = new bool[sz_elem];
	stack = new int[sz_elem];

	for (int i = 0; i < sz_elem; i++)
	{
		used[i] = false;
	}

	currentComp = 0;
	int currentStack = -1;
	int currentEl = -1;
	for (int i = 0; i < sz_row; i++)
	{
		if (!used[i])
		{
			currentStack++;
			stack[currentStack] = i;
			used[i] = true;

			while (currentStack != -1)
			{
				currentEl = stack[currentStack];
				currentStack--;
				for (int j = csr_ia[currentEl]; j < csr_ia[currentEl + 1]; j++)
				{
					if (!used[csr_ja[j]])
					{
						used[csr_ja[j]] = true;
						currentStack++;
						stack[currentStack] = csr_ja[j];
					}
				}
			}
			currentComp++;
		}
	}

	if (currentComp == 1)
	{
		ExistenceOfIsolatedSubmatrices = true;
	}

	delete[] used;
	delete[] stack;
	return 0;
}

void MatrixCSR::Create_out_html()
{
  
  ofstream stream1((path_to_matrix+".out.html").c_str());
  stream1 << "<h1 style=\"color: #5e9ca0;\">Matrix characteristics for " <<  path_to_matrix << " </h1> " << endl;
  stream1 << "<h2><strong><span class=\"VIiyi\" lang=\"en\"><span class=\"JLqJ4b ChMk0b\" data-language-for-alternatives=\"en\" data-language-to-translate-into=\"ru\" data-phrase-index=\"0\">Matrix portrait:</span></span></strong></h2> " << endl;
  stream1 <<  "<p><strong><span class=\"VIiyi\" lang=\"en\"><span class=\"JLqJ4b ChMk0b\" data-language-for-alternatives=\"en\" data-language-to-translate-into=\"ru\" data-phrase-index=\"0\"><img src=\"" <<  path_to_picture1 << "\" alt=\"\" /></span></span></strong></p>" << endl;
  stream1 << "<p>&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;</p> " << endl;
  stream1 << "<h2 style=\"color: #2e6c80;\">Matrix details:</h2>" << endl; 
  stream1 << "<table class=\"editorDemoTable\" style=\"width: 497px; border-style: solid; border-color: blue; float: left; height: 393px;\">" << endl;

  std::vector<std::string> names;
  names.push_back("Number of rows  ");   names.push_back(to_string(sz_row));
  names.push_back("Number of columns  "); names.push_back(to_string(sz_col));
  names.push_back("Int  "); 
  if(is_init_data)names.push_back("Yes");
  else names.push_back("No");
  names.push_back("Symmetrical ");
  if(is_symmetrical)names.push_back("Yes");
  else names.push_back("No");
  names.push_back("MinElement  ");   names.push_back(to_string(MinElement));
  names.push_back("MaxElement  ");  names.push_back(to_string(MaxElement));
  names.push_back("MinModElement  ");  names.push_back(to_string(MinModElement));
  names.push_back("MaxModElement  ");  names.push_back(to_string(MaxModElement));
  names.push_back("MaxDiag  ");      names.push_back(to_string(MaxDiag));
  names.push_back("MinDiag  ");     names.push_back(to_string(MinDiag));
  names.push_back("MaxModDiag  ");  names.push_back(to_string(MaxModDiag));
  names.push_back("MinModDiag  ");  names.push_back(to_string(MinModDiag));

  for(int i = 0; i < names.size(); i = i +2){
    stream1 << " <tr style=\"height: 36px;\"> " << endl;
    stream1 << " <td style=\"width: 207.783px; height: 36px;\"> " << names[i] << "        </td> " << endl;
    stream1 << " <td style=\"width: 272.517px; height: 36px;\"> " <<  names[i+1] << "</td> </tr> " << endl;
  }
  stream1 << "</tbody> </table> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p><p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> " << endl;
  
  
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
