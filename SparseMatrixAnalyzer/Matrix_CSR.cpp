#include "Matrix_CSR.h"

using namespace std;

void MatrixCSR::Create_log()
{
	ofstream ofp;
	ofp.open("Log.txt");
	ofp << "Time of Read  " << TimeRead << endl;
	ofp << "Time of Convertation " << TimeConvertation << endl;
	ofp << "Time of Transportation " << TimeTransportation << endl;
	ofp << "Time of Calculation " << TimeCalculation << endl;
	ofp << "Time of Write " << TimeWrite << endl;
	
	ofp.clear();
	ofp.close();
}

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

int MatrixCSR::CSRtoCSR_t()
{
	csr_t_exist = true;
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

	countComponents = 0;
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

				for (int j = csr_t_ia[currentEl]; j < csr_t_ia[currentEl + 1]; j++)
				{
					if (!used[csr_t_ja[j]])
					{
						used[csr_t_ja[j]] = true;
						currentStack++;
						stack[currentStack] = csr_t_ja[j];
					}
				}
			}
			countComponents++;
		}
	}

	if (countComponents == 1)
	{
		ExistenceOfIsolatedSubmatrices = true;
	}

	delete[] used;
	delete[] stack;
	return 0;
}

int MatrixCSR::CalculateParameters()
{
	CheckExistenceOfIsolatedSubmatrices();
	double MinRow;
	double MaxRow;
	double MinModRow;
	double MaxModRow;
	MinRows = new double[sz_row];
	MaxRows = new double[sz_row];
	MinModRows = new double[sz_row];
	MaxModRows = new double[sz_row];

	MinElement = DBL_MAX;
	MaxElement = -DBL_MAX;
	MinModElement = DBL_MAX;
	MaxModElement = -DBL_MAX;
	MinDiag = DBL_MAX;
	MaxDiag = -DBL_MAX;
	MinModDiag = DBL_MAX;
	MaxModDiag = -DBL_MAX;
	countUpElements = 0;
	countDownElements = 0;
	countDiag = 0;

	if (sz_row != sz_col)
	{
		is_symmetrical = false;
		is_structural_symmetrical = false;
	}

	for (int i = 0; i < sz_row; i++)
	{
		MinRow = DBL_MAX;;
		MaxRow = -DBL_MAX;
		MinModRow = DBL_MAX;
		MaxModRow = -DBL_MAX;
		double AbsDiag = 0;
		double SumAbsRows = 0;

		if (csr_ia[i + 1] != csr_t_ia[i + 1])
		{
			is_symmetrical = false;
			is_structural_symmetrical = false;
		}

		for (int j = csr_ia[i]; j < csr_ia[i + 1]; j++)
		{
			double val = csr_aa[j];
			double modVal = fabs(csr_aa[j]);

			if (i == csr_ja[j])
			{
				countDiag++;
				AbsDiag = modVal;
			}
			else
			{
				SumAbsRows += modVal;
				if (i < csr_ja[j])
				{
					countUpElements++;
				}
				else
				{
					countDownElements++;
				}
			}
			if (val != 0)
			{
				Nonzeros++;

				if (val < MinRow)
				{
					MinRow = val;
					if (val < MinElement)
					{
						MinElement = val;
					}
				}
				else
				{
					if (val > MaxRow)
					{
						MaxRow = val;

						if (val > MaxElement)
						{
							MaxElement = val;
						}
					}
				}

				if (modVal < MinModRow)
				{
					MinModRow = modVal;

					if (modVal < MinModElement)
					{
						MinModElement = modVal;
					}
				}
				else
				{
					if (modVal > MaxModRow)
					{
						MaxModRow = modVal;

						if (modVal > MaxModElement)
						{
							MaxModElement = modVal;
						}
					}
				}

				if (i == csr_ja[j])
				{
					if (val < MinDiag)
					{
						MinDiag = val;
					}
					else
					{
						if (val > MaxDiag)
						{
							MaxDiag = val;
						}
					}

					if (modVal < MinModDiag)
					{
						MinModDiag = modVal;
					}
					else
					{
						if (modVal > MaxModDiag)
						{
							MaxModDiag = modVal;
						}
					}
				}

				
			}

			if (csr_ja[j] != csr_t_ja[j])
			{
				is_symmetrical = false;
				is_structural_symmetrical = false;
			}
			else
			{
				if (csr_aa[j] != csr_t_aa[j])
				{
					is_structural_symmetrical = false;
				}
			}
		}

		if ((AbsDiag == 0 && SumAbsRows > 0) || AbsDiag < SumAbsRows)
		{
			DiagonallyDominantMatrix = false;
		}

		MinRows[i] = MinRow;
		MaxRows[i] = MaxRow;
		MinModRows[i] = MinModRow;
		MaxModRows[i] = MaxModRow;
	}
	 
	return 0;
}

void MatrixCSR::Create_out_html()
{
	path_to_matrix = "out";
	path_to_picture1 = "portrait.png";
	plot(csr_ia, csr_ja, csr_aa, sz_row, MaxElement, MinElement);
	char buffer[255];
	
	ofstream stream1((path_to_matrix+".out.html").c_str());
	stream1 << scientific << setprecision(15);
	stream1 << "<h1 style=\"color: #5e9ca0;\">Matrix characteristics for " <<  path_to_matrix << " </h1> " << endl;
	stream1 << "<h2><strong><span class=\"VIiyi\" lang=\"en\"><span class=\"JLqJ4b ChMk0b\" data-language-for-alternatives=\"en\" data-language-to-translate-into=\"ru\" data-phrase-index=\"0\">Matrix portrait:</span></span></strong></h2> " << endl;
	stream1 <<  "<p><strong><span class=\"VIiyi\" lang=\"en\"><span class=\"JLqJ4b ChMk0b\" data-language-for-alternatives=\"en\" data-language-to-translate-into=\"ru\" data-phrase-index=\"0\"><img src=\"" <<  path_to_picture1 << "\" width=\"400\" height=\"400\" alt=\"\" /></span></span></strong></p>" << endl;
	stream1 << "<p>&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;</p> " << endl;
	stream1 << "<h2 style=\"color: #2e6c80;\">Matrix details:</h2>" << endl;
	stream1 << "<table class=\"editorDemoTable\" style=\"width: 497px; border-style: solid; border-color: blue; float: left; height: 393px;\">" << endl;

	std::vector<std::string> names;
	int namesSize = 0;
	names.push_back("Num Rows  ");   names.push_back(to_string(sz_row));
	names.push_back("Num Cols  "); names.push_back(to_string(sz_col));
	names.push_back("Nonzeros   ");   names.push_back(to_string(Nonzeros));
	names.push_back("Pattern Entries  "); names.push_back(to_string(sz_elem));
	names.push_back("Pattern Symmetry");
	if (is_structural_symmetrical)names.push_back("Yes");
	else names.push_back("No");

	names.push_back("Numeric Symmetry");
	if(is_symmetrical)names.push_back("Yes");
	else names.push_back("No");

	names.push_back("Diagonally Dominant");
	if (DiagonallyDominantMatrix) names.push_back("Yes");
	else names.push_back("No");

	names.push_back("Num Components  ");  names.push_back(to_string(countComponents));
	names.push_back("Num Upper triangular  ");  names.push_back(to_string(countUpElements));
	names.push_back("Num Lower triangular  ");  names.push_back(to_string(countDownElements));
	names.push_back("Num Diag  ");  names.push_back(to_string(countDiag));
	namesSize = names.size();

	for (int i = 0; i < names.size(); i = i + 2)
	{
		stream1 << " <tr style=\"height: 36px;\"> " << endl;
		stream1 << " <td style=\"width: 207.783px; height: 36px;\"> " << names[i] << "        </td> " << endl;
		stream1 << " <td style=\"width: 272.517px; height: 36px;\"> " << names[i + 1] << "</td> </tr> " << endl;
	}

	sprintf_s(buffer, "%.2e", MinElement);
	names.push_back("Min/Max Value  ");   names.push_back(buffer);
	sprintf_s(buffer, "%.2e", MaxElement);
	names.push_back(buffer);
	sprintf_s(buffer, "%.2e", MinModElement);
	names.push_back("Min/Max Abs Value  ");  names.push_back(buffer);
	sprintf_s(buffer, "%.2e", MaxModElement);
    names.push_back(buffer);
	sprintf_s(buffer, "%.2e", MinDiag);
	names.push_back("Min/Max Diag  ");     names.push_back(buffer);
	sprintf_s(buffer, "%.2e", MaxDiag);
	names.push_back(buffer);
	
	sprintf_s(buffer, "%.2e", MinModDiag);
	names.push_back("Min/Max Abs Diag  ");  names.push_back(buffer);
	sprintf_s(buffer, "%.2e", MaxModDiag);
	names.push_back(buffer);
	
	for(int i = namesSize; i < names.size(); i = i + 3)
	{
		stream1 << " <tr style=\"height: 36px;\"> " << endl;
		stream1 << " <td style=\"width: 207.783px; height: 36px;\"> " << names[i] << "        </td> " << endl;
		stream1 << " <td style=\"width: 272.517px; height: 36px;\"> " <<  names[i+1] <<" / "<< names[i + 2] << "</td> </tr> " << endl;
	}
	stream1 << "</tbody> </table> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p><p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> <p>&nbsp;</p> " << endl;
}

typedef unsigned int uint;

void MatrixCSR::plot(int* ptr, int* Y, double* data, int n, double max, double min)
{
	const char* filename = "portrait.png";

	unsigned int x, y,
		mod = 1;
	int size = n,
		i = 0;

	if (size > 1000) {
		n = 1000;
		mod = size / 1000;
	}
	unsigned int width = n, height = n;
	unsigned char* image = new unsigned char[width * height * 4];
	sort(csr_t_aa, csr_t_aa + sz_elem);
	double AltitudeTop = 255 * (csr_t_aa[(int)(sz_elem/3)] - min) / (max - min);
	double AltitudeDown = 255 * (csr_t_aa[(int)(sz_elem*(2/3))] - min) / (max - min);
	for (int row = 0; row < size; row++) 
	{
		for (int k = 0; k < (ptr[row + 1] - ptr[row]); k++)
		{
			x = Y[i] / mod;
			y = (row) / mod;
			if (x < width && y < height)
			{
				image[4 * width * y + 4 * x + 0] = 0;
				image[4 * width * y + 4 * x + 1] = 0;
				image[4 * width * y + 4 * x + 2] = 0;
				image[4 * width * y + 4 * x + 3] = 255;

				double Amplitude = 255 * (csr_aa[csr_ja[k]] - min) / (max - min);

				if (Amplitude > AltitudeTop)
				{
					if (image[4 * width * y + 4 * x + 2] < Amplitude)
					{
						image[4 * width * y + 4 * x + 2] = Amplitude;
					}
				}
				else
				{
					if (Amplitude > AltitudeDown)
					{
						if (image[4 * width * y + 4 * x + 1] < Amplitude)
						{
							image[4 * width * y + 4 * x + 1] = Amplitude;
						}
					}
					else
					{
						if (image[4 * width * y + 4 * x + 0] < Amplitude)
						{
							image[4 * width * y + 4 * x + 0] = Amplitude;
						}
					}
				}

				i++;
			}
		}
	}

	encodeTwoSteps(filename, image, width, height);
	delete[] image;
}

void MatrixCSR::encodeTwoSteps(const char* filename, const unsigned char* image, unsigned width, unsigned height) 
{
	unsigned char* png;
	size_t pngsize;

	unsigned error = lodepng_encode32(&png, &pngsize, image, width, height);
	if (!error) lodepng_save_file(png, pngsize, filename);

	/*if there's an error, display it*/
	if (error) printf("error %u: %s\n", error, lodepng_error_text(error));

	free(png);
}

void MatrixCSR::Clear()
{
	delete [] rowA;
	delete [] colA;
	delete [] valA;
	delete [] csr_ia;
	delete [] csr_ja;
	delete [] csr_aa;
	delete[] MinRows;
	delete[] MaxRows;
	delete[] MinModRows;
	delete[] MaxModRows;
	is_init_data = false;
}