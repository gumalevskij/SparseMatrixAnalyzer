#pragma once
#include "Matrix_CSR.h"

int main()
{
	MatrixCSR mtxCSR;
	auto start = std::chrono::steady_clock::now();
	cout << "Start Reading" << endl;
	mtxCSR.ReadSortMtx();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	mtxCSR.TimeRead = elapsed_seconds.count();
	cout << "End Reading" << endl;

	cout << "Start Convertating" << endl;
	start = std::chrono::steady_clock::now();
	mtxCSR.ConvertMatrixMtxToCSR();
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	mtxCSR.TimeConvertation = elapsed_seconds.count();
	cout << "End Convertating" << endl;

	cout << "Start Transposing" << endl;
	start = std::chrono::steady_clock::now();
    mtxCSR.CSRtoCSR_t();
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	mtxCSR.TimeTransportation = elapsed_seconds.count();
	cout << "End Transposing" << endl;

	cout << "Start Calculating" << endl;
	start = std::chrono::steady_clock::now();
    mtxCSR.CalculateParameters();
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	mtxCSR.TimeCalculation = elapsed_seconds.count();
	cout << "End Calculating" << endl;

	cout << "Start Writing" << endl;
	start = std::chrono::steady_clock::now();
    mtxCSR.Create_out_html();
	end = std::chrono::steady_clock::now();
	elapsed_seconds = end - start;
	mtxCSR.TimeWrite = elapsed_seconds.count();
	cout << "End Writing" << endl;

	mtxCSR.Create_log();
	mtxCSR.Clear();
	return 1;
}