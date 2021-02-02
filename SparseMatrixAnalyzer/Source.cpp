#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "Matrix_CSR.h"

int main()
{
	MatrixCSR mtxCSR;
	mtxCSR.ReadSortMtx();
	mtxCSR.ConvertMatrixMtxToCSR();
	mtxCSR.Write_crs();
	mtxCSR.CalculateProperties();
	mtxCSR.Clear();
	return 1;
}

