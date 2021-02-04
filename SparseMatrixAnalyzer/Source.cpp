#pragma once
#include "Matrix_CSR.h"

int main()
{
	MatrixCSR mtxCSR;
	mtxCSR.ReadSortMtx();
	mtxCSR.ConvertMatrixMtxToCSR();
	mtxCSR.CheckExistenceOfIsolatedSubmatrices();
	mtxCSR.Write_crs();
	//mtxCSR.WriteProperties();
	mtxCSR.Clear();
	return 1;
}

