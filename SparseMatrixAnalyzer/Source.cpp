#pragma once

#include "Matrix_CSR.h"

int main()
{
	MatrixCSR mtxCSR;
	mtxCSR.ReadMtx();
	mtxCSR.ConvertMatrixMtxToCCS();
	mtxCSR.CCStoCRS();
	mtxCSR.Write_ccs();
	mtxCSR.Write_crs();
	mtxCSR.Clear();
}

