#pragma once
#include "Matrix_CSR.h"

int main()
{
	MatrixCSR mtxCSR;
	mtxCSR.ReadSortMtx();
	mtxCSR.ConvertMatrixMtxToCSR();
    mtxCSR.CSRtoCSR_t();
    mtxCSR.CalculateParameters();
	mtxCSR.Write_csr();
    mtxCSR.Write_csr_t();
    mtxCSR.path_to_matrix = "out";
    mtxCSR.path_to_picture1 = "685_bus.png";

    mtxCSR.Create_out_html();
	mtxCSR.Clear();
	return 1;
}

