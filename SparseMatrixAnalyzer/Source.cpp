#pragma once
#include "Matrix_CSR.h"

int main()
{
	MatrixCSR mtxCSR;
	mtxCSR.ReadSortMtx();
	mtxCSR.ConvertMatrixMtxToCSR();
	mtxCSR.CheckExistenceOfIsolatedSubmatrices();
	mtxCSR.Write_crs();
    mtxCSR.path_to_matrix = "out";
    mtxCSR.path_to_picture1 = "685_bus.png";
    mtxCSR.sz_row = 0;
    mtxCSR.sz_col = 0;
    mtxCSR.is_init_data = 0;
    mtxCSR.is_symmetrical = 0;
    mtxCSR.MinElement = 0;
    mtxCSR.MaxElement = 0;
    mtxCSR.MinModElement = 0;
    mtxCSR.MaxModElement = 0;
    mtxCSR.MaxDiag = 0;
    mtxCSR.MinDiag = 0;
    mtxCSR.MaxModDiag = 0;
    mtxCSR.MinModDiag = 0;
    mtxCSR.Create_out_html();
	mtxCSR.Clear();
	return 1;
}

