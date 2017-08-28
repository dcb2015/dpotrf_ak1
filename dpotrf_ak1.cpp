// dpotrf_ak1.cpp -- This program computes the Cholesky factorization of a real symmetric positive definite matrix.
//
// The main routine is a translation of the FORTRAN routine
// DPOTRF.F, from Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.
//
// This translation makes a copy of the original input array and works with the copy.
// This way, the original matrix is saved.
//
// This translation uses the upper triangular part of the matrix, computing and outputing the upper
// half of the factorization such that A = U^T * U
//
// References:
//
// 
// 
// 
// 
//
// 
// 
// 
// 
//
// 
// 
// 
// 
// 
// 
//
// To distinguish this version of fmin from other translations,
// an '_ak1' suffix has been appended to its name.
//
// 27 August 2017
//
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cctype>
#include <vector>

using namespace std;

typedef vector<vector<double> > C2DArray;

void Echo2DArray(const C2DArray a2);

void Echo2DArray(const C2DArray a2) {
	//Routine to output a 2D Array of type double to the console
	cout << "\n";
	for (unsigned int i = 0; i < a2.size(); i++) {
		for (unsigned int j = 0; j < a2[0].size(); j++) {
			cout << a2[i][j] << " ";
		} // End for j
		cout << "\n";
	} // End for i
	cout << "\n";
	return;
} // End Echo2DArray

void dpotrf_ak1(const int N, C2DArray& A_Matrix, int *info);

void dpotrf_ak1(const int N, C2DArray& A_Matrix, int *info){

	// For this program, it is assumed that UPLO = 1, so it will be hard-coded to use the upper
	// triangular part of the matrix

	// On output:
	// info = 0:		successful exit
	// info = k > 0:	the leading minor of order k is not positive definite,
	//					and the factorization could not be completed.

	// Use unblocked code: call DPOTF2

	int i, j, jj, jy;
	double ajj, ddot;
	bool ffFlag = 0;

	// BEGIN DPOTF2

	for (j = 0; j < N; ++j){

		ajj = A_Matrix[j][j];
		ddot = 0.0;

		if (ffFlag){

			// DDOT
			i = j;
			do {
				--i;
				ddot += A_Matrix[i][j] * A_Matrix[i][j];
			} while (i > 0);
			// End DDOT

			ajj -= ddot;

		} // End if ffFlag

		//else  ffFlag = 1;

		//for (i = 0; i < j; ++i)  ddot += A_Matrix[i][j] * A_Matrix[i][j];  // DDOT

		//cout << "j = " << j << "      " << "ddot = " << ddot << endl;

		if (ajj <= 0){
			A_Matrix[j][j] = ajj;
			*info = j + 1;
			return;
		} // End if ajj <= 0

		A_Matrix[j][j] = ajj = sqrt(ajj);

		// Compute elements J+1:N of row J

		if (j < (N - 1)){

			// BEGIN DGEMV

			//if (j > 0) {

				// Start the operations. In this version, the elements of A are accessed sequentially
				// with one pass through A.

				// Since this program assumes we are working with the upper diagonal matrix of A,
				// form y = alpha*A^T * x + y

				//if (j1Flag){
				if (ffFlag){

					//jy = j + 1;
					//for (jj = 0; jj < (N - 1 - j); ++jj){
					for (jj = j + 1; jj < N; ++jj){
						ddot = 0.0;
						for (i = 0; i < j; ++i){
							//ddot += A_Matrix[i][jy] * A_Matrix[i][j];  // A * X
							ddot += A_Matrix[i][jj] * A_Matrix[i][j];  // A * X
						} // End for i
						//A_Matrix[j][jy] -= ddot;
						A_Matrix[j][jj] -= ddot;
						//++jy;
					} // End for jj

				} // End if ffFlag
				else ffFlag = 1;

			//} // End if j > 0

			// END DGEMV

			// DSCAL
			for (i = (j + 1); i < N; ++i) A_Matrix[j][i] /= ajj;

		} // End if (j < (N - 1))

	} // End for j

	// END DPOTF2

	return;
} // End dpotrf_ak1

int main() {
	char rflag;			//Readiness flag

	cout << "                       dpotrf_ak1 (27 August 2017)" << endl;
	cout << "======================================================================" << endl << endl;
	cout << "This program computes the Cholesky factorization of a real symmetric matrix [A]." << endl << endl;
	cout << "The upper half of the factorization is computed such that A = U^T * U." << endl;
	cout << "The strictly lower triangular part of A is not referenced." << endl << endl;
	cout << "The dimension of this matrix, as well as the entries of the matrix should have" << endl;
	cout << "been saved beforehand in a file named inReSymA.txt" << endl << endl;
	cout << "inReSymA.txt should be in the same folder as the dpotrf executable." << endl << endl;
	cout << "The data is assumed to be of type double. Variables used within this program" << endl;
	cout << "are type double." << endl << endl;
	cout << "The output is written to the file outCholFac.txt." << endl << endl;
	cout << "The results are calculated to double precision--" << DBL_DIG << " decimal places." << endl << endl;
	cout << "Everything ready? If yes, press y." << endl;
	cout << "Otherwise, press any other key." << endl;
	cin >> rflag;

	if (toupper(rflag) == 'Y') {

		C2DArray A, B; // A and B Matrices
		int i, j, info = 0, mDim; // Array index variables, flag, and matrix dimension
		bool erFlag = false;  // Error flag

		ifstream in("inReSymA.txt", ios::in);

		if (!in) {
			cout << "\nCannot open the input file." << endl;
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;
		}

		in >> mDim; //Input the Matrix dimension from the file, N x N

		if (mDim < 1){
			in.close(); //Close the input file before terminating
			cout << "\nInvalid dimension entered. Program terminated." << endl;
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;
		}

		ofstream out("outCholFac.txt", ios::out);

		if (!out) {
			in.close(); //Close the input file before terminating
			cout << "\nCannot open the output file. Program terminated." << endl;
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;
		}

		// Beginning of try block, if vector re-sizing unsuccessful
		try {

			// Resize the arrays to the appropriate sizes
			A.resize(mDim); // N rows
			B.resize(mDim); // N rows

			for (i = 0; i < mDim; ++i){
				A[i].resize(mDim); // N columns
				B[i].resize(mDim); // N columns
			} // End for i

		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions

			in.close();
			out.close();
			cerr << "\nIn catch block, so an exception occurred: " << xa.what() << endl;
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;

		} // End of catch block

		for (i = 0; i < mDim; ++i){ //Input the A Matrix from the file
			for (j = 0; j < mDim; ++j){
				in >> A[i][j];
			}//End for j
		}//End for i

		in.close(); //Close the input file

		Echo2DArray(A);

		// Confirm the symmetry of the matrix
		i = mDim - 1;
		while ((i > 0) && (!erFlag)){

			j = i;
			do {
				--j;
				if (fabs(A[i][j] - A[j][i]) > DBL_EPSILON){
					erFlag = true;
					break;
				} // End if
			} while (j > 0); // End do-while

			--i;

		} // End while (i > 0)

		if (erFlag) {
			cout << "Non-symmetry in matrix detected. No further action taken. Program terminated." << endl << endl;
			out.close();
			cout << "\nEnter any key to continue." << endl;
			cin >> rflag;
			return 0;
		} // End if erFlag

		cout << "Matrix seems to be symmetric." << endl << endl;

		// Copy A matrix into B and use B as the working matrix

		for (i = 0; i < mDim; ++i){
			for (j = i; j < mDim; ++j) B[i][j] = A[i][j];
		}//End for i

		dpotrf_ak1(mDim, B, &info);

		if (info != 0){
			cout << "The leading minor of order " << info << "is not positive definite," << endl;
			cout << "and the factorization could not be completed." << endl;
		}
		else { // else info == 0

			cout << "The factorization was completed successfully." << endl;
			out.precision(DBL_DIG);

			out << "The factored matrix follows:" << endl << endl;

			for (i = 0; i < mDim; ++i){
				for (j = 0; j < mDim; ++j) out << B[i][j] << "   ";
				out << endl;
			} // End for i
		} // End else info == 0

		out.close();

	}	//End if ready
	else cout << "Not ready. Try again when ready with information.\n";

	cout << "\nEnter any key to continue." << endl;
	cin >> rflag;
	return 0;
}	//End main program
