#include <iostream> 
// #include <vector>
// #include <math.h>
#include "MATRIX_LIB.h"

// Stream from file
#include <fstream>

//
#include <Windows.h> // For timing


using std::cin;
using std::cout;
using std::vector;

int test_arraySize(float* array_in){
	return ( sizeof(array_in)/sizeof(float) );
}

void printMatrix(const vector<vector<float>> &M){
	// cout << "M:\n";
	size_t m = M.size();
	size_t n = M[0].size();
	//
	cout << "\n";
	for (size_t i = 0; i < m; ++i){
		for (size_t j = 0; j < n; ++j){
			cout << M[i][j] << "\t";
		}
		cout << "\n";
	}
	cout << "\n";
}
void printMatrix(Mat M){
	// cout << "M:\n";
	size_t m = M.nRow();
	size_t n = M.nCol();
	//
	cout << "\n";
	for (size_t i = 0; i < m; ++i){
		for (size_t j = 0; j < n; ++j){
			cout << M.data[i][j] << "\t";
		}
		cout << "\n";
	}
	cout << "\n";
}

// Test Fast Matrix Inversion
//
// Important:
// The problem is that the blocks may be non-invertible for an invertible matrix.
//
Mat MatrixInversion_fast(Mat M){
	size_t N = M.nRow(); // Suppose that M.nCol() is the same as M.nRow()
	size_t N_2 = N/2;
	size_t N_rest = N - N_2;
	
	Mat M_inv;
	
	// Termination
	if (N <= 1){ // Actually it wont be less than 1
		M_inv.resize(1,1, 1.0/(M.data[0][0]));
		// M_inv = 1.0/(M.data[0][0]);
		return M_inv;
	}
	
	// Blocks
	// A = M.getPart( 0,(N_2-1), 0,(N_2-1) );
	// B = M.getPart( 0,(N_2-1), N_2,(N-1) );
	// C = M.getPart( N_2,(N-1), 0,(N_2-1) );
	// D = M.getPart( N_2,(N-1), N_2,(N-1) );
	
	// Blocks
	Mat A_inv;
	Mat MB;
	Mat C;
	Mat MC;
	Mat E_inv;
	Mat E_inv_MC;
	
	// Other case: N > 1
	// Step 1: Calculate A^-1
	A_inv = MatrixInversion_fast( M.getPart(0,(N_2-1), 0,(N_2-1)) );
	
	// Step 2: Calculate MB and C
	MB = A_inv.dot( M.getPart( 0,(N_2-1), N_2,(N-1) ) );
	C = M.getPart( N_2,(N-1), 0,(N_2-1) );
	
	// Step 3: Calculate E_inv and MC
	E_inv = MatrixInversion_fast( (M.getPart( N_2,(N-1), N_2,(N-1) ).minus(C.dot(MB)) ) );
	MC = (C.dot( A_inv ));
	
	// Step 4: Calculate M_inv
	E_inv_MC = E_inv.dot(MC); 
	M_inv.setPart( ( A_inv.plus(MB.dot(E_inv_MC)) ), 0, 0);
	M_inv.setPart( E_inv_MC.times(-1.0), N_2, 0);
	M_inv.setPart( MB.dot(E_inv).times(-1.0), 0, N_2);
	M_inv.setPart( E_inv, N_2, N_2);
	
	// Return the result
	return M_inv;
}


int main(){
	
	size_t n;
	cin >> n; // Matrix dimension
	
	// Matrices
	Mat M(n, n, 0.0); // The matrix
	Mat M_inv(n, n, 0.0); // The inversion of M
	Mat M_times_M_inv(n, n, 0.0); // The resulting matrix of ( M * M_inv )
	
	Mat M_t(n,n,0.0); // M.', the transpose of M
	Mat M_sum(n, n, 0.0); // M_sum = M + M_times_M_inv
	
		
	// Assign elements
	for (size_t i = 0; i < n; ++i){
		for (size_t j = 0; j < n; ++j){
			cin >> M.data[i][j];
		}
	}
	
	// Timing
	unsigned long dwStart = 0;
	unsigned long dwEnd = 0;
	unsigned long dwDuration;
	//
	float time_1;
	
	// The repeating time
	size_t num_repeat = 300; // 1000;
	
	// Timing
	//////////////////////////////////////////////
	dwStart = GetTickCount();
	//---------------TimingStart-----------------//
	
	
	// Calculate the matrix inversion
	M_inv = M.inverse();
	// M_inv = MatrixInversion_fast(M);
	/*
	for (size_t i = 0; i < num_repeat-1; ++i){
		// Calculate the matrix inversion
		M_inv = M.inverse();
		// M_inv = MatrixInversion_fast(M);
	}
	*/
	
	
	
	
	
	
	
	// Check if M*M_inv is I
	M_times_M_inv = M.dot(M_inv);
	/*
	for (size_t i = 0; i < num_repeat-1; ++i){
		// Check if M*M_inv is I
		M_times_M_inv = M.dot(M_inv);
	}
	*/
	
	
	
	
	
	
	
	
	// Test for transposing
	// M_t = M.T();	
	/*
	for (size_t i = 0; i < num_repeat-1; ++i){
		// Transpose
		M_t = M.T();	
	}
	*/
	
	
	
	
	// Test for summation
	// plus
	// M_sum = M.plus(M_times_M_inv);
	
	/*
	M_sum = M;
	//
	for (size_t i = 0; i < num_repeat; ++i){
		// plus
		// M_sum = M.plus(M_times_M_inv);
		// minus
		// M_sum = M_sum.minus(M_times_M_inv);
		
		
		// Increament	
		// M_sum.increase(M_times_M_inv);
		// M_sum.decrease(M_times_M_inv);
		
		
		
		// Scale up
		// M_sum.scaleUp(4.0);
		// M_sum.scaleUp(0.25);
		
		
		// times
		M_sum = M_sum.times(4.0);
		M_sum = M_sum.times(0.25);
		M_sum = M.times(M_times_M_inv);
		 
	}
	// M_sum.scaleUp(4.0);
	// M_sum = M.times(4.0);
	*/
	
	
	//---------------TimingEnd-----------------//
	dwEnd = GetTickCount();
	dwDuration = dwEnd-dwStart;
	//
	time_1 = double(dwDuration)/double(num_repeat);
	cout << "Average time for single execution: " << time_1 << " ms\n";
	////////////////////////////////////////////// End Timing
	
	
	// Printing results
	bool enable_printingResult = (n <= 10); // true;
	//
	if (enable_printingResult){
		// Print the results
		cout << "\nResult:\n\n";
		// M
		cout << "M:\n";
		printMatrix(M.data);
		// M_inv
		cout << "M_inv:\n";
		printMatrix(M_inv.data);
		// M_inv
		cout << "M_times_M_inv:\n";
		printMatrix(M_times_M_inv.data);
		
		/////////////
		// M_t
		cout << "M_t:\n";
		printMatrix(M_t.data);
		// M_sum
		cout << "M_sum:\n";
		printMatrix(M_sum.data);
	}
	
	////////////////////////
	Mat A;
	Mat B;
	Mat I;
	
	float A_in[] = {1, 2, 0, 4, 5, 6, 0, 8, 9};
	float B_in[] = {1, 0, -1};
	float I_in[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
	A.assign(A_in, 3,3);
	B.assign(B_in, 3, 1);
	I.assign(I_in, 3,3);
	
	cout << "A:\n";
	printMatrix(A);
	cout << "B:\n";
	printMatrix(B);
	cout << "I:\n";
	printMatrix(I);	
	
	
	
	// Test the assignment using stringstream
	//------------------------------------------------//
	Mat AA;

	std::string str_AA;
	
	// Include the parameters as a c-string during complilation time (for embeded system, it's "not suitable" to open a file.)
	/*
	char c_str[] = {
							#include "M.txt"	
						};
	*/
	
	/*
	// Stream from file
	std::ifstream F_handle;
	F_handle.open("M.txt");
	//
	int fileSize = F_handle.gcount();
	cout << "fileSize = " << fileSize << "\n";
	
	F_handle >> str_AA;
	//
	F_handle.close();
	//
	*/
	
	str_AA = "2\n 1 2\n 3 4";
	//
	cout << "str_AA: " << str_AA << "\n";
	
	
	// Creat the stringstream
	std::stringstream ss_AA;
	ss_AA << str_AA; // Assign the string into stream
	
	// Get the first element as the size of the square matrix M
	size_t AA_nRow;
	ss_AA >> AA_nRow;
	
	//
	cout << "AA_nRow = " << AA_nRow << "\n";
	
	// Use the member function "assign(ss, nRow, nCol)"
	AA.assign(ss_AA, AA_nRow, AA_nRow);
	
	// Print the result
	cout << "AA from ss_AA:\n";
	cout << AA.print();
	//
	//------------------------------------------------//
	// end Test the assignment using stringstream
	
	
	//
	Mat C;

	C = A.dot(I);
	cout << "C:\n";
	printMatrix(C);
	
	C = A;
	C.setPart(A.getPart(0,1,0,1), 1,1);	
	cout << "C.setPart(A.getPart(0,1,0,1), 1,1):\n";
	printMatrix(C);
	
	C = A;
	C.setPart(A, 1,2);	
	cout << "C.setPart(A, 1,2):\n";
	printMatrix(C);
	
	C = A.getPart(0,1,1,2);	
	cout << "C = A.getPart(0,1,1,2):\n";
	printMatrix(C);
	
	// Out of range
	C = A.getPart(1,4,1,4);	
	cout << "C = A.getPart(1,4,1,4):\n";
	printMatrix(C);
	//
	
	C = A.getCol(1);
	cout << "C = A.getCol(1):\n";
	printMatrix(C);
	
	C = A.getRow(1);
	cout << "C = A.getRow(1):\n";
	printMatrix(C);
	
	//
	C = A.cat(A, true); // [A, A]
	cout << "C = A.cat(A, true):\n";
	printMatrix(C);
	//
	C = A.cat(A, false); // [A; A]
	cout << "C = A.cat(A, false):\n";
	printMatrix(C);
	
	// Reshape
	C = A;
	C.reshape(1,9,true);
	cout << "C.reshape(1,9,true):\n";
	printMatrix(C);	
	
	C = A;
	C.reshape(2,5,true);
	cout << "C.reshape(2,5,true):\n";
	printMatrix(C);		
	
	C = A;
	C.reshape(2,5,false);
	cout << "C.reshape(2,5,false):\n";
	printMatrix(C);		
	
	C = A;
	C.reshape(2,4,false);
	cout << "C.reshape(2,4,false):\n";
	printMatrix(C);		
	
	C = A;
	
	// Integer power
	C = A.intPower(0);
	cout << "C:\n";
	printMatrix(C);	
	
	C = A.intPower(1);
	cout << "C:\n";
	printMatrix(C);	
	
	C = A.intPower(2);
	cout << "C:\n";
	printMatrix(C);	

	C = A.intPower(-1);
	cout << "C:\n";
	printMatrix(C);	

	C = A.intPower(-2);
	cout << "C:\n";
	printMatrix(C);	
	//
	
	/*
	// The function "Mat_multiply_Mat()" is invisible in the global scope
	// C.data = Mat_multiply_Mat( A.data, B.data ); // This is not going to work since Mat_multiply_Mat() is invisible here
	C.data = MP::Mat_multiply_Mat( A.data, B.data ); // This is able to work
	C.syncSizeInData();
	cout << "C:\n";
	printMatrix(C);
	*/
	
	//
	cout << "A.nRow() = " << A.nRow() << "\n";
	cout << "A.nCol() = " << A.nCol() << "\n";
	cout << "B.nRow() = " << B.nRow() << "\n";
	cout << "B.nCol() = " << B.nCol() << "\n";
	
	cout << "A*B:\n";
	printMatrix(A.dot(B));
	
	cout << "A*A*B:\n";
	printMatrix( A.dot(A.dot(B)) ) ;	

	cout << "B.'*B:\n";
	printMatrix( (B.T()).dot(B) ) ;		
	
	cout << "B*B.':\n";
	printMatrix( B.dot(false, B, true) ) ;		
	
	
	// Test for operators
	////////////////////////////
	Mat D;
	
	C = A;
	
	
	//
	cout << "=============================\n";
	//
	cout << "A:\n";
	printMatrix(A);
	//
	cout << "C:\n";
	printMatrix(C);
	//
	
	//
	cout << "-----------------------------\n";
	
	D = A;
	//
	D = 3.0;
	cout << "D = 3.0:\n";
	printMatrix(D);
	
	vector<float> vec_in(5, 100);
	D = vec_in;
	cout << "D = vec_in:\n";
	printMatrix(D);
	
	vector<vector<float> > MatP_in(2, vector<float>(3,10));
	D = MatP_in;
	cout << "D = MatP_in:\n";
	printMatrix(D);
	
	//
	cout << "-----------------------------\n";
	
	float scale_1;
	
	D = A;
	D.resize(1,1);
	D *= 10;
	cout << "D:\n";
	cout << D.print();
	
	scale_1 = float(D);
	cout << "scale_1 = float(D): scale_1 = " << scale_1 << "\n\n";
	
	D = B.cat(B.cat(B,false),false);
	cout << "D = B.cat(B.cat(B,false),false):\n";
	cout << D.print();
	
	vector<float> vec_1;	
	
	// No supporting for implicit conversion!!
	// vec_1 = D;
	vec_1 = vector<float>(D); // D = B.cat(B.cat(B,false),false);
	cout << "vec_1 = vector<float>(D):\n\n"; 
	for (size_t i = 0; i < vec_1.size(); ++i){
		cout << vec_1[i] << "\n";
	}
	cout << "\n\n";
	
	D = A;
	
	//
	cout << "-----------------------------\n";
	
	cout << "\n";
	
	cout << "A[0][1] = " << A[0][1] << "\n\n";
	
	C = A;
	
	cout << "C[0][1] = 100:\nC:\n";
	C[0][1] = 100;
	cout << C.print();
	
	C = A;
	
	cout << "\n";

	//
	cout << "-----------------------------\n";
	
	cout << "\n";
	cout << "A == B : " << (A == B) << "\n";
	C = A;
	cout << "A == C : " << (A == C) << "\n";
	cout << "A == (2.0 * A) : " << ( A == (2.0 * A) ) << "\n";
	cout << "\n";
	
	//
	cout << "-----------------------------\n";
	
	cout << "\n";
	
	cout << "A^0:\n";
	printMatrix(A^0);

	cout << "A^1:\n";
	printMatrix(A^1);

	cout << "A^2:\n";
	printMatrix(A^2);
	
	cout << "A^-1:\n";
	printMatrix(A^-1);

	cout << "A^-2:\n";
	printMatrix(A^-2);
	
	cout << "\n";

	//
	cout << "----------Non-member-operators:--------------\n";
	

	C = A;
	
	//
	cout << "-----------------------------\n";
	D = -A;
	cout << "D:\n";
	printMatrix(D);
	
		
	//	
	cout << "-----------------------------\n";
	
	D = 1.0 + A;
	cout << "D:\n";
	printMatrix(D);
	
	D = A + 1.0;
	cout << "D:\n";
	printMatrix(D);
	
	D = A + C;
	cout << "D:\n";
	printMatrix(D);
	
	//
	cout << "-----------------------------\n";
	
	D = 1.0 - A;
	cout << "D:\n";
	printMatrix(D);
	
	D = A - 1.0;
	cout << "D:\n";
	printMatrix(D);		
	
	D = A - C;
	cout << "D:\n";
	printMatrix(D);		
	
	//
	cout << "-----------------------------\n";
	
	D = 2.0 * A;
	cout << "D:\n";
	printMatrix(D);
	
	D = A * 2.0;
	cout << "D:\n";
	printMatrix(D);		
	
	D = A * C;
	cout << "D:\n";
	printMatrix(D);		
	
	//
	cout << "-----------------------------\n";
	
	D = 1.0 / A;
	cout << "D:\n";
	printMatrix(D);
	
	D = A / 2.0;
	cout << "D:\n";
	printMatrix(D);		
	
	D = A / C;
	cout << "D:\n";
	printMatrix(D);		
	
	//
	cout << "-----------------------------\n";
	
	D = A;
	D += 1.0;
	cout << "D:\n";
	printMatrix(D);
	
	D = A;
	D += A;
	cout << "D:\n";
	printMatrix(D);		
	
	//
	cout << "-----------------------------\n";
	
	D = A;
	D -= 1.0;
	cout << "D:\n";
	printMatrix(D);		
	
	D = A;
	D -= A;
	cout << "D:\n";
	printMatrix(D);		
	
	//
	cout << "-----------------------------\n";
	
	D = A;
	D *= 2.0;
	cout << "D:\n";
	printMatrix(D);		
	
	D = A;
	D *= 1.0/A;
	cout << "D:\n";
	printMatrix(D);		
	
	//
	cout << "-----------------------------\n";
	
	cout << "B.'*A*B:\n";
	printMatrix((B.T())*A*B);	
	
	//
	cout << "-----------------------------\n";
	
	Mat Z;
	Mat O;
	Mat I_46;
	
	Z.zeros(3,5);
	O.ones(5,7);
	I_46.eye(4,6);
	
	cout << "Z:\n";
	printMatrix(Z);		
	cout << "O:\n";
	printMatrix(O);
	cout << "I_46:\n";
	printMatrix(I_46);

	cout << "B:\n";
	printMatrix(B);	
	cout << "D:\n";
	printMatrix(D);
		
	D.diag(B.cat(B,true));
	cout << "D.diag(B.cat(B,true)):\n";
	printMatrix(D);
	
	D.diag(B.cat(B,false));
	cout << "D.diag(B.cat(B,false)):\n";
	printMatrix(D);
	
	D = A;
	
	//
	cout << "-----------------------------\n";
	
	D = A;
	cout << "D.print():\n";
	cout << D.print();
	
	cout << "D.inverse().print():\n";
	cout << D.inverse().print();
	
	//
	cout << "-----------------------------\n";
	//d
	float array_test[] = {1, 2, 3, 4, 5};
	cout << "( sizeof(array_in)/sizeof(float) ) = " << ( sizeof(array_test)/sizeof(float) ) << "\n";
	cout << "\ntest_arraySize() = " << test_arraySize(array_test) << "\n";
	
	
		
	return 1;
}


