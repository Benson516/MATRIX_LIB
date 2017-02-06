#ifndef MATRIX_LIB_H
#define MATRIX_LIB_H
//
#include <vector>
#include <math.h>

// For printing matrix out
////////////////////////
#include <sstream>
#include <string>
//////////////////////// end For printing matrix out


#include "MATRIX_PRIMITIVE.h"

using std::vector;
// using std::string;

//
namespace MP = MATRIX_PRIMITIVE;
//

/*
// Debugging
//////////////////
#include <iostream>
using std::cout;
////////////////// end Debugging
*/

class Mat{
public:
    // Matrix data
    vector<vector<float> > data; // data is a m by n matrix


    // Create a nRow_in by nCol_in matrix
    Mat(void); // No pre-assignments, empty matrix
    Mat(size_t nRow_in, size_t nCol_in); //
    Mat(size_t nRow_in, size_t nCol_in, float element_in); // Pre-assign elements

    // Public methods
    //----------------------------------//
    // Size ---
    // "Get" the size of the matrix
    size_t nRow();
    size_t nCol();
    size_t nElement(); // Number of elements
    // "Set" and "get" the new size of the matrix
    size_t nRow(size_t nRow_in);
    size_t nCol(size_t nCol_in);
    // "Set" and "get" the new size of the matrix with assignment to new elements
    size_t nRow(size_t nRow_in, float element_in);
    size_t nCol(size_t nCol_in, float element_in);
    // Resize on both row and column
    void resize(size_t nRow_in, size_t nCol_in);
    void resize(size_t nRow_in, size_t nCol_in, float element_in); // Assign the new element_in to new elements
    // Reshape the matrix, keep all the elements but change the shpae of the matrix (eg. 1 by 6 -> 2 by 3 or 3 by 2)
    void reshape(size_t nRow_in, size_t n_Col_in, bool is_Row_first); // If is_Row_first, keep the elements' order the same in the row direction; otherwise, keep the order the same in column direction.
    // Synchronize the _nRow and _nCol with the true size of data
    void syncSizeInData();

    // Assignment ---
    void assign(std::stringstream &ss_in, size_t nRow_in, size_t nCol_in); // The stringstream ss_in contains a nRow_in by nCol_in matrix.
    void assign(float* Matrix_in, size_t nRow_in, size_t nCol_in); // From primitive 1-D array in c++, Matrix_in is a nRow_in by nCol_in matrix.
    void assign(vector<float> &vec_in, bool is_RowVec); // From 1-D vector
    void assign(vector<vector<float> > &MatData_in); // A = MatData_in, assign a primitive matrix directly
    // Partial asignment
    void setPart(Mat M_part, size_t m_from, size_t n_from); // The block starts from (m_from, n_from)
    // Spetial assignments
    void zeros(size_t nRow_in, size_t nCol_in); // zeros(m,n)
    void ones(size_t nRow_in, size_t nCol_in); // ones(m,n)
    void eye(size_t nRow_in, size_t nCol_in); // Identity matrix, eye(m,n)
    void diag(Mat M_diag); // Transform the row/column vector to a diagonal matrix and assign into this matrix

    // Get elements ---
    float at(size_t m, size_t n); // Get the element at (m,n)
    Mat getPart(size_t m_from, size_t m_to, size_t n_from, size_t n_to); // Get partial, M_part = M(m_from:m_to, n_from:n_to)
    Mat getCol(size_t n); // Get a specific column vector in 2-D matrix form
    Mat getRow(size_t m); // Get a specific row vector in 2-D matrix form
    vector<float> getColVec(size_t n); // Return a c++ vector, M(:,n)
    vector<float> getRowVec(size_t m); // Return a c++ vector, M(m,:)

    // Print out matrix ---
    std::string print(void); // Print this matrix out as string



    // Operations =====

    // Comparison
    bool is_equal(Mat M_in);

    // Self-operation (affecting on this matrix) ---
    void scaleUp(float scale); // M *= scale
    void scaleUp_Mat(Mat M_right); // M = M.times(M_right), element-wise multiplication
    //
    void increase(float scale); // M += scale, for each element in M
    void decrease(float scale); // M -= scale, for each element in M
    //
    void increase(Mat M_right); // M += M_right
    void decrease(Mat M_right); // M -= M_right


    // Single-operation, output another matrix ---
    Mat T(void); // Transpose, return the transposed version of this matrix
    Mat inverse(void); // Inverse, return the inversion of this matrix
    Mat intPower(int power); // M^power, M^0 -> I, M^-1 -> M_inverse

    // Duo-operation, output another matrix ---
    // Plus
    Mat plus(float scale); // (M + scale), for each element in M
    Mat minus(float scale); // (M - scale), for each element in M
    Mat minus(float scale, bool is_reversed); // is_reversed -> (scale - M), for each element in M
    Mat plus(Mat M_right);
    Mat minus(Mat M_right);
    // Concatenation
    Mat cat_below(Mat M_in); // Below this matrix, [A; B]
    Mat cat_right(Mat M_in); // Right-side of this matrix, [A, B]
    Mat cat(Mat M_in, bool is_horizontal); // is_horizontal --> cat_Right(); otherwise --> cat_Below()


    // Scalar/Element-wise multiplication
    Mat times(float scale); // Scalar multiplication
    Mat times(Mat M_right); // Element-wise multiplication

    // Matrix multiplication
    // Note: this matrix is the "left" one
    Mat dot(Mat M_right); // Similar to the nomenclature of numpy in Python
    Mat dot(bool Transpose_left, Mat M_right, bool Transpose_right); // Extended version for conbining the ability of transpose of both mtrices

    // Operator overloading
    //----------------------------//
    // A <- b, assignment
    Mat& operator = (float scale); // Assign a real number as 1 by 1 matrix
    Mat& operator = (vector<float> &colVec_in); // A = vec_in, assign as a column vector
    Mat& operator = (vector<vector<float> > &MatData_in); // A = MatData_in, assign a primitive matrix directly
    // b <- A, get value, using implicit type conversion
    // Note: no implicit conversion to float, since it would be ambiguous in other operator overriding
    explicit operator float (); // Get the first element as float, good for a 1 by 1 matrix
    explicit operator std::vector<float> (); // Get the column vector as std::vector<float> in c++
    // A[]
    vector<float>& operator [] (size_t m); // Indexing the m-th row, equivalent to data[m]
    // A == B
    bool operator == (Mat const& B); // is_equal()
    // A^z, z \in Z
    Mat operator ^ (int power); // A^power, this->intPower()

    //----------------------------//
    // end Operator overloading
private:

    size_t _nRow; // Number of rows
    size_t _nCol; // Number of columns


    // Private methods
    // void get_maxBound(size_t &M_max, size_t &N_max, size_t M_new, size_t N_new); // Cauculate the approproate size to contain all the matrices

};

// Operator overloading
//---------------------------------------//


// -A
Mat operator - (Mat A);
// A + B
Mat operator + (float a, Mat B);
Mat operator + (Mat A, float b);
Mat operator + (Mat A, Mat B);
// A - B
Mat operator - (float a, Mat B);
Mat operator - (Mat A, float b);
Mat operator - (Mat A, Mat B);
// A * B
Mat operator * (float a, Mat B);
Mat operator * (Mat A, float b);
Mat operator * (Mat A, Mat B); // Matrix multiplication
// A/B, including matrix inversion (eg, (1.0/B) <--> B^-1)
Mat operator / (float a, Mat B); // a*(B^-1), for that B to be inversed
Mat operator / (Mat A, float b); // A*(1/b), scalar multiplication of 1/b
Mat operator / (Mat A, Mat B); // A*(B^-1), for that B to be inversed
// A += B
Mat& operator += (Mat &A, float b);
Mat& operator += (Mat &A, Mat B);
// A -= B
Mat& operator -= (Mat &A, float b);
Mat& operator -= (Mat &A, Mat B);
// A *= B
Mat& operator *= (Mat &A, float b);
Mat& operator *= (Mat &A, Mat B); // Matrix multiplication

//---------------------------------------//
// end Operator overloading

#endif
