#include "MATRIX_LIB.h"

// No pre-assignments, empty matrix
Mat::Mat(void):
        _nRow(0), _nCol(0),
        data(0, vector<float>())
{

}
// Pre-assign sizes
Mat::Mat(size_t nRow_in, size_t nCol_in):
        _nRow(nRow_in), _nCol(nCol_in),
        data(nRow_in, vector<float>(nCol_in))
{

}
// Pre-assign sizes and elements
Mat::Mat(size_t nRow_in, size_t nCol_in, float element_in):
        _nRow(nRow_in), _nCol(nCol_in),
        data(nRow_in, vector<float>(nCol_in, element_in))
{

}

// Public methods
//----------------------------------//
// Size ---
// "Get" the size of the matrix
size_t Mat::nRow() {
    return this->_nRow;
}
size_t Mat::nCol() {
    return this->_nCol;
}
size_t Mat::nElement() { // Number of elements
    return (_nRow*_nCol);
}
// const functions
size_t Mat::nRow() const {
    return this->_nRow;
}
size_t Mat::nCol() const {
    return this->_nCol;
}
size_t Mat::nElement() const { // Number of elements
    return (_nRow*_nCol);
}
// "Set" and "get" the new size of the matrix
size_t Mat::nRow(size_t nRow_in){
    _nRow = nRow_in;
    // Resize the matrix
    data.resize(_nRow, vector<float>(_nCol));
    return _nRow;
}
size_t Mat::nCol(size_t nCol_in){
    _nCol = nCol_in;
    // Resize the matrix
    for (size_t i = 0; i < _nRow; ++i){
        data[i].resize(_nCol);
    }
    return _nCol;
}
// "Set" and "get" the new size of the matrix with assignment to new elements
size_t Mat::nRow(size_t nRow_in, float element_in){
    _nRow = nRow_in;
    // Resize the matrix
    data.resize(_nRow, vector<float>(_nCol, element_in));
    return _nRow;
}
size_t Mat::nCol(size_t nCol_in, float element_in){
    _nCol = nCol_in;
    // Resize the matrix
    for (size_t i = 0; i < _nRow; ++i){
        data[i].resize(_nCol, element_in);
    }
    return _nCol;
}
// Resize on both row and column
void Mat::resize(size_t nRow_in, size_t nCol_in){
    // Check the size
    if (_nRow == nRow_in && _nCol == nCol_in){
        return; // Nothing to do
    }

    // Change the nCol first to resize the old row vector
    _nCol = nCol_in;
    // Resize the matrix
    for (size_t i = 0; i < _nRow; ++i){
        data[i].resize(_nCol);
    }
    // Then extend the row
    _nRow = nRow_in;
    data.resize(_nRow, vector<float>(_nCol));
}
void Mat::resize(size_t nRow_in, size_t nCol_in, float element_in){ // Assign the new element_in to new elements
    // Check the size
    if (_nRow == nRow_in && _nCol == nCol_in){
        return; // Nothing to do
    }

    // Change the nCol first to resize the old row vector
    _nCol = nCol_in;
    // Resize the matrix
    for (size_t i = 0; i < _nRow; ++i){
        data[i].resize(_nCol, element_in);
    }
    // Then extend the row
    _nRow = nRow_in;
    data.resize(_nRow, vector<float>(_nCol, element_in));
}
// Reshape the matrix, keep all the elements but change the shpae of the matrix (eg. 1 by 6 -> 2 by 3 or 3 by 2)
void Mat::reshape(size_t nRow_in, size_t nCol_in, bool is_Row_first){ // If is_Row_first, keep the elements' order the same in the row direction; otherwise, keep the order the same in column direction.
    // Note: nRow_in * nCol_in should be better be the same as this->nElement(), or some elemnets will be lost or be added (add zeros, actually)
    //
    vector<vector<float> > data_new(nRow_in, vector<float>(nCol_in));
    // Mat M_new(nRow_in, nCol_in);

    // The index for the old matrix
    size_t io = 0;
    size_t jo = 0;
    //
    if (is_Row_first){ // Row go first
        //
        for (size_t i = 0; i < nRow_in; ++i){
            for (size_t j = 0; j < nCol_in; ++j){ // Row go first
                //-------------------//
                // In-range check:
                if (io < this->_nRow && jo < this->_nCol){ // In range
                    data_new[i][j] = this->data[io][jo];
                    // Row go first
                    jo++;
                    // Change row
                    if (jo >= this->_nCol){
                        jo = 0;
                        io++;
                    }
                    //
                }else{ // Out of range
                    data_new[i][j] = 0.0;
                }
                //-------------------//
            }
        }
    }else{ // Column go first
        //
        for (size_t j = 0; j < nCol_in; ++j){
            for (size_t i = 0; i < nRow_in; ++i){ // Column go first
                //-------------------//
                // In-range check:
                if (io < this->_nRow && jo < this->_nCol){ // In range
                    data_new[i][j] = this->data[io][jo];
                    // Column go fist
                    io++;
                    // Change row
                    if (io >= this->_nRow){
                        io = 0;
                        jo++;
                    }
                    //
                }else{ // Out of range
                    data_new[i][j] = 0.0;
                }
                //-------------------//
            }
        }
    }
    //
    // Update this matrix
    // *this = M_new;
    /*
    vector<vector<float>> *temp_ptr = &(this->data);
    //
    &(this->data) = &data_new;
    this->syncSizeInData();
    //
    delete temp_ptr;
    */

    this->data = data_new;
    this->syncSizeInData();
    //
    return;
}
// Synchronize the _nRow and _nCol with the true size of data
void Mat::syncSizeInData(){
    //
    this->_nRow = this->data.size();
    //
    if (this->_nRow == 0){ // Empty
        this->_nCol = 0;
    }else{
        this->_nCol = this->data[0].size();
    }
}

// Assignment ---
void Mat::assign(std::stringstream &ss_in, size_t nRow_in, size_t nCol_in){ // The stringstream ss_in contains a nRow_in by nCol_in matrix.
    // Resize
    this->resize(nRow_in, nCol_in);
    //
    for (size_t i = 0; i < this->_nRow; ++i){
        for (size_t j = 0; j < this->_nCol; ++j){
            ss_in >> this->data[i][j];
        }
    }
}
// From primitive array in c++
void Mat::assign(float* Matrix_in, size_t nRow_in, size_t nCol_in){ // From primitive 1-D array in c++, Matrix_in is a nRow_in by nCol_in matrix.
    // Resize
    this->resize(nRow_in, nCol_in);
    // 
    // Note: the pointer Matrix_in is a copied one, operation on it will not make any effect on the pointer outside this scope.
    //
    for (size_t i = 0; i < this->_nRow; ++i){
        for (size_t j = 0; j < this->_nCol; ++j){
            // this->data[i][j] = Matrix_in[i][j];
            this->data[i][j] = *Matrix_in;
            Matrix_in++;
        }
    }
}
void Mat::assign(const vector<float> &vec_in, bool is_RowVec){ // From 1-D vector
    size_t n = vec_in.size();
    //
    if (is_RowVec){ // Row vector
        //
        this->resize(1,n);
        //
        for (size_t j = 0; j < _nCol; ++j){
            this->data[0][j] = vec_in[j];
        }

        //
    }else{ // Column vector
        //
        this->resize(n,1);
        //
        for (size_t i = 0; i < _nRow; ++i){
            this->data[i][0] = vec_in[i];
        }
        //
    }
}
void Mat::assign(const vector<vector<float> > &MatData_in){ // A = MatData_in, assign a primitive matrix directly
    this->data = MatData_in;
    this->syncSizeInData();
}
// Partial asignment
void Mat::setPart(const Mat &M_part, size_t m_from, size_t n_from){ // The block starts from (m_from, n_from)
    size_t m_to, n_to;
    // The index of the last element of the inserted block
    m_to = m_from + M_part.nRow() - 1;
    n_to = n_from + M_part.nCol() - 1;

    //-------------------------------//
    // Note: if the region of the insertion block exceed the region of this matrix,
    //       resize to a biger one and insert zero at the blank sides.

    // The size of the new matrix
    size_t M,N;
    //
    M = this->_nRow;
    N = this->_nCol;
    //
    // Check if block exceeds the region
    if (m_to >= M ){
        M = m_to + 1;
    }
    if(n_to >= N){
        N = n_to + 1;
    }
    // Deside if the resize action is needed
    this->resize(M, N, 0.0);

    // Start the insertion
    //
    // Initialization
    size_t ip = 0;
    size_t jp = 0;
    for (size_t i = m_from; i <= m_to; ++i ){
        jp = 0;
        for (size_t j = n_from; j <= n_to; ++j){
            // M_part.data[ip][jp] = data[i][j]; // <- get
            this->data[i][j] = M_part.data[ip][jp]; // <- insert
            //
            jp++;
        }
        ip++;
    }

}
// Spetial assignments
void Mat::zeros(size_t nRow_in, size_t nCol_in){ // zeros(m,n)
    resize(0,0);
    resize(nRow_in,nCol_in, 0.0);
}
void Mat::ones(size_t nRow_in, size_t nCol_in){ // ones(m,n)
    resize(0,0);
    resize(nRow_in,nCol_in, 1.0);
}
void Mat::eye(size_t nRow_in, size_t nCol_in){ // Identity matrix, eye(m,n)
    zeros(nRow_in, nCol_in);
    //
    size_t K;
    //
    // min()
    if (nRow_in < nCol_in){
        K = nRow_in;
    }else{
        K = nCol_in;
    }
    //
    for (size_t i = 0; i < K; ++i){
        data[i][i] = 1.0;
    }
    //
}
void Mat::diag(const Mat &M_diag){ // Transform the row/column vector to a diagonal matrix and assign into this matrix
    size_t M = M_diag.nRow();
    size_t N = M_diag.nCol();
    size_t nE = M_diag.nElement();
    //
    this->zeros(nE,nE); // nE by nE matrix
    //
    size_t id = 0;
    size_t jd = 0;
    for (size_t i = 0; i < nE; ++i){
        this->data[i][i] = M_diag.data[id][jd];
        jd++;
        if (jd >= N){
            jd = 0;
            id++;
        }
    }
}


// Get elements ---
float Mat::at(size_t m, size_t n){ // Get the element at (m,n)
    return data[m][n];
}
Mat Mat::getPart(size_t m_from, size_t m_to, size_t n_from, size_t n_to){ // Get partial, M_part = M(m_from:m_to, n_from:n_to)
    int M = m_to - m_from + 1;
    int N = n_to - n_from + 1;
    //
    if (M < 0 || N < 0){
        Mat Empty;
        return Empty;
    }
    //
    Mat M_part(M,N);

    //
    // Initialization
    size_t ip = 0;
    size_t jp = 0;
    //
    for (size_t i = m_from; i <= m_to; ++i ){
        jp = 0;
        for (size_t j = n_from; j <= n_to; ++j){
            // In-range check
            if (i < this->_nRow && j < this->_nCol){
                M_part.data[ip][jp] = this->data[i][j];
            }else{ // Out of bound, set elements to zeros
                M_part.data[ip][jp] = 0.0;
            }
            //
            jp++;
        }
        ip++;
    }
    return M_part;
}
Mat Mat::getCol(size_t n){ // Get a specific column vector in 2-D matrix form
    return getPart(0, (this->_nRow-1), n, n);
}
Mat Mat::getRow(size_t m){ // Get a specific row vector in 2-D matrix form
    return getPart(m, m, 0, (this->_nCol-1));
}
vector<float> Mat::getColVec(size_t n){ // Return a c++ vector, M(:,n)
    vector<float> ColVec(this->_nRow);
    for (size_t i = 0; i < this->_nRow; ++i){
        ColVec[i] = this->data[i][n];
    }
    return ColVec;
}
vector<float> Mat::getRowVec(size_t m){ // Return a c++ vector, M(m,:)
    return this->data[m];
}

// Print out matrix ---
std::string Mat::print(void){ // Print this matrix out as string
    /*
    std::string str_out;
    //
    std::string endl("\n"); // The end line character

    //
    str_out += endl;
    for (size_t i = 0; i < this->_nRow; ++i){
        for (size_t j = 0; j < this->_nCol; ++j){
            str_out += std::to_string( this->data[i][j] ) + "\t";
        }
        str_out += endl;
    }
    str_out += endl;
    //
    return str_out;
    */

    // Using stringstream
    std::stringstream ss_out;
    //
    ss_out << std::endl;
    //
    for (size_t i = 0; i < this->_nRow; ++i){
        for (size_t j = 0; j < this->_nCol; ++j){
            ss_out << this->data[i][j] << "\t";
        }
        ss_out << std::endl;
    }
    ss_out << std::endl;
    //
    return ss_out.str();
}


// Operations =====
// Comparison
bool Mat::is_equal(const Mat &M_in){
    //
    if (M_in.nRow() != this->_nRow || M_in.nCol() != this->_nCol){
        return false;
    }

    // They are of the same size.
    for (size_t i = 0; i < this->_nRow; ++i){
        for (size_t j = 0; j < this->_nCol; ++j){
            if ( this->data[i][j] != M_in.data[i][j] ){
                return false;
            }
        }
    }
    return true;
}

// Self-operation (affecting on this matrix) ---
void Mat::scaleUp(float scale){ // M *= scale
    //
    for (size_t i = 0; i < _nRow; ++i){
        MP::Get_VectorScaleUp(data[i], scale);
    }
}
void Mat::scaleUp_Mat(const Mat &M_right){ // M = M.times(M_right), element-wise multiplication
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            data[i][j] *= M_right.data[i][j];
        }
    }
}
void Mat::increase(float scale){ // M += scale, for each element in M
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            data[i][j] += scale;
        }
    }
}
void Mat::decrease(float scale){ // M -= scale, for each element in M
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            data[i][j] -= scale;
        }
    }
}
void Mat::increase(const Mat &M_right){ // M += M_right
    //
    for (size_t i = 0; i < _nRow; ++i){
        MP::Get_VectorIncrement(data[i], M_right.data[i], false); // +=
    }
}
void Mat::decrease(const Mat &M_right){ // M -= M_right
    //
    for (size_t i = 0; i < _nRow; ++i){
        MP::Get_VectorIncrement(data[i], M_right.data[i], true); // -=
    }
}
// Single-operation, output another matrix ---
// Transpose
Mat Mat::T(void){ // Return the transposed version of this matrix
    //
    Mat M_t(_nCol, _nRow);
    /*
    static Mat M_t;
    M_t.resize(_nCol, _nRow);
    */
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            M_t.data[j][i] = data[i][j];
        }
    }
    return M_t;
}
// Inverse
Mat Mat::inverse(void){ // Return the inversion of this matrix
    // The output matrix

    Mat M_inv;

    if (this->_nRow != this->_nCol){
        M_inv.resize(0,0); // Empty matrix
        return M_inv;
    }


    // M_inv.resize( _nCol, _nRow); // Maybe we will impliment the psodo inverse next time

    //
    M_inv.data = MP::MatrixInversion(this->data);

    // Update the _nRow and _nCol
    M_inv.syncSizeInData();

    return M_inv;
}
Mat Mat::intPower(int power){ // M^power, M^0 -> I, M^-1 -> M_inverse
    // The output matrix
    Mat M_out;

    // Spetial case M^1 -> M
    if (power == 1){
        M_out = *this;
        return M_out;
    }
    // Spetial case M^0 -> I
    if (power == 0){
        M_out.eye(_nRow,_nCol);
        return M_out;
    }

    // Check for square matrix
    if (this->_nRow != this->_nCol){
        M_out.resize(0,0); // Error, return empty matrix
        return M_out;
    }
    //-----------------------------------------------//

    // Spetial case M^-1 -> this->inverse()
    if (power == -1){
        // M^-1, special case for speed up the frequent useage
        return this->inverse();
    }
    //

    // Other cases
    //
    if (power > 1){ // power > 1
        M_out = *this;
        for (size_t i = 0; i < (power-1); ++i){
            M_out = this->dot(M_out);
        }
    }else if (power < -1){ // power < -1
        Mat M_inv = this->inverse();
        M_out = M_inv; // this->inverse()
        for (size_t i = 0; i < (-power - 1); ++i){
            M_out = M_inv.dot(M_out);
        }
    }else{
        // Nothing, power \in {-1, 0, 1} has been processed
    }
    //
    return M_out;

}

// Duo-operation, output another matrix ---
// Plus
Mat Mat::plus(float scale){ // (M + scale), for each element in M
    // Create the output matrix
    Mat M_out(_nRow,_nCol);
    /*
    static Mat M_out;
    M_out.resize(_nRow,_nCol);
    */
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            M_out.data[i][j] = this->data[i][j] + scale;
        }
    }
    return M_out;
}
Mat Mat::minus(float scale){ // (M - scale), for each element in M
    // Create the output matrix
    Mat M_out(_nRow,_nCol);
    /*
    static Mat M_out;
    M_out.resize(_nRow,_nCol);
    */
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            M_out.data[i][j] = this->data[i][j] - scale;
        }
    }
    return M_out;
}
Mat Mat::minus(float scale, bool is_reversed){ // is_reversed -> (scale - M), for each element in M
    // Create the output matrix
    Mat M_out(_nRow,_nCol);
    /*
    static Mat M_out;
    M_out.resize(_nRow,_nCol);
    */
    //
    if (is_reversed){ // Reversed
        for (size_t i = 0; i < _nRow; ++i){
            for (size_t j = 0; j < _nCol; ++j){
                M_out.data[i][j] = scale - this->data[i][j]; // (scale - M)
            }
        }
    }else{
        for (size_t i = 0; i < _nRow; ++i){
            for (size_t j = 0; j < _nCol; ++j){
                M_out.data[i][j] = this->data[i][j] - scale; // (M - scale)
            }
        }
    }
    return M_out;
}
Mat Mat::plus(const Mat &M_right){
    // Create the output matrix
    Mat M_out(_nRow,_nCol);
    /*
    static Mat M_out;
    M_out.resize(_nRow,_nCol);
    */
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            M_out.data[i][j] = this->data[i][j] + M_right.data[i][j];
        }
    }
    return M_out;
}
Mat Mat::minus(const Mat &M_right){
    // Create the output matrix
    Mat M_out(_nRow,_nCol);
    /*
    static Mat M_out;
    M_out.resize(_nRow,_nCol);
    */
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            M_out.data[i][j] = this->data[i][j] - M_right.data[i][j];
        }
    }
    return M_out;
}
// Concatenation
Mat Mat::cat_below(const Mat &M_in){ // Below this matrix, [A; B]
    Mat M_cat;
    // Using setPart(), the size will be automatically changed.
    //------------------------------//
    // Put this matrix at (0, 0)
    M_cat.setPart(*this, 0, 0);
    // Put M_in at (_nRow, 0)
    M_cat.setPart(M_in, this->_nRow, 0);
    //------------------------------//
    return M_cat;
}
Mat Mat::cat_right(const Mat &M_in){ // Right-side of this matrix, [A, B]
    Mat M_cat;
    // Using setPart(), the size will be automatically changed.
    //------------------------------//
    // Put this matrix at (0, 0)
    M_cat.setPart(*this, 0, 0);
    // Put M_in at (0, _nCol)
    M_cat.setPart(M_in, 0, this->_nCol);
    //------------------------------//
    return M_cat;
}
Mat Mat::cat(const Mat &M_in, bool is_horizontal){ // is_horizontal --> cat_Right(); otherwise --> cat_Below()
    //
    if(is_horizontal){ // [A, B]
        return this->cat_right(M_in);
    }else{ // [A; B]
        return this->cat_below(M_in);
    }
}

// Scalar/Element-wise multiplication
Mat Mat::times(float scale){ // Scalar multiplication
    // Create the output matrix
    Mat M_out(_nRow,_nCol);
    //
    /*
    static Mat M_out;
    M_out.resize(_nRow, _nCol);
    */


    if(scale == -1.0){
        //
        for (size_t i = 0; i < _nRow; ++i){
            for (size_t j = 0; j < _nCol; ++j){
                M_out.data[i][j] = -(this->data[i][j]);
            }
        }
    }else{
        //
        for (size_t i = 0; i < _nRow; ++i){
            for (size_t j = 0; j < _nCol; ++j){
                M_out.data[i][j] = this->data[i][j] * scale;
            }
        }
    }
    return M_out;
}
Mat Mat::times(const Mat &M_right){ // Element-wise multiplication
    // Create the output matrix
    Mat M_out(_nRow,_nCol);
    /*
    static Mat M_out;
    M_out.resize(_nRow,_nCol);
    */
    //
    for (size_t i = 0; i < _nRow; ++i){
        for (size_t j = 0; j < _nCol; ++j){
            M_out.data[i][j] = this->data[i][j] * M_right.data[i][j];
        }
    }
    return M_out;
}


// Matrix multiplication
// Note: this matrix is the "left" one
Mat Mat::dot(const Mat &M_right){ // Similar to the nomenclature of numpy in Python

    // Size check
    // mxn = (mxk)*(kxn)
    size_t M, K, N;
    //
    M = _nRow;
    K = _nCol;
    N = M_right.nCol();

    // The output matrix
    Mat M_out(M, N);
    /*
    static Mat M_out;
    M_out.resize(M, N);
    */

    // cout << "here in\n";
    float* ptr_data = NULL;
    // Using indexing
    for (size_t i = 0; i < M; ++i){ // row in m_out
        //
        // cout << "#i: " << i << "\n";
        //
        for (size_t j = 0; j < N; ++j){ // column in m_out
            //
            // cout << "#j: " << j << "\n";
            //
            ptr_data = &(M_out.data[i][j]);
            //
            *ptr_data = 0.0;
            for (size_t k = 0; k < K; ++k){
                if (data[i][k] != 0.0 && M_right.data[k][j] != 0.0)
                    *ptr_data += data[i][k]*M_right.data[k][j];
            }
        }
    }

    return M_out;
}
Mat Mat::dot(bool Transpose_left, const Mat &M_right, bool Transpose_right){ // Extended version for conbining the ability of transpose of both mtrices

    // Size check
    // mxn = (mxk)*(kxn)
    size_t M,K,N;
    // Left transpose
    if (Transpose_left){
        M = _nCol;
        K = _nRow;
    }else{
        M = _nRow;
        K = _nCol;
    }
    // Right transpose
    if (Transpose_right){
        N = M_right.nRow();
    }else{
        N = M_right.nCol();
    }

    // The output matrix
    Mat M_out(M, N);
    /*
    static Mat M_out;
    M_out.resize(M, N);
    */

    //
    float* ptr_data = NULL;
    //

    // Check the conditions of transpotations
    if(Transpose_left){
        if(Transpose_right){ // Both transposed
            //
            // Using indexing
            for (size_t i = 0; i < M; ++i){ // row in m_out
                for (size_t j = 0; j < N; ++j){ // column in m_out
                    // M_out.data[i][j] = 0.0;
                    ptr_data = &(M_out.data[i][j]);
                    *ptr_data = 0.0;
                    //
                    for (size_t k = 0; k < K; ++k){
                        if (data[k][i] != 0.0 && M_right.data[j][k] != 0.0) // (i,k) -> (k,i), (k,j) -> (j,k)
                            *ptr_data += data[k][i]*M_right.data[j][k]; // (i,k) -> (k,i), (k,j) -> (j,k)
                    }
                }
            }
            //
        }else{ // Only left transpose
            //
            // Using indexing
            for (size_t i = 0; i < M; ++i){ // row in m_out
                for (size_t j = 0; j < N; ++j){ // column in m_out
                    // M_out.data[i][j] = 0.0;
                    ptr_data = &(M_out.data[i][j]);
                    *ptr_data = 0.0;
                    //
                    for (size_t k = 0; k < K; ++k){
                        if (data[k][i] != 0.0 && M_right.data[k][j] != 0.0) // (i,k) -> (k,i)
                            *ptr_data += data[k][i]*M_right.data[k][j]; // (i,k) -> (k,i)
                    }
                }
            }
            //
        }
    }else{
        if(Transpose_right){ // Only right transpose
            //
            // Using indexing
            for (size_t i = 0; i < M; ++i){ // row in m_out
                for (size_t j = 0; j < N; ++j){ // column in m_out
                    // M_out.data[i][j] = 0.0;
                    ptr_data = &(M_out.data[i][j]);
                    *ptr_data = 0.0;
                    //
                    for (size_t k = 0; k < K; ++k){
                        if (data[i][k] != 0.0 && M_right.data[j][k] != 0.0) // (k,j) -> (j,k)
                            *ptr_data += data[i][k]*M_right.data[j][k]; // (k,j) -> (j,k)
                    }
                }
            }
            //
        }else{ // Normal
            //
            // Using indexing
            for (size_t i = 0; i < M; ++i){ // row in m_out
                for (size_t j = 0; j < N; ++j){ // column in m_out
                    // M_out.data[i][j] = 0.0;
                    ptr_data = &(M_out.data[i][j]);
                    *ptr_data = 0.0;
                    //
                    for (size_t k = 0; k < K; ++k){
                        if (data[i][k] != 0.0 && M_right.data[k][j] != 0.0)
                            *ptr_data += data[i][k]*M_right.data[k][j];
                    }
                }
            }
            //
        }
    }

    return M_out;
}
// Operator overloading
//---------------------------------------//
// Member
// A <- b, assignment
Mat& Mat::operator = (float scale){ // Assign a real number as 1 by 1 matrix
    // this->resize(1,1, scale);
    this->resize(1,1);
    this->data[0][0] = scale;
    return *this;
}
Mat& Mat::operator = (const vector<float> &colVec_in){ // A = vec_in, assign as a column vector
    this->assign(colVec_in, false); // A column vector
    return *this;
}
Mat& Mat::operator = (const vector<vector<float> > &MatData_in){ // A = MatData_in, assign a primitive matrix directly
    this->assign(MatData_in);
    return *this;
}
// b <- A, get value, using implicit type conversion
// Note: no implicit conversion to float, since it would be ambiguous in other operator overriding
Mat::operator float (){ // Get the first element as float, good for a 1 by 1 matrix
    return this->data[0][0];
}
Mat::operator std::vector<float> (){ // Get the column vector as std::vector<float> in c++
    return this->getColVec(0);
}
// A[]
vector<float>& Mat::operator [] (size_t m){ // Indexing the m-th row, equivalent to data[m]
    return this->data[m];
}
// A == B
bool Mat::operator == (Mat const& B){ // is_equal()
    return this->is_equal(B);
}
// A^z, z \in Z
Mat Mat::operator ^ (int power){ // A^power, this->intPower()
    return this->intPower(power);
}


// Non-member
//
// -A
Mat operator - (const Mat &A){
    return Mat(A).times(-1.0);
}
// A + B
Mat operator + (float a, const Mat &B){
    return Mat(B).plus(a);
}
Mat operator + (const Mat &A, float b){
    return Mat(A).plus(b);
}
Mat operator + (const Mat &A, const Mat &B){
    return Mat(A).plus(B);
}
// A - B
Mat operator - (float a, const Mat &B){
    return Mat(B).minus(a, true); // (a - B)
}
Mat operator - (const Mat &A, float b){
    return Mat(A).minus(b);
}
Mat operator - (const Mat &A, const Mat &B){
    return Mat(A).minus(B);
}
// A * B
Mat operator * (float a, const Mat &B){
    return Mat(B).times(a);
}
Mat operator * (const Mat &A, float b){
    return Mat(A).times(b);
}
Mat operator * (const Mat &A, const Mat &B){
    // Matrix multiplication
    return Mat(A).dot(B);
}
// A/B, including matrix inversion (eg, (1.0/B) <--> B^-1)
Mat operator / (float a, const Mat &B){ // a*(B^-1), for that B to be inversed
    if (a == 1.0){
        return Mat(B).inverse();
    }
    return ( ( Mat(B).inverse()).times(a) );
}
Mat operator / (const Mat &A, float b){ // A*(1/b), scalar multiplication of 1/b
    return Mat(A).times(1.0/b);
}
Mat operator / (const Mat &A, const Mat &B){ // A*(B^-1), for that B to be inversed
    return ( Mat(A).dot(Mat(B).inverse()) );
}
// A += B
Mat& operator += (Mat &A, float b){
    A.increase(b);
    return A;
}
Mat& operator += (Mat &A, const Mat &B){
    A.increase(B);
    return A;
}
// A -= B
Mat& operator -= (Mat &A, float b){
    A.decrease(b);
    return A;
}
Mat& operator -= (Mat &A, const Mat &B){
    A.decrease(B);
    return A;
}
// A *= B
Mat& operator *= (Mat &A, float b){
    A.scaleUp(b);
    return A;
}
Mat& operator *= (Mat &A, const Mat &B){ // Matrix multiplication
    A = A.dot(B);
    return A;
}

//---------------------------------------//
// end Operator overloading
