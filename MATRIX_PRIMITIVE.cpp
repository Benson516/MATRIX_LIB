#include "MATRIX_PRIMITIVE.h"

//
// using namespace MATRIX_PRIMITIVE;
//

// Namespace MATRIX_PRIMITIVE
//////////////////////////
namespace MATRIX_PRIMITIVE
{
// Utilities
///////////////////////////
void Mat_multiply_Vec(vector<float> &v_out, vector<vector<float> > &m_left, vector<float> &v_right){ // v_out = m_left*v_right
    vector<float>::iterator it_out;
    vector<float>::iterator it_m_row;
    vector<float>::iterator it_v;
    //
    it_out = v_out.begin();
    for (size_t i = 0; i < m_left.size(); ++i){
        *it_out = 0.0;
        it_m_row = m_left[i].begin();
        it_v = v_right.begin();
        for (size_t j = 0; j < m_left[i].size(); ++j){
            // *it_out += m_left[i][j] * v_right[j];
            if (*it_m_row != 0.0 && *it_v != 0.0){
                (*it_out) += (*it_m_row) * (*it_v);
            }else{
                // (*it_out) += 0.0
            }
            // (*it_out) += (*it_m_row) * (*it_v);
            //
            it_m_row++;
            it_v++;
        }
        it_out++;
    }
}
vector<float> Mat_multiply_Vec(vector<vector<float> > &m_left, vector<float> &v_right){ // v_out = m_left*v_right
    static vector<float> v_out;
    // Size check
    if (v_out.size() != m_left.size()){
        v_out.resize(m_left.size());
    }
    // Iterators
    vector<float>::iterator it_out;
    vector<float>::iterator it_m_row;
    vector<float>::iterator it_v;
    //
    it_out = v_out.begin();
    for (size_t i = 0; i < m_left.size(); ++i){
        *it_out = 0.0;
        it_m_row = m_left[i].begin();
        it_v = v_right.begin();
        for (size_t j = 0; j < m_left[i].size(); ++j){
            // *it_out += m_left[i][j] * v_right[j];
            if (*it_m_row != 0.0 && *it_v != 0.0){
                (*it_out) += (*it_m_row) * (*it_v);
            }else{
                // (*it_out) += 0.0
            }
            // (*it_out) += (*it_m_row) * (*it_v);
            //
            it_m_row++;
            it_v++;
        }
        it_out++;
    }
    return v_out;
}

// New function for matrix multiply matrix
vector<vector<float> > Mat_multiply_Mat(vector<vector<float> > &m_left, vector<vector<float> > &m_right){ // m_out = m_left*m_right
    static vector<vector<float> > m_out;

	// mxn = (mxk)*(kxn)
    size_t M, K, N;
    //
    M = m_left.size();
    K = m_left[0].size();
    N = m_right[0].size();

    // Size check
    if (m_out.size() != M || m_out[0].size() != N ){
        m_out.resize(M, vector<float>( N ) );
    }

    // Using indexing
    for (size_t i = 0; i < M; ++i){ // row in m_out
        for (size_t j = 0; j < N; ++j){ // column in m_out
        	m_out[i][j] = 0.0;
			for (size_t k = 0; k < K; ++k){
				if (m_left[i][k] != 0.0 && m_right[k][j] != 0.0)
					m_out[i][j] += m_left[i][k]*m_right[k][j];
			}
		}
    }

    return m_out;
}

vector<float> Get_VectorPlus(const vector<float> &v_a, const vector<float> &v_b, bool is_minus) // v_a + (or -) v_b
{
    static vector<float> v_c;
    // Size check
    if (v_c.size() != v_a.size()){
        v_c.resize(v_a.size());
    }
    //
    for (size_t i = 0; i < v_a.size(); ++i){
        if (is_minus){
            v_c[i] = v_a[i] - v_b[i];
        }else{
            v_c[i] = v_a[i] + v_b[i];
        }
    }
    return v_c;
}
vector<float> Get_VectorScalarMultiply(const vector<float> &v_a, float scale) // scale*v_a
{
    static vector<float> v_c;
    // Size check
    if (v_c.size() != v_a.size()){
        v_c.resize(v_a.size());
    }
    // for pure negative
    if (scale == -1.0){
        for (size_t i = 0; i < v_a.size(); ++i){
            v_c[i] = -v_a[i];
        }
        return v_c;
    }
    // else
    for (size_t i = 0; i < v_a.size(); ++i){
        v_c[i] = scale*v_a[i];

    }
    return v_c;
}

// Important!
// New function for scale-up a vector
// v_a *= scale
void Get_VectorScaleUp(vector<float> &v_a, float scale) // v_a *= scale
{
    // for pure negative
    if (scale == -1.0){
        for (size_t i = 0; i < v_a.size(); ++i){
            v_a[i] = -v_a[i];
        }
        //
    }else{
        // else
        for (size_t i = 0; i < v_a.size(); ++i){
            v_a[i] *= scale;
        }
        //
    }
}

// Increment
void Get_VectorIncrement(vector<float> &v_a, const vector<float> &v_b, bool is_minus){ // v_a += (or -=) v_b
    // Size check
    if (v_a.size() != v_b.size()){
        v_a.resize(v_b.size());
    }
    //
    if (is_minus){ // -=
        for (size_t i = 0; i < v_b.size(); ++i){
            v_a[i] -= v_b[i];
        }
    }else{ // +=
        for (size_t i = 0; i < v_b.size(); ++i){
            v_a[i] += v_b[i];
        }
    }

}
/////////////////////////// end Utilities


// Matrix inversion, using Gauss method
/////////////////////////////////////////
bool SolveSingularityOnDiag(vector<vector<float> > &M, vector<vector<float> > &M_inv, size_t i){
    // Return value: is_singular
    //-----------------------------------//
    // true: the matrix is singular
    // false: the singularity is not yet been found or the matrix is invertible
    //-----------------------------------//


    // For ith-row
    //-----------------------------------//
	// If the ith element on diagonal is zero, find a row with it ith elemnt non-zero and add to the current row (ith row)
	// For both M and M_inv

	size_t n = M.size();
	// Search for other non-zero elements
	float max_absValue = 0.0;
	size_t idx_max = i;

    // The following block of code is wrong.
    /*
	for (size_t j = 0; j < n; ++j){
		if (j != i){ // Other than i
			float absValue = fabs(M[j][i]);
			if ( absValue > max_absValue){
				max_absValue = absValue;
				idx_max = j;
				//
				// break; // Once found a non-zero element, break it (Not going to find the maximum one)
				//
			}
		}
	}
    */

    // It's sufficient and "should" only looking downward
    for (size_t j = (i+1); j < n; ++j){
        // Below ith-row
        float absValue = fabs(M[j][i]);
        if ( absValue > max_absValue){
            max_absValue = absValue;
            idx_max = j;
            //
            // break; // Once found a non-zero element, break it (Not going to find the maximum one)
            //
        }
        //
	}

    //
    if (idx_max == i){ // The matrix is singular !!
        return true; // is_singular
    }

	// Add that row to the current row
	Get_VectorIncrement(M[i], M[idx_max], false); // +=
	Get_VectorIncrement(M_inv[i], M_inv[idx_max], false); // +=
	//
    return false; // may be non-singular
}
vector<vector<float> > MatrixInversion(vector<vector<float> > M){
    //
	size_t n = M.size(); // The size of the square matrix

	// Check if M is a square matrix
	if (n != M[0].size()){
		return M;
	}

    // Output matrix
    vector<vector<float> > M_inv(n,vector<float>(n, 0.0));

	// Initialize the identity matrix M_inv
	for (size_t i = 0; i < n; ++i){
		M_inv[i][i] = 1.0;
	}

	//
	// cout << "M_inv:\n";
	// printMatrix(M_inv);


	// A row vector
	vector<float> row_left;
	vector<float> row_right;
	float M_diag = 1.0;

	/*
	// Check if each element on diagonal is not zero
	// If it is zero, find a row with a non-zero element on that column and add that row to the current row
	for (size_t i = 0; i < n; ++i){
		M_diag = M[i][i];
		if (M_diag == 0.0){
			// Search for other non-zero elements
			// SolveSingularityOnDiag(M, M_inv, i);

			// Search for other non-zero elements
			float max_absValue = 0.0;
			size_t idx_max = i;
			for (size_t j = 0; j < n; ++j){
				if (j != i){ // Other that i
					float absValue = fabs(M[j][i]);
					if ( absValue > max_absValue){
						max_absValue = absValue;
						idx_max = j;
					}
				}
			}
			// Add that row to the current row
			Get_VectorIncrement(M[i], M[idx_max], false); // +=
			Get_VectorIncrement(M_inv[i], M_inv[idx_max], false); // +=
			//

		}else{
			// Fine! Nothing to do.
		}
	}
	*/

	//
	/*
	cout << "M:\n";
	printMatrix(M);
	//
	cout << "M_inv:\n";
	printMatrix(M_inv);
	*/
	//

	/*
	// Scale each row vector for normalizing each elements on the diagonal to 1.0
	for (size_t i = 0; i < n; ++i){
		float ratio;
		ratio = 1.0/M[i][i];
		//
		Get_VectorScaleUp(M[i], ratio);
		Get_VectorScaleUp(M_inv[i], ratio);
		//
		// M[i] = Get_VectorScalarMultiply(M[i], ratio);
		// M_inv[i] = Get_VectorScalarMultiply(M_inv[i], ratio);
		//
	}
	*/

	/*
	//
	cout << "M:\n";
	printMatrix(M);
	//
	cout << "M_inv:\n";
	printMatrix(M_inv);
	//
	*/

    //
    // The flag for singularity
    bool is_singular = false;
    //

	// Eliminate all the other elements except those on the diagonal
	for (size_t i = 0; i < n; ++i){ // For each element on diagonal
		float ratio;

		// For solving the problem of that diagonal element is zero
		if (M[i][i] == 0.0){
			is_singular |= SolveSingularityOnDiag(M, M_inv, i);
		}

        //
        if (is_singular){
            // The matrix is singular,
            // stop this meaningless process and return.
            // break;
            return M_inv;
        }

		// Normalize the diagonal element in M to 1.0
		if (M[i][i] != 1.0){
			ratio = 1.0/M[i][i];
			Get_VectorScaleUp(M[i], ratio);
			Get_VectorScaleUp(M_inv[i], ratio);
		}
		//

		//
		row_left = M[i];
		row_right = M_inv[i];
		//
		// M_diag = M[i][i]; // This has been normalized, which is 1.0
		//
		for (size_t j = 0; j < n; ++j){ // For the row other than than the current row
			if (j != i){ // Not going to do anything with that row itself
				//
				// ratio = M[j][i]/M_diag;
				ratio = M[j][i]; // Because M_diag = 1.0
				Get_VectorIncrement(M[j], Get_VectorScalarMultiply(row_left, ratio), true); // -=
				Get_VectorIncrement(M_inv[j], Get_VectorScalarMultiply(row_right, ratio), true); // -=
                //
                /*
                Get_VectorIncrement(M[j], Get_VectorScalarMultiply(M[i], ratio), true); // -=
                Get_VectorIncrement(M_inv[j], Get_VectorScalarMultiply(M_inv[i], ratio), true); // -=
                */
			}
		}

		/*
		//
		cout << "Diagonal element i = " << i << "\n";
		//
		cout << "M:\n";
		printMatrix(M);
		//
		cout << "M_inv:\n";
		printMatrix(M_inv);
		*/
	}
	return M_inv;

}
///////////////////////////////////////// end Matrix inversion, using Gauss method


}
////////////////////////// end Namespace MATRIX_PRIMITIVE
