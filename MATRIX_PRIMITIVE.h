#ifndef MATRIX_PRIMITIVE_H
#define MATRIX_PRIMITIVE_H
//
#include <vector>
#include <math.h>

using std::vector;

// Namespace MATRIX_PRIMITIVE
//////////////////////////
namespace MATRIX_PRIMITIVE
{

// Utilities
///////////////////////////
void Mat_multiply_Vec(vector<float> &v_out, vector<vector<float> > &m_left, vector<float> &v_right); // v_out = m_left*v_right
vector<float> Mat_multiply_Vec(vector<vector<float> > &m_left, vector<float> &v_right); // v_out = m_left*v_right

// New function for matrix multiply matrix
vector<vector<float> > Mat_multiply_Mat(vector<vector<float> > &m_left, vector<vector<float> > &m_right); // m_out = m_left*m_right
vector<float> Get_VectorPlus(const vector<float> &v_a, const vector<float> &v_b, bool is_minus); // v_a + (or -) v_b
vector<float> Get_VectorScalarMultiply(const vector<float> &v_a, float scale); // scale*v_a

// Important!
// New function for scale-up a vector
// v_a *= scale
void Get_VectorScaleUp(vector<float> &v_a, float scale); // v_a *= scale

// Increment
void Get_VectorIncrement(vector<float> &v_a, const vector<float> &v_b, bool is_minus); // v_a += (or -=) v_b
/////////////////////////// end Utilities


// Matrix inversion, using Gauss method
/////////////////////////////////////////
bool SolveSingularityOnDiag(vector<vector<float> > &M, vector<vector<float> > &M_inv, size_t i);

vector<vector<float> > MatrixInversion(vector<vector<float> > M); // Note: we need a copy of M, don't use reference
///////////////////////////////////////// end Matrix inversion, using Gauss method

}
////////////////////////// end Namespace MATRIX_PRIMITIVE


#endif
