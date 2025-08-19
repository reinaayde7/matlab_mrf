//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// mtimes.cpp
//
// Code generation for function 'mtimes'
//

// Include files
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "blas.h"
#include "coder_array.h"
#include <cstddef>

// Variable Definitions
static emlrtRSInfo u_emlrtRSI{
    142,      // lineNo
    "mtimes", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\eml\\+coder\\+internal\\+"
    "blas\\mtimes.m" // pathName
};

static emlrtRSInfo v_emlrtRSI{
    178,           // lineNo
    "mtimes_blas", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\eml\\+coder\\+internal\\+"
    "blas\\mtimes.m" // pathName
};

static emlrtRTEInfo ab_emlrtRTEI{
    218,      // lineNo
    20,       // colNo
    "mtimes", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\eml\\+coder\\+internal\\+"
    "blas\\mtimes.m" // pName
};

static emlrtRTEInfo bb_emlrtRTEI{
    140,      // lineNo
    5,        // colNo
    "mtimes", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\eml\\+coder\\+internal\\+"
    "blas\\mtimes.m" // pName
};

// Function Definitions
namespace coder {
namespace internal {
namespace blas {
void mtimes(const emlrtStack &sp, const real_T A[9], const array<real_T, 2U> &B,
            array<real_T, 2U> &C)
{
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  emlrtStack b_st;
  emlrtStack st;
  real_T alpha1;
  real_T beta1;
  char_T TRANSA1;
  char_T TRANSB1;
  st.prev = &sp;
  st.tls = sp.tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  if (B.size(1) == 0) {
    C.set_size(&bb_emlrtRTEI, &sp, 3, 0);
  } else {
    st.site = &u_emlrtRSI;
    b_st.site = &v_emlrtRSI;
    TRANSB1 = 'N';
    TRANSA1 = 'N';
    alpha1 = 1.0;
    beta1 = 0.0;
    m_t = (ptrdiff_t)3;
    n_t = (ptrdiff_t)B.size(1);
    k_t = (ptrdiff_t)3;
    lda_t = (ptrdiff_t)3;
    ldb_t = (ptrdiff_t)3;
    ldc_t = (ptrdiff_t)3;
    C.set_size(&ab_emlrtRTEI, &b_st, 3, B.size(1));
    dgemm(&TRANSA1, &TRANSB1, &m_t, &n_t, &k_t, &alpha1, (real_T *)&A[0],
          &lda_t, &(((array<real_T, 2U> *)&B)->data())[0], &ldb_t, &beta1,
          &(C.data())[0], &ldc_t);
  }
}

} // namespace blas
} // namespace internal
} // namespace coder

// End of code generation (mtimes.cpp)
