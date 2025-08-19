//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// throt.cpp
//
// Code generation for function 'throt'
//

// Include files
#include "throt.h"
#include "rt_nonfinite.h"
#include "warning.h"
#include "mwmathutil.h"
#include <algorithm>

// Variable Definitions
static emlrtRSInfo o_emlrtRSI{
    5,       // lineNo
    "throt", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\throt.m" // pathName
};

static emlrtRSInfo p_emlrtRSI{
    31,    // lineNo
    "inv", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\matfun\\inv.m" // pathName
};

static emlrtRSInfo q_emlrtRSI{
    42,          // lineNo
    "checkcond", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\matfun\\inv.m" // pathName
};

static emlrtRSInfo r_emlrtRSI{
    46,          // lineNo
    "checkcond", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\matfun\\inv.m" // pathName
};

static emlrtMCInfo c_emlrtMCI{
    53,        // lineNo
    19,        // colNo
    "flt2str", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\flt2str.m" // pName
};

static emlrtRSInfo gb_emlrtRSI{
    53,        // lineNo
    "flt2str", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\shared\\coder\\coder\\lib\\+coder\\+"
    "internal\\flt2str.m" // pathName
};

// Function Declarations
static void b_emlrt_marshallIn(const emlrtStack &sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, char_T ret[14]);

static const mxArray *b_sprintf(const emlrtStack &sp, const mxArray *m1,
                                const mxArray *m2, emlrtMCInfo &location);

static void emlrt_marshallIn(const emlrtStack &sp,
                             const mxArray *a__output_of_sprintf_,
                             const char_T *identifier, char_T y[14]);

static void emlrt_marshallIn(const emlrtStack &sp, const mxArray *u,
                             const emlrtMsgIdentifier *parentId, char_T y[14]);

// Function Definitions
static void b_emlrt_marshallIn(const emlrtStack &sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, char_T ret[14])
{
  static const int32_T dims[2]{1, 14};
  emlrtCheckBuiltInR2012b((emlrtConstCTX)&sp, msgId, src, "char", false, 2U,
                          (const void *)&dims[0]);
  emlrtImportCharArrayR2015b((emlrtConstCTX)&sp, src, &ret[0], 14);
  emlrtDestroyArray(&src);
}

static const mxArray *b_sprintf(const emlrtStack &sp, const mxArray *m1,
                                const mxArray *m2, emlrtMCInfo &location)
{
  const mxArray *pArrays[2];
  const mxArray *m;
  pArrays[0] = m1;
  pArrays[1] = m2;
  return emlrtCallMATLABR2012b((emlrtConstCTX)&sp, 1, &m, 2, &pArrays[0],
                               "sprintf", true, &location);
}

static void emlrt_marshallIn(const emlrtStack &sp,
                             const mxArray *a__output_of_sprintf_,
                             const char_T *identifier, char_T y[14])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = const_cast<const char_T *>(identifier);
  thisId.fParent = nullptr;
  thisId.bParentIsCell = false;
  emlrt_marshallIn(sp, emlrtAlias(a__output_of_sprintf_), &thisId, y);
  emlrtDestroyArray(&a__output_of_sprintf_);
}

static void emlrt_marshallIn(const emlrtStack &sp, const mxArray *u,
                             const emlrtMsgIdentifier *parentId, char_T y[14])
{
  b_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

void throt(const emlrtStack &sp, real_T phi, real_T theta, real_T Rth[9])
{
  static const int32_T iv[2]{1, 6};
  static const char_T rfmt[6]{'%', '1', '4', '.', '6', 'e'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *m;
  const mxArray *y;
  real_T Rz[9];
  real_T a[9];
  real_T b_a[9];
  real_T x[9];
  real_T Rx_tmp;
  real_T b_Rx_tmp;
  real_T n1x;
  real_T t1;
  real_T t2;
  int32_T itmp;
  int32_T p1;
  int32_T p2;
  boolean_T exitg1;
  st.prev = &sp;
  st.tls = sp.tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  t1 = muDoubleScalarSin(-theta);
  n1x = muDoubleScalarCos(-theta);
  Rz[0] = n1x;
  Rz[3] = -t1;
  Rz[6] = 0.0;
  Rz[1] = t1;
  Rz[4] = n1x;
  Rz[7] = 0.0;
  Rz[2] = 0.0;
  Rz[5] = 0.0;
  Rz[8] = 1.0;
  Rx_tmp = muDoubleScalarSin(phi);
  b_Rx_tmp = muDoubleScalarCos(phi);
  st.site = &o_emlrtRSI;
  std::copy(&Rz[0], &Rz[9], &x[0]);
  p1 = 0;
  p2 = 3;
  itmp = 6;
  if (muDoubleScalarAbs(t1) > muDoubleScalarAbs(n1x)) {
    p1 = 3;
    p2 = 0;
    x[0] = t1;
    x[1] = n1x;
    x[3] = n1x;
    x[4] = -t1;
    x[6] = 0.0;
    x[7] = 0.0;
  }
  x[1] /= x[0];
  x[2] /= x[0];
  x[4] -= x[1] * x[3];
  x[5] -= x[2] * x[3];
  x[7] -= x[1] * x[6];
  x[8] -= x[2] * x[6];
  if (muDoubleScalarAbs(x[5]) > muDoubleScalarAbs(x[4])) {
    itmp = p2;
    p2 = 6;
    t1 = x[1];
    x[1] = x[2];
    x[2] = t1;
    t1 = x[4];
    x[4] = x[5];
    x[5] = t1;
    t1 = x[7];
    x[7] = x[8];
    x[8] = t1;
  }
  x[5] /= x[4];
  x[8] -= x[5] * x[7];
  t1 = (x[1] * x[5] - x[2]) / x[8];
  t2 = -(x[1] + x[7] * t1) / x[4];
  a[p1] = ((1.0 - x[3] * t2) - x[6] * t1) / x[0];
  a[p1 + 1] = t2;
  a[p1 + 2] = t1;
  t1 = -x[5] / x[8];
  t2 = (1.0 - x[7] * t1) / x[4];
  a[p2] = -(x[3] * t2 + x[6] * t1) / x[0];
  a[p2 + 1] = t2;
  a[p2 + 2] = t1;
  t1 = 1.0 / x[8];
  t2 = -x[7] * t1 / x[4];
  a[itmp] = -(x[3] * t2 + x[6] * t1) / x[0];
  a[itmp + 1] = t2;
  a[itmp + 2] = t1;
  b_st.site = &p_emlrtRSI;
  n1x = 0.0;
  p1 = 0;
  exitg1 = false;
  while ((!exitg1) && (p1 < 3)) {
    t1 = (muDoubleScalarAbs(Rz[3 * p1]) + muDoubleScalarAbs(Rz[3 * p1 + 1])) +
         muDoubleScalarAbs(Rz[3 * p1 + 2]);
    if (muDoubleScalarIsNaN(t1)) {
      n1x = rtNaN;
      exitg1 = true;
    } else {
      if (t1 > n1x) {
        n1x = t1;
      }
      p1++;
    }
  }
  t2 = 0.0;
  p1 = 0;
  exitg1 = false;
  while ((!exitg1) && (p1 < 3)) {
    t1 = (muDoubleScalarAbs(a[3 * p1]) + muDoubleScalarAbs(a[3 * p1 + 1])) +
         muDoubleScalarAbs(a[3 * p1 + 2]);
    if (muDoubleScalarIsNaN(t1)) {
      t2 = rtNaN;
      exitg1 = true;
    } else {
      if (t1 > t2) {
        t2 = t1;
      }
      p1++;
    }
  }
  t1 = 1.0 / (n1x * t2);
  if ((n1x == 0.0) || (t2 == 0.0) || (t1 == 0.0)) {
    c_st.site = &q_emlrtRSI;
    coder::internal::warning(c_st);
  } else if (muDoubleScalarIsNaN(t1) || (t1 < 2.2204460492503131E-16)) {
    char_T str[14];
    c_st.site = &r_emlrtRSI;
    y = nullptr;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a(&c_st, 6, m, &rfmt[0]);
    emlrtAssign(&y, m);
    b_y = nullptr;
    m = emlrtCreateDoubleScalar(t1);
    emlrtAssign(&b_y, m);
    d_st.site = &gb_emlrtRSI;
    emlrt_marshallIn(d_st, b_sprintf(d_st, y, b_y, c_emlrtMCI),
                     "<output of sprintf>", str);
    c_st.site = &r_emlrtRSI;
    coder::internal::warning(c_st, str);
  }
  x[0] = 1.0;
  x[3] = 0.0;
  x[6] = 0.0;
  x[1] = 0.0;
  x[4] = b_Rx_tmp;
  x[7] = -Rx_tmp;
  x[2] = 0.0;
  x[5] = Rx_tmp;
  x[8] = b_Rx_tmp;
  for (p1 = 0; p1 < 3; p1++) {
    t1 = a[p1];
    n1x = a[p1 + 3];
    t2 = a[p1 + 6];
    for (p2 = 0; p2 < 3; p2++) {
      b_a[p1 + 3 * p2] =
          (t1 * x[3 * p2] + n1x * x[3 * p2 + 1]) + t2 * x[3 * p2 + 2];
    }
    t1 = b_a[p1];
    n1x = b_a[p1 + 3];
    t2 = b_a[p1 + 6];
    for (p2 = 0; p2 < 3; p2++) {
      Rth[p1 + 3 * p2] =
          (t1 * Rz[3 * p2] + n1x * Rz[3 * p2 + 1]) + t2 * Rz[3 * p2 + 2];
    }
  }
}

// End of code generation (throt.cpp)
