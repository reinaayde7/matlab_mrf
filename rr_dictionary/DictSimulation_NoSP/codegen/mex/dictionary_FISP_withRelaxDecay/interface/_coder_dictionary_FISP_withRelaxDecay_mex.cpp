//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// _coder_dictionary_FISP_withRelaxDecay_mex.cpp
//
// Code generation for function '_coder_dictionary_FISP_withRelaxDecay_mex'
//

// Include files
#include "_coder_dictionary_FISP_withRelaxDecay_mex.h"
#include "_coder_dictionary_FISP_withRelaxDecay_api.h"
#include "dictionary_FISP_withRelaxDecay_data.h"
#include "dictionary_FISP_withRelaxDecay_initialize.h"
#include "dictionary_FISP_withRelaxDecay_terminate.h"
#include "rt_nonfinite.h"
#include <stdexcept>

void emlrtExceptionBridge();
void emlrtExceptionBridge()
{
  throw std::runtime_error("");
}
// Function Definitions
void dictionary_FISP_withRelaxDecay_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                                int32_T nrhs,
                                                const mxArray *prhs[11])
{
  emlrtStack st{
      nullptr, // site
      nullptr, // tls
      nullptr  // prev
  };
  const mxArray *outputs[3];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  // Check for proper number of arguments.
  if (nrhs != 11) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 11, 4,
                        30, "dictionary_FISP_withRelaxDecay");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 30,
                        "dictionary_FISP_withRelaxDecay");
  }
  // Call the function.
  c_dictionary_FISP_withRelaxDeca(prhs, nlhs, outputs);
  // Copy over outputs to the caller.
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&dictionary_FISP_withRelaxDecay_atexit);
  // Module initialization.
  dictionary_FISP_withRelaxDecay_initialize();
  try { // Dispatch the entry-point.
    dictionary_FISP_withRelaxDecay_mexFunction(nlhs, plhs, nrhs, prhs);
    // Module termination.
    dictionary_FISP_withRelaxDecay_terminate();
  } catch (...) {
    emlrtCleanupOnException((emlrtCTX *)emlrtRootTLSGlobal);
    throw;
  }
}

emlrtCTX mexFunctionCreateRootTLS()
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, nullptr, 1,
                           (void *)&emlrtExceptionBridge, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

// End of code generation (_coder_dictionary_FISP_withRelaxDecay_mex.cpp)
