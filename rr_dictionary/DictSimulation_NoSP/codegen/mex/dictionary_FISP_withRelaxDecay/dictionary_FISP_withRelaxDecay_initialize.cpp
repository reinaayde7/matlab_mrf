//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// dictionary_FISP_withRelaxDecay_initialize.cpp
//
// Code generation for function 'dictionary_FISP_withRelaxDecay_initialize'
//

// Include files
#include "dictionary_FISP_withRelaxDecay_initialize.h"
#include "_coder_dictionary_FISP_withRelaxDecay_mex.h"
#include "dictionary_FISP_withRelaxDecay_data.h"
#include "rt_nonfinite.h"

// Function Declarations
static void dictionary_FISP_withRelaxDecay_once();

// Function Definitions
static void dictionary_FISP_withRelaxDecay_once()
{
  mex_InitInfAndNan();
}

void dictionary_FISP_withRelaxDecay_initialize()
{
  emlrtStack st{
      nullptr, // site
      nullptr, // tls
      nullptr  // prev
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, nullptr);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    dictionary_FISP_withRelaxDecay_once();
  }
}

// End of code generation (dictionary_FISP_withRelaxDecay_initialize.cpp)
