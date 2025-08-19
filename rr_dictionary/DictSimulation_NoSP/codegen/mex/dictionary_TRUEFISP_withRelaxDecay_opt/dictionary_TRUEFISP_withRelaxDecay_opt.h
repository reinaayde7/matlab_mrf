//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// dictionary_TRUEFISP_withRelaxDecay_opt.h
//
// Code generation for function 'dictionary_TRUEFISP_withRelaxDecay_opt'
//

#pragma once

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include "emlrt.h"
#include "mex.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Function Declarations
void dictionary_TRUEFISP_withRelaxDecay_opt(
    const emlrtStack *sp, const coder::array<real_T, 2U> &r,
    const coder::array<real_T, 1U> &flip, const coder::array<real_T, 1U> &tr,
    real_T te, const coder::array<real_T, 1U> &phase, real_T ti,
    real_T phasetwist, real_T spins, real_T frames, real_T Inv,
    real_T waitingDuration, coder::array<real32_T, 2U> &Mx,
    coder::array<real32_T, 2U> &My, coder::array<real32_T, 2U> &Mz);

// End of code generation (dictionary_TRUEFISP_withRelaxDecay_opt.h)
