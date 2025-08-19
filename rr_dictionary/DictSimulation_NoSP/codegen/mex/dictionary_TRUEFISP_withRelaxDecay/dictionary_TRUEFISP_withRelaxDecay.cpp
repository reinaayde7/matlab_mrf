//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// dictionary_TRUEFISP_withRelaxDecay.cpp
//
// Code generation for function 'dictionary_TRUEFISP_withRelaxDecay'
//

// Include files
#include "dictionary_TRUEFISP_withRelaxDecay.h"
#include "dictionary_TRUEFISP_withRelaxDecay_data.h"
#include "eml_int_forloop_overflow_check.h"
#include "mean.h"
#include "mtimes.h"
#include "rt_nonfinite.h"
#include "throt.h"
#include "coder_array.h"
#include "mwmathutil.h"
#include <emmintrin.h>

// Variable Definitions
static emlrtRSInfo emlrtRSI{
    37,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo b_emlrtRSI{
    38,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo c_emlrtRSI{
    46,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo d_emlrtRSI{
    50,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo e_emlrtRSI{
    63,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo f_emlrtRSI{
    66,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo g_emlrtRSI{
    71,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo h_emlrtRSI{
    72,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo i_emlrtRSI{
    73,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo j_emlrtRSI{
    80,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo k_emlrtRSI{
    94,                                   // lineNo
    "dictionary_TRUEFISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pathName
};

static emlrtRSInfo l_emlrtRSI{
    38,       // lineNo
    "repmat", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\elmat\\repmat.m" // pathName
};

static emlrtRSInfo m_emlrtRSI{
    74,       // lineNo
    "repmat", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\elmat\\repmat.m" // pathName
};

static emlrtRSInfo
    s_emlrtRSI{
        94,                  // lineNo
        "eml_mtimes_helper", // fcnName
        "C:\\Program "
        "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_"
        "helper.m" // pathName
    };

static emlrtDCInfo emlrtDCI{
    27,                                   // lineNo
    17,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    4               // checkKind
};

static emlrtDCInfo b_emlrtDCI{
    27,                                   // lineNo
    17,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtRTEInfo emlrtRTEI{
    42,                                   // lineNo
    17,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtECInfo emlrtECI{
    2,                                    // nDims
    50,                                   // lineNo
    21,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtECInfo b_emlrtECI{
    2,                                    // nDims
    66,                                   // lineNo
    17,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtECInfo c_emlrtECI{
    2,                                    // nDims
    80,                                   // lineNo
    15,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtECInfo d_emlrtECI{
    2,                                    // nDims
    94,                                   // lineNo
    25,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtECInfo e_emlrtECI{
    -1,                                   // nDims
    94,                                   // lineNo
    9,                                    // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo b_emlrtRTEI{
    58,                   // lineNo
    23,                   // colNo
    "assertValidSizeArg", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\assertValidSizeArg.m" // pName
};

static emlrtRTEInfo c_emlrtRTEI{
    64,                   // lineNo
    15,                   // colNo
    "assertValidSizeArg", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\eml\\+coder\\+"
    "internal\\assertValidSizeArg.m" // pName
};

static emlrtDCInfo c_emlrtDCI{
    16,                                   // lineNo
    12,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtDCInfo d_emlrtDCI{
    16,                                   // lineNo
    12,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    4               // checkKind
};

static emlrtDCInfo e_emlrtDCI{
    17,                                   // lineNo
    12,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtDCInfo f_emlrtDCI{
    18,                                   // lineNo
    12,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtDCInfo g_emlrtDCI{
    16,                                   // lineNo
    1,                                    // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtDCInfo h_emlrtDCI{
    17,                                   // lineNo
    1,                                    // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtDCInfo i_emlrtDCI{
    18,                                   // lineNo
    1,                                    // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtBCInfo emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    33,                                   // lineNo
    12,                                   // colNo
    "r",                                  // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtBCInfo b_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    34,                                   // lineNo
    12,                                   // colNo
    "r",                                  // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtBCInfo c_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    35,                                   // lineNo
    12,                                   // colNo
    "r",                                  // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtBCInfo d_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    63,                                   // lineNo
    28,                                   // colNo
    "flip",                               // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtBCInfo e_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    63,                                   // lineNo
    37,                                   // colNo
    "phase",                              // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtBCInfo f_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    71,                                   // lineNo
    16,                                   // colNo
    "Mx",                                 // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtDCInfo j_emlrtDCI{
    71,                                   // lineNo
    16,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtBCInfo g_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    71,                                   // lineNo
    31,                                   // colNo
    "Mx",                                 // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtBCInfo h_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    72,                                   // lineNo
    16,                                   // colNo
    "My",                                 // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtDCInfo k_emlrtDCI{
    72,                                   // lineNo
    16,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtBCInfo i_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    72,                                   // lineNo
    31,                                   // colNo
    "My",                                 // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtBCInfo j_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    73,                                   // lineNo
    16,                                   // colNo
    "Mz",                                 // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtDCInfo l_emlrtDCI{
    73,                                   // lineNo
    16,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    1               // checkKind
};

static emlrtBCInfo k_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    73,                                   // lineNo
    31,                                   // colNo
    "Mz",                                 // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtBCInfo l_emlrtBCI{
    -1,                                   // iFirst
    -1,                                   // iLast
    78,                                   // lineNo
    34,                                   // colNo
    "tr",                                 // aName
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m", // pName
    0               // checkKind
};

static emlrtRTEInfo e_emlrtRTEI{
    16,                                   // lineNo
    1,                                    // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo f_emlrtRTEI{
    17,                                   // lineNo
    1,                                    // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo g_emlrtRTEI{
    18,                                   // lineNo
    1,                                    // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo h_emlrtRTEI{
    69,       // lineNo
    28,       // colNo
    "repmat", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\elmat\\repmat.m" // pName
};

static emlrtRTEInfo i_emlrtRTEI{
    46,                                   // lineNo
    37,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo j_emlrtRTEI{
    49,                                   // lineNo
    17,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo k_emlrtRTEI{
    50,                                   // lineNo
    17,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo l_emlrtRTEI{
    63,                                   // lineNo
    41,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo m_emlrtRTEI{
    65,                                   // lineNo
    13,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo n_emlrtRTEI{
    66,                                   // lineNo
    13,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo o_emlrtRTEI{
    71,                                   // lineNo
    41,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo p_emlrtRTEI{
    72,                                   // lineNo
    41,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo q_emlrtRTEI{
    73,                                   // lineNo
    41,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo r_emlrtRTEI{
    79,                                   // lineNo
    13,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo s_emlrtRTEI{
    80,                                   // lineNo
    13,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo t_emlrtRTEI{
    93,                                   // lineNo
    9,                                    // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

static emlrtRTEInfo u_emlrtRTEI{
    94,                                   // lineNo
    25,                                   // colNo
    "dictionary_TRUEFISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_TRUEFISP_with"
    "RelaxDecay.m" // pName
};

// Function Declarations
static void plus(const emlrtStack &sp, coder::array<real_T, 2U> &in1,
                 const coder::array<real_T, 2U> &in2,
                 const coder::array<real_T, 2U> &in3);

static void plus(const emlrtStack &sp, coder::array<real_T, 2U> &in1,
                 const coder::array<real_T, 2U> &in2);

// Function Definitions
static void plus(const emlrtStack &sp, coder::array<real_T, 2U> &in1,
                 const coder::array<real_T, 2U> &in2,
                 const coder::array<real_T, 2U> &in3)
{
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in1.set_size(&s_emlrtRTEI, &sp, 3, in1.size(1));
  if (in3.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(1);
  }
  in1.set_size(&s_emlrtRTEI, &sp, in1.size(0), loop_ub);
  stride_0_1 = (in2.size(1) != 1);
  stride_1_1 = (in3.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (int32_T i{0}; i < loop_ub; i++) {
    __m128d r;
    __m128d r1;
    r = _mm_loadu_pd(&in2[3 * aux_0_1]);
    r1 = _mm_loadu_pd(&in3[3 * aux_1_1]);
    _mm_storeu_pd(&in1[3 * i], _mm_add_pd(r, r1));
    in1[3 * i + 2] = in2[3 * aux_0_1 + 2] + in3[3 * aux_1_1 + 2];
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

static void plus(const emlrtStack &sp, coder::array<real_T, 2U> &in1,
                 const coder::array<real_T, 2U> &in2)
{
  coder::array<real_T, 2U> b_in1;
  int32_T aux_0_1;
  int32_T aux_1_1;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)&sp);
  if (in2.size(1) == 1) {
    loop_ub = in1.size(1);
  } else {
    loop_ub = in2.size(1);
  }
  b_in1.set_size(&u_emlrtRTEI, &sp, 3, loop_ub);
  stride_0_1 = (in1.size(1) != 1);
  stride_1_1 = (in2.size(1) != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (int32_T i{0}; i < loop_ub; i++) {
    __m128d r;
    __m128d r1;
    r = _mm_loadu_pd(&in1[3 * aux_0_1]);
    r1 = _mm_loadu_pd(&in2[3 * aux_1_1]);
    _mm_storeu_pd(&b_in1[3 * i], _mm_add_pd(r, r1));
    b_in1[3 * i + 2] = in1[3 * aux_0_1 + 2] + in2[3 * aux_1_1 + 2];
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  in1.set_size(&u_emlrtRTEI, &sp, 3, in1.size(1));
  loop_ub = b_in1.size(1);
  in1.set_size(&u_emlrtRTEI, &sp, in1.size(0), b_in1.size(1));
  for (int32_T i{0}; i < loop_ub; i++) {
    in1[3 * i] = b_in1[3 * i];
    in1[3 * i + 1] = b_in1[3 * i + 1];
    in1[3 * i + 2] = b_in1[3 * i + 2];
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)&sp);
}

void dictionary_TRUEFISP_withRelaxDecay(
    const emlrtStack *sp, const coder::array<real_T, 2U> &r,
    const coder::array<real_T, 1U> &flip, const coder::array<real_T, 1U> &tr,
    const real_T te[2], const coder::array<real_T, 1U> &phase, real_T ti,
    real_T spins, real_T frames, real_T Inv, real_T waitingDuration,
    coder::array<real32_T, 2U> &Mx, coder::array<real32_T, 2U> &My,
    coder::array<real32_T, 2U> &Mz)
{
  coder::array<real_T, 3U> M0;
  coder::array<real_T, 2U> Bti;
  coder::array<real_T, 2U> M;
  coder::array<real_T, 2U> b_M;
  coder::array<real_T, 2U> b_r;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  real_T E2;
  real_T Rz_tmp;
  real_T b_n;
  real_T nt;
  int32_T iv[2];
  int32_T b_loop_ub;
  int32_T b_ntilecols;
  int32_T b_outsize_idx_1;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T loop_ub;
  int32_T ntilecols;
  int32_T outsize_idx_1;
  boolean_T b;
  boolean_T b1;
  boolean_T b_overflow;
  boolean_T overflow;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  // Rudy: phase twisting is commented out since bSSFP does not add spoilers.
  // We commented them out to save sim time.
  nt = flip.size(0);
  if (flip.size(0) > frames) {
    nt = frames;
  }
  // df=0;
  Rz_tmp = nt * 2.0;
  if (!(Rz_tmp >= 0.0)) {
    emlrtNonNegativeCheckR2012b(Rz_tmp, &d_emlrtDCI, (emlrtConstCTX)sp);
  }
  E2 = static_cast<int32_T>(muDoubleScalarFloor(Rz_tmp));
  if (Rz_tmp != E2) {
    emlrtIntegerCheckR2012b(Rz_tmp, &c_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = r.size(0);
  Mx.set_size(&e_emlrtRTEI, sp, static_cast<int32_T>(Rz_tmp), r.size(0));
  if (Rz_tmp != E2) {
    emlrtIntegerCheckR2012b(Rz_tmp, &g_emlrtDCI, (emlrtConstCTX)sp);
  }
  loop_ub = static_cast<int32_T>(Rz_tmp) * r.size(0);
  for (i1 = 0; i1 < loop_ub; i1++) {
    Mx[i1] = 0.0F;
  }
  if (Rz_tmp != E2) {
    emlrtIntegerCheckR2012b(Rz_tmp, &e_emlrtDCI, (emlrtConstCTX)sp);
  }
  My.set_size(&f_emlrtRTEI, sp, static_cast<int32_T>(Rz_tmp), r.size(0));
  if (Rz_tmp != E2) {
    emlrtIntegerCheckR2012b(Rz_tmp, &h_emlrtDCI, (emlrtConstCTX)sp);
  }
  loop_ub = static_cast<int32_T>(Rz_tmp) * r.size(0);
  for (i1 = 0; i1 < loop_ub; i1++) {
    My[i1] = 0.0F;
  }
  if (Rz_tmp != E2) {
    emlrtIntegerCheckR2012b(Rz_tmp, &f_emlrtDCI, (emlrtConstCTX)sp);
  }
  Mz.set_size(&g_emlrtRTEI, sp, static_cast<int32_T>(Rz_tmp), r.size(0));
  if (Rz_tmp != E2) {
    emlrtIntegerCheckR2012b(Rz_tmp, &i_emlrtDCI, (emlrtConstCTX)sp);
  }
  loop_ub = static_cast<int32_T>(Rz_tmp) * r.size(0);
  for (i1 = 0; i1 < loop_ub; i1++) {
    Mz[i1] = 0.0F;
  }
  //  no need for phase twisting in trueFISP
  //  phirange = linspace(-phasetwist/2,phasetwist/2,spins);
  //  Z = zeros(3,3,spins);
  //  for s = 1:spins
  //      Z(:,:,s) = zrot(phirange(s));
  //  end
  if (!(spins >= 0.0)) {
    emlrtNonNegativeCheckR2012b(spins, &emlrtDCI, (emlrtConstCTX)sp);
  }
  if (spins != static_cast<int32_T>(muDoubleScalarFloor(spins))) {
    emlrtIntegerCheckR2012b(spins, &b_emlrtDCI, (emlrtConstCTX)sp);
  }
  if (r.size(0) - 1 >= 0) {
    b = !muDoubleScalarIsInf(spins);
    if (b && (!(spins > 2.147483647E+9))) {
      b1 = true;
    } else {
      b1 = false;
    }
    outsize_idx_1 = static_cast<int32_T>(spins);
    ntilecols = static_cast<int32_T>(spins) * 3;
    overflow = (ntilecols > 2147483646);
    if (b && (!(spins > 2.147483647E+9))) {
      b = true;
    } else {
      b = false;
    }
    if (spins <= 0.0) {
      b_n = 0.0;
    } else {
      b_n = spins;
    }
    b_outsize_idx_1 = static_cast<int32_T>(spins);
    b_ntilecols = static_cast<int32_T>(spins);
    b_overflow = (static_cast<int32_T>(spins) > 2147483646);
    i2 = static_cast<int32_T>(nt);
  }
  if (i - 1 >= 0) {
    b_loop_ub = static_cast<int32_T>(spins);
  }
  for (int32_T n{0}; n < i; n++) {
    real_T E1;
    real_T T1;
    real_T T2;
    real_T c_E2;
    real_T c_Rz_tmp;
    real_T d_Rz_tmp;
    real_T df;
    real_T phi;
    int32_T ibtile;
    int32_T jtilecol;
    // Rudy debug
    if (n + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(n + 1, 1, i, &emlrtBCI, (emlrtConstCTX)sp);
    }
    T1 = r[n] / 1000.0;
    if (n + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(n + 1, 1, i, &b_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    T2 = r[n + r.size(0)] / 1000.0;
    if (n + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(n + 1, 1, i, &c_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    df = r[n + r.size(0) * 2];
    st.site = &emlrtRSI;
    b_st.site = &l_emlrtRSI;
    if (!b1) {
      emlrtErrorWithMessageIdR2018a(
          &b_st, &b_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
          "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
    }
    if (spins <= 0.0) {
      Rz_tmp = 0.0;
    } else {
      Rz_tmp = spins;
    }
    if (!(Rz_tmp * 3.0 <= 2.147483647E+9)) {
      emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
                                    "Coder:MATLAB:pmaxsize",
                                    "Coder:MATLAB:pmaxsize", 0);
    }
    M0.set_size(&h_emlrtRTEI, &st, 3, outsize_idx_1, 3);
    b_st.site = &m_emlrtRSI;
    if (overflow) {
      c_st.site = &n_emlrtRSI;
      coder::check_forloop_overflow_error(c_st);
    }
    for (jtilecol = 0; jtilecol < ntilecols; jtilecol++) {
      ibtile = jtilecol * 3;
      M0[ibtile] = 0.0;
      M0[ibtile + 1] = 0.0;
      M0[ibtile + 2] = 1.0;
    }
    st.site = &b_emlrtRSI;
    b_st.site = &l_emlrtRSI;
    if (!b) {
      emlrtErrorWithMessageIdR2018a(
          &b_st, &b_emlrtRTEI, "Coder:MATLAB:NonIntegerInput",
          "Coder:MATLAB:NonIntegerInput", 4, 12, MIN_int32_T, 12, MAX_int32_T);
    }
    if (!(b_n <= 2.147483647E+9)) {
      emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
                                    "Coder:MATLAB:pmaxsize",
                                    "Coder:MATLAB:pmaxsize", 0);
    }
    M.set_size(&h_emlrtRTEI, &st, 3, b_outsize_idx_1);
    b_st.site = &m_emlrtRSI;
    if (b_overflow) {
      c_st.site = &n_emlrtRSI;
      coder::check_forloop_overflow_error(c_st);
    }
    for (jtilecol = 0; jtilecol < b_ntilecols; jtilecol++) {
      ibtile = jtilecol * 3;
      M[ibtile] = 0.0;
      M[ibtile + 1] = 0.0;
      M[ibtile + 2] = 0.0;
    }
    phi = 6.2831853071795862 * df * waitingDuration;
    E1 = muDoubleScalarExp(-waitingDuration / T1);
    c_E2 = muDoubleScalarExp(-waitingDuration / T2);
    c_Rz_tmp = muDoubleScalarSin(phi);
    d_Rz_tmp = muDoubleScalarCos(phi);
    for (int32_T iEcho{0}; iEcho < 2; iEcho++) {
      __m128d r1;
      __m128d r2;
      real_T a[9];
      real_T b_E2[9];
      real_T b_Rz_tmp[9];
      // Rudy debug
      // iEcho
      emlrtForLoopVectorCheckR2021a(1.0, 1.0, nt, mxDOUBLE_CLASS,
                                    static_cast<int32_T>(nt), &emlrtRTEI,
                                    (emlrtConstCTX)sp);
      for (int32_T k{0}; k < i2; k++) {
        real_T b_E1;
        int32_T c_loop_ub;
        // Rudy debug
        //  disp(k)  %Rudy debug
        if (static_cast<uint32_T>(k) + 1U == 1U) {
          if (Inv == 1.0) {
            st.site = &c_emlrtRSI;
            b_st.site = &c_emlrtRSI;
            throt(b_st, 3.1415926535897931, 0.0, a);
            ibtile = M0.size(1);
            Bti.set_size(&i_emlrtRTEI, &st, 3, M0.size(1));
            for (i1 = 0; i1 < ibtile; i1++) {
              Bti[3 * i1] = M0[3 * i1 + 3 * M0.size(1) * iEcho];
              Bti[3 * i1 + 1] = M0[(3 * i1 + 3 * M0.size(1) * iEcho) + 1];
              Bti[3 * i1 + 2] = M0[(3 * i1 + 3 * M0.size(1) * iEcho) + 2];
            }
            b_st.site = &s_emlrtRSI;
            coder::internal::blas::mtimes(b_st, a, Bti, M);
          }
          //
          // 	Function simulates free precession and decay
          // 	over a time interval T, given relaxation times T1 and T2
          // 	and off-resonance df.  Times in ms, off-resonance in Hz.
          phi = 6.2831853071795862 * df * ti;
          //  Resonant precession, radians.
          b_E1 = muDoubleScalarExp(-ti / T1);
          E2 = muDoubleScalarExp(-ti / T2);
          Rz_tmp = muDoubleScalarSin(phi);
          phi = muDoubleScalarCos(phi);
          loop_ub = static_cast<int32_T>(spins);
          Bti.set_size(&j_emlrtRTEI, sp, 3, static_cast<int32_T>(spins));
          for (i1 = 0; i1 < loop_ub; i1++) {
            Bti[3 * i1] = 0.0;
            Bti[3 * i1 + 1] = 0.0;
            Bti[3 * i1 + 2] = 1.0 - b_E1;
          }
          st.site = &d_emlrtRSI;
          a[0] = E2;
          a[3] = 0.0;
          a[6] = 0.0;
          a[1] = 0.0;
          a[4] = E2;
          a[7] = 0.0;
          a[2] = 0.0;
          a[5] = 0.0;
          a[8] = b_E1;
          b_Rz_tmp[0] = phi;
          b_Rz_tmp[3] = -Rz_tmp;
          b_Rz_tmp[6] = 0.0;
          b_Rz_tmp[1] = Rz_tmp;
          b_Rz_tmp[4] = phi;
          b_Rz_tmp[7] = 0.0;
          b_Rz_tmp[2] = 0.0;
          b_Rz_tmp[5] = 0.0;
          b_Rz_tmp[8] = 1.0;
          for (i1 = 0; i1 < 3; i1++) {
            Rz_tmp = a[i1];
            E2 = a[i1 + 3];
            phi = a[i1 + 6];
            for (ibtile = 0; ibtile < 3; ibtile++) {
              b_E2[i1 + 3 * ibtile] = (Rz_tmp * b_Rz_tmp[3 * ibtile] +
                                       E2 * b_Rz_tmp[3 * ibtile + 1]) +
                                      phi * b_Rz_tmp[3 * ibtile + 2];
            }
          }
          b_st.site = &s_emlrtRSI;
          coder::internal::blas::mtimes(b_st, b_E2, M, b_r);
          if ((b_r.size(1) != static_cast<int32_T>(spins)) &&
              ((b_r.size(1) != 1) && (static_cast<int32_T>(spins) != 1))) {
            emlrtDimSizeImpxCheckR2021b(b_r.size(1),
                                        static_cast<int32_T>(spins), &emlrtECI,
                                        (emlrtConstCTX)sp);
          }
          if (b_r.size(1) == Bti.size(1)) {
            M.set_size(&k_emlrtRTEI, sp, 3, b_r.size(1));
            loop_ub = 3 * b_r.size(1);
            ibtile = (loop_ub / 2) << 1;
            jtilecol = ibtile - 2;
            for (i1 = 0; i1 <= jtilecol; i1 += 2) {
              r1 = _mm_loadu_pd(&b_r[i1]);
              r2 = _mm_loadu_pd(&Bti[i1]);
              _mm_storeu_pd(&M[i1], _mm_add_pd(r1, r2));
            }
            for (i1 = ibtile; i1 < loop_ub; i1++) {
              M[i1] = b_r[i1] + Bti[i1];
            }
          } else {
            st.site = &d_emlrtRSI;
            plus(st, M, b_r, Bti);
          }
          // Rudy 240627: spoiling not used in TRUEFISP + commented out to save
          // sim time
          //  spoiling
          //  for s=1:spins
          //      M(:,s)=Z(:,:,s)*M(:,s);
          //  end
        }
        //  alpha pulse
        st.site = &e_emlrtRSI;
        if ((static_cast<int32_T>(static_cast<uint32_T>(k) + 1U) < 1) ||
            (static_cast<int32_T>(static_cast<uint32_T>(k) + 1U) >
             flip.size(0))) {
          emlrtDynamicBoundsCheckR2012b(
              static_cast<int32_T>(static_cast<uint32_T>(k) + 1U), 1,
              flip.size(0), &d_emlrtBCI, &st);
        }
        if ((static_cast<int32_T>(static_cast<uint32_T>(k) + 1U) < 1) ||
            (static_cast<int32_T>(static_cast<uint32_T>(k) + 1U) >
             phase.size(0))) {
          emlrtDynamicBoundsCheckR2012b(
              static_cast<int32_T>(static_cast<uint32_T>(k) + 1U), 1,
              phase.size(0), &e_emlrtBCI, &st);
        }
        b_st.site = &e_emlrtRSI;
        throt(b_st, flip[k], phase[k], a);
        Bti.set_size(&l_emlrtRTEI, &st, 3, M.size(1));
        loop_ub = M.size(0) * M.size(1) - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
          Bti[i1] = M[i1];
        }
        b_st.site = &s_emlrtRSI;
        coder::internal::blas::mtimes(b_st, a, Bti, M);
        //
        // 	Function simulates free precession and decay
        // 	over a time interval T, given relaxation times T1 and T2
        // 	and off-resonance df.  Times in ms, off-resonance in Hz.
        Rz_tmp = te[iEcho];
        phi = 6.2831853071795862 * df * Rz_tmp;
        //  Resonant precession, radians.
        b_E1 = muDoubleScalarExp(-Rz_tmp / T1);
        E2 = muDoubleScalarExp(-Rz_tmp / T2);
        Rz_tmp = muDoubleScalarSin(phi);
        phi = muDoubleScalarCos(phi);
        loop_ub = static_cast<int32_T>(spins);
        Bti.set_size(&m_emlrtRTEI, sp, 3, static_cast<int32_T>(spins));
        for (i1 = 0; i1 < loop_ub; i1++) {
          Bti[3 * i1] = 0.0;
          Bti[3 * i1 + 1] = 0.0;
          Bti[3 * i1 + 2] = 1.0 - b_E1;
        }
        st.site = &f_emlrtRSI;
        a[0] = E2;
        a[3] = 0.0;
        a[6] = 0.0;
        a[1] = 0.0;
        a[4] = E2;
        a[7] = 0.0;
        a[2] = 0.0;
        a[5] = 0.0;
        a[8] = b_E1;
        b_Rz_tmp[0] = phi;
        b_Rz_tmp[3] = -Rz_tmp;
        b_Rz_tmp[6] = 0.0;
        b_Rz_tmp[1] = Rz_tmp;
        b_Rz_tmp[4] = phi;
        b_Rz_tmp[7] = 0.0;
        b_Rz_tmp[2] = 0.0;
        b_Rz_tmp[5] = 0.0;
        b_Rz_tmp[8] = 1.0;
        for (i1 = 0; i1 < 3; i1++) {
          Rz_tmp = a[i1];
          E2 = a[i1 + 3];
          phi = a[i1 + 6];
          for (ibtile = 0; ibtile < 3; ibtile++) {
            b_E2[i1 + 3 * ibtile] = (Rz_tmp * b_Rz_tmp[3 * ibtile] +
                                     E2 * b_Rz_tmp[3 * ibtile + 1]) +
                                    phi * b_Rz_tmp[3 * ibtile + 2];
          }
        }
        b_st.site = &s_emlrtRSI;
        coder::internal::blas::mtimes(b_st, b_E2, M, b_r);
        if ((b_r.size(1) != static_cast<int32_T>(spins)) &&
            ((b_r.size(1) != 1) && (static_cast<int32_T>(spins) != 1))) {
          emlrtDimSizeImpxCheckR2021b(b_r.size(1), static_cast<int32_T>(spins),
                                      &b_emlrtECI, (emlrtConstCTX)sp);
        }
        if (b_r.size(1) == Bti.size(1)) {
          M.set_size(&n_emlrtRTEI, sp, 3, b_r.size(1));
          c_loop_ub = 3 * b_r.size(1);
          ibtile = (c_loop_ub / 2) << 1;
          jtilecol = ibtile - 2;
          for (i1 = 0; i1 <= jtilecol; i1 += 2) {
            r1 = _mm_loadu_pd(&b_r[i1]);
            r2 = _mm_loadu_pd(&Bti[i1]);
            _mm_storeu_pd(&M[i1], _mm_add_pd(r1, r2));
          }
          for (i1 = ibtile; i1 < c_loop_ub; i1++) {
            M[i1] = b_r[i1] + Bti[i1];
          }
        } else {
          st.site = &f_emlrtRSI;
          plus(st, M, b_r, Bti);
        }
        // Rudy: rephase signal (ADC phase == RF phase)
        //  M = zrot(-phase(k))*M;
        c_loop_ub = M.size(1);
        b_M.set_size(&o_emlrtRTEI, sp, 1, M.size(1));
        for (i1 = 0; i1 < c_loop_ub; i1++) {
          b_M[i1] = M[3 * i1];
        }
        Rz_tmp = nt * ((static_cast<real_T>(iEcho) + 1.0) - 1.0) +
                 (static_cast<real_T>(k) + 1.0);
        i1 = static_cast<int32_T>(muDoubleScalarFloor(Rz_tmp));
        if (Rz_tmp != i1) {
          emlrtIntegerCheckR2012b(Rz_tmp, &j_emlrtDCI, (emlrtConstCTX)sp);
        }
        ibtile = static_cast<int32_T>(Rz_tmp);
        if ((Rz_tmp < 1.0) || (ibtile > Mx.size(0))) {
          emlrtDynamicBoundsCheckR2012b(static_cast<int32_T>(Rz_tmp), 1,
                                        Mx.size(0), &f_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        if (n + 1 > Mx.size(1)) {
          emlrtDynamicBoundsCheckR2012b(n + 1, 1, Mx.size(1), &g_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        st.site = &g_emlrtRSI;
        Mx[(static_cast<int32_T>(Rz_tmp) + Mx.size(0) * n) - 1] =
            static_cast<real32_T>(coder::mean(st, b_M));
        b_M.set_size(&p_emlrtRTEI, sp, 1, M.size(1));
        for (jtilecol = 0; jtilecol < c_loop_ub; jtilecol++) {
          b_M[jtilecol] = M[3 * jtilecol + 1];
        }
        if (ibtile != i1) {
          emlrtIntegerCheckR2012b(Rz_tmp, &k_emlrtDCI, (emlrtConstCTX)sp);
        }
        if ((Rz_tmp < 1.0) || (ibtile > My.size(0))) {
          emlrtDynamicBoundsCheckR2012b(static_cast<int32_T>(Rz_tmp), 1,
                                        My.size(0), &h_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        if (n + 1 > My.size(1)) {
          emlrtDynamicBoundsCheckR2012b(n + 1, 1, My.size(1), &i_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        st.site = &h_emlrtRSI;
        My[(static_cast<int32_T>(Rz_tmp) + My.size(0) * n) - 1] =
            static_cast<real32_T>(coder::mean(st, b_M));
        b_M.set_size(&q_emlrtRTEI, sp, 1, M.size(1));
        for (jtilecol = 0; jtilecol < c_loop_ub; jtilecol++) {
          b_M[jtilecol] = M[3 * jtilecol + 2];
        }
        if (ibtile != i1) {
          emlrtIntegerCheckR2012b(Rz_tmp, &l_emlrtDCI, (emlrtConstCTX)sp);
        }
        if ((Rz_tmp < 1.0) || (ibtile > Mz.size(0))) {
          emlrtDynamicBoundsCheckR2012b(static_cast<int32_T>(Rz_tmp), 1,
                                        Mz.size(0), &j_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        if (n + 1 > Mz.size(1)) {
          emlrtDynamicBoundsCheckR2012b(n + 1, 1, Mz.size(1), &k_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        st.site = &i_emlrtRSI;
        Mz[(static_cast<int32_T>(Rz_tmp) + Mz.size(0) * n) - 1] =
            static_cast<real32_T>(coder::mean(st, b_M));
        // Rudy: dephase signal (RF oscillation)
        //  M = zrot(phase(k))*M;
        if ((static_cast<int32_T>(static_cast<uint32_T>(k) + 1U) < 1) ||
            (static_cast<int32_T>(static_cast<uint32_T>(k) + 1U) >
             tr.size(0))) {
          emlrtDynamicBoundsCheckR2012b(
              static_cast<int32_T>(static_cast<uint32_T>(k) + 1U), 1,
              tr.size(0), &l_emlrtBCI, (emlrtConstCTX)sp);
        }
        E2 = tr[k] - te[iEcho];
        //
        // 	Function simulates free precession and decay
        // 	over a time interval T, given relaxation times T1 and T2
        // 	and off-resonance df.  Times in ms, off-resonance in Hz.
        phi = 6.2831853071795862 * df * E2;
        //  Resonant precession, radians.
        b_E1 = muDoubleScalarExp(-E2 / T1);
        E2 = muDoubleScalarExp(-E2 / T2);
        Rz_tmp = muDoubleScalarSin(phi);
        phi = muDoubleScalarCos(phi);
        Bti.set_size(&r_emlrtRTEI, sp, 3, static_cast<int32_T>(spins));
        for (i1 = 0; i1 < loop_ub; i1++) {
          Bti[3 * i1] = 0.0;
          Bti[3 * i1 + 1] = 0.0;
          Bti[3 * i1 + 2] = 1.0 - b_E1;
        }
        st.site = &j_emlrtRSI;
        a[0] = E2;
        a[3] = 0.0;
        a[6] = 0.0;
        a[1] = 0.0;
        a[4] = E2;
        a[7] = 0.0;
        a[2] = 0.0;
        a[5] = 0.0;
        a[8] = b_E1;
        b_Rz_tmp[0] = phi;
        b_Rz_tmp[3] = -Rz_tmp;
        b_Rz_tmp[6] = 0.0;
        b_Rz_tmp[1] = Rz_tmp;
        b_Rz_tmp[4] = phi;
        b_Rz_tmp[7] = 0.0;
        b_Rz_tmp[2] = 0.0;
        b_Rz_tmp[5] = 0.0;
        b_Rz_tmp[8] = 1.0;
        for (i1 = 0; i1 < 3; i1++) {
          Rz_tmp = a[i1];
          E2 = a[i1 + 3];
          phi = a[i1 + 6];
          for (ibtile = 0; ibtile < 3; ibtile++) {
            b_E2[i1 + 3 * ibtile] = (Rz_tmp * b_Rz_tmp[3 * ibtile] +
                                     E2 * b_Rz_tmp[3 * ibtile + 1]) +
                                    phi * b_Rz_tmp[3 * ibtile + 2];
          }
        }
        b_st.site = &s_emlrtRSI;
        coder::internal::blas::mtimes(b_st, b_E2, M, b_r);
        if ((b_r.size(1) != static_cast<int32_T>(spins)) &&
            ((b_r.size(1) != 1) && (static_cast<int32_T>(spins) != 1))) {
          emlrtDimSizeImpxCheckR2021b(b_r.size(1), static_cast<int32_T>(spins),
                                      &c_emlrtECI, (emlrtConstCTX)sp);
        }
        if (b_r.size(1) == Bti.size(1)) {
          M.set_size(&s_emlrtRTEI, sp, 3, b_r.size(1));
          loop_ub = 3 * b_r.size(1);
          ibtile = (loop_ub / 2) << 1;
          jtilecol = ibtile - 2;
          for (i1 = 0; i1 <= jtilecol; i1 += 2) {
            r1 = _mm_loadu_pd(&b_r[i1]);
            r2 = _mm_loadu_pd(&Bti[i1]);
            _mm_storeu_pd(&M[i1], _mm_add_pd(r1, r2));
          }
          for (i1 = ibtile; i1 < loop_ub; i1++) {
            M[i1] = b_r[i1] + Bti[i1];
          }
        } else {
          st.site = &j_emlrtRSI;
          plus(st, M, b_r, Bti);
        }
        //  Rudy 240627: spoiler for FISP taken away
        //  for s=1:spins
        //      M(:,s)=Z(:,:,s)*M(:,s);
        //  end
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
      }
      //  timepoints
      //  relaxation time during the waiting time
      //
      // 	Function simulates free precession and decay
      // 	over a time interval T, given relaxation times T1 and T2
      // 	and off-resonance df.  Times in ms, off-resonance in Hz.
      //  Resonant precession, radians.
      Bti.set_size(&t_emlrtRTEI, sp, 3, static_cast<int32_T>(spins));
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        Bti[3 * i1] = 0.0;
        Bti[3 * i1 + 1] = 0.0;
        Bti[3 * i1 + 2] = 1.0 - E1;
      }
      st.site = &k_emlrtRSI;
      a[0] = c_E2;
      a[3] = 0.0;
      a[6] = 0.0;
      a[1] = 0.0;
      a[4] = c_E2;
      a[7] = 0.0;
      a[2] = 0.0;
      a[5] = 0.0;
      a[8] = E1;
      b_Rz_tmp[0] = d_Rz_tmp;
      b_Rz_tmp[3] = -c_Rz_tmp;
      b_Rz_tmp[6] = 0.0;
      b_Rz_tmp[1] = c_Rz_tmp;
      b_Rz_tmp[4] = d_Rz_tmp;
      b_Rz_tmp[7] = 0.0;
      b_Rz_tmp[2] = 0.0;
      b_Rz_tmp[5] = 0.0;
      b_Rz_tmp[8] = 1.0;
      for (i1 = 0; i1 < 3; i1++) {
        Rz_tmp = a[i1];
        E2 = a[i1 + 3];
        phi = a[i1 + 6];
        for (ibtile = 0; ibtile < 3; ibtile++) {
          b_E2[i1 + 3 * ibtile] =
              (Rz_tmp * b_Rz_tmp[3 * ibtile] + E2 * b_Rz_tmp[3 * ibtile + 1]) +
              phi * b_Rz_tmp[3 * ibtile + 2];
        }
      }
      b_st.site = &s_emlrtRSI;
      coder::internal::blas::mtimes(b_st, b_E2, M, b_r);
      if ((b_r.size(1) != static_cast<int32_T>(spins)) &&
          ((b_r.size(1) != 1) && (static_cast<int32_T>(spins) != 1))) {
        emlrtDimSizeImpxCheckR2021b(b_r.size(1), static_cast<int32_T>(spins),
                                    &d_emlrtECI, (emlrtConstCTX)sp);
      }
      if (b_r.size(1) == Bti.size(1)) {
        loop_ub = 3 * b_r.size(1);
        b_r.set_size(&u_emlrtRTEI, sp, 3, b_r.size(1));
        ibtile = (loop_ub / 2) << 1;
        jtilecol = ibtile - 2;
        for (i1 = 0; i1 <= jtilecol; i1 += 2) {
          r1 = _mm_loadu_pd(&b_r[i1]);
          r2 = _mm_loadu_pd(&Bti[i1]);
          _mm_storeu_pd(&b_r[i1], _mm_add_pd(r1, r2));
        }
        for (i1 = ibtile; i1 < loop_ub; i1++) {
          b_r[i1] = b_r[i1] + Bti[i1];
        }
      } else {
        st.site = &k_emlrtRSI;
        plus(st, b_r, Bti);
      }
      iv[0] = 3;
      loop_ub = M0.size(1);
      iv[1] = M0.size(1);
      emlrtSubAssignSizeCheckR2012b(&iv[0], 2, b_r.size(), 2, &e_emlrtECI,
                                    (emlrtCTX)sp);
      for (i1 = 0; i1 < loop_ub; i1++) {
        M0[3 * i1 + 3 * M0.size(1) * (iEcho + 1)] = b_r[3 * i1];
        M0[(3 * i1 + 3 * M0.size(1) * (iEcho + 1)) + 1] = b_r[3 * i1 + 1];
        M0[(3 * i1 + 3 * M0.size(1) * (iEcho + 1)) + 2] = b_r[3 * i1 + 2];
      }
      //  disp('done1')
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b((emlrtConstCTX)sp);
      }
    }
    // end
    //  disp('done2')
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  //  T1 T2 entry
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

// End of code generation (dictionary_TRUEFISP_withRelaxDecay.cpp)
