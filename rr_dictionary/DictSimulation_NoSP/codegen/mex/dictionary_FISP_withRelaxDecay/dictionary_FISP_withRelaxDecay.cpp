//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// dictionary_FISP_withRelaxDecay.cpp
//
// Code generation for function 'dictionary_FISP_withRelaxDecay'
//

// Include files
#include "dictionary_FISP_withRelaxDecay.h"
#include "assertValidSizeArg.h"
#include "dictionary_FISP_withRelaxDecay_data.h"
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
    18,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo b_emlrtRSI{
    33,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo c_emlrtRSI{
    34,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo d_emlrtRSI{
    42,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo e_emlrtRSI{
    46,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo f_emlrtRSI{
    58,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo g_emlrtRSI{
    61,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo h_emlrtRSI{
    64,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo i_emlrtRSI{
    65,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo j_emlrtRSI{
    66,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo k_emlrtRSI{
    72,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo l_emlrtRSI{
    84,                               // lineNo
    "dictionary_FISP_withRelaxDecay", // fcnName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pathName
};

static emlrtRSInfo m_emlrtRSI{
    38,       // lineNo
    "repmat", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\elmat\\repmat.m" // pathName
};

static emlrtRSInfo n_emlrtRSI{
    74,       // lineNo
    "repmat", // fcnName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\elmat\\repmat.m" // pathName
};

static emlrtRSInfo
    t_emlrtRSI{
        94,                  // lineNo
        "eml_mtimes_helper", // fcnName
        "C:\\Program "
        "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\ops\\eml_mtimes_"
        "helper.m" // pathName
    };

static emlrtRTEInfo emlrtRTEI{
    20,                               // lineNo
    9,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtBCInfo emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    21,                               // lineNo
    30,                               // colNo
    "phirange",                       // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtDCInfo emlrtDCI{
    23,                               // lineNo
    17,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    1           // checkKind
};

static emlrtECInfo emlrtECI{
    2,                                // nDims
    46,                               // lineNo
    21,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo b_emlrtRTEI{
    49,                               // lineNo
    23,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtBCInfo b_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    50,                               // lineNo
    34,                               // colNo
    "Z",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo c_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    50,                               // lineNo
    41,                               // colNo
    "M",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtECInfo b_emlrtECI{
    2,                                // nDims
    61,                               // lineNo
    17,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtECInfo c_emlrtECI{
    2,                                // nDims
    72,                               // lineNo
    15,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo c_emlrtRTEI{
    73,                               // lineNo
    19,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtBCInfo d_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    74,                               // lineNo
    30,                               // colNo
    "Z",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo e_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    74,                               // lineNo
    37,                               // colNo
    "M",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtECInfo d_emlrtECI{
    2,                                // nDims
    84,                               // lineNo
    12,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo d_emlrtRTEI{
    33,         // lineNo
    37,         // colNo
    "linspace", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\elmat\\linspace.m" // pName
};

static emlrtBCInfo f_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    21,                               // lineNo
    11,                               // colNo
    "Z",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo g_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    50,                               // lineNo
    25,                               // colNo
    "M",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo h_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    74,                               // lineNo
    21,                               // colNo
    "M",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtDCInfo b_emlrtDCI{
    14,                               // lineNo
    12,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    1           // checkKind
};

static emlrtDCInfo c_emlrtDCI{
    14,                               // lineNo
    12,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    4           // checkKind
};

static emlrtDCInfo d_emlrtDCI{
    14,                               // lineNo
    1,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    1           // checkKind
};

static emlrtDCInfo e_emlrtDCI{
    19,                               // lineNo
    15,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    1           // checkKind
};

static emlrtDCInfo f_emlrtDCI{
    19,                               // lineNo
    15,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    4           // checkKind
};

static emlrtBCInfo i_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    29,                               // lineNo
    12,                               // colNo
    "r",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo j_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    30,                               // lineNo
    12,                               // colNo
    "r",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo k_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    31,                               // lineNo
    12,                               // colNo
    "r",                              // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo l_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    58,                               // lineNo
    28,                               // colNo
    "flip",                           // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo m_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    58,                               // lineNo
    37,                               // colNo
    "phase",                          // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo n_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    69,                               // lineNo
    34,                               // colNo
    "tr",                             // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo o_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    64,                               // lineNo
    20,                               // colNo
    "Mx",                             // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo p_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    64,                               // lineNo
    22,                               // colNo
    "Mx",                             // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo q_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    65,                               // lineNo
    20,                               // colNo
    "My",                             // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo r_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    65,                               // lineNo
    22,                               // colNo
    "My",                             // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo s_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    66,                               // lineNo
    20,                               // colNo
    "Mz",                             // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtBCInfo t_emlrtBCI{
    -1,                               // iFirst
    -1,                               // iLast
    66,                               // lineNo
    22,                               // colNo
    "Mz",                             // aName
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m", // pName
    0           // checkKind
};

static emlrtRTEInfo g_emlrtRTEI{
    14,                               // lineNo
    1,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo h_emlrtRTEI{
    15,                               // lineNo
    1,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo i_emlrtRTEI{
    16,                               // lineNo
    1,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo j_emlrtRTEI{
    49,         // lineNo
    20,         // colNo
    "linspace", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\elmat\\linspace.m" // pName
};

static emlrtRTEInfo k_emlrtRTEI{
    18,                               // lineNo
    1,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo l_emlrtRTEI{
    19,                               // lineNo
    5,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo m_emlrtRTEI{
    69,       // lineNo
    28,       // colNo
    "repmat", // fName
    "C:\\Program "
    "Files\\MATLAB\\R2024a\\toolbox\\eml\\lib\\matlab\\elmat\\repmat.m" // pName
};

static emlrtRTEInfo n_emlrtRTEI{
    45,                               // lineNo
    17,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo o_emlrtRTEI{
    46,                               // lineNo
    17,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo p_emlrtRTEI{
    58,                               // lineNo
    41,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo q_emlrtRTEI{
    60,                               // lineNo
    13,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo r_emlrtRTEI{
    61,                               // lineNo
    13,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo s_emlrtRTEI{
    64,                               // lineNo
    32,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo t_emlrtRTEI{
    65,                               // lineNo
    32,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo u_emlrtRTEI{
    66,                               // lineNo
    32,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo v_emlrtRTEI{
    71,                               // lineNo
    13,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo w_emlrtRTEI{
    72,                               // lineNo
    13,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo x_emlrtRTEI{
    83,                               // lineNo
    9,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo y_emlrtRTEI{
    84,                               // lineNo
    9,                                // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
};

static emlrtRTEInfo cb_emlrtRTEI{
    84,                               // lineNo
    12,                               // colNo
    "dictionary_FISP_withRelaxDecay", // fName
    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-main_d240531\\rr_"
    "dictionary\\DictSimulation_NoSP\\dictionary_FISP_withRela"
    "xDecay.m" // pName
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
  in1.set_size(&w_emlrtRTEI, &sp, 3, in1.size(1));
  if (in3.size(1) == 1) {
    loop_ub = in2.size(1);
  } else {
    loop_ub = in3.size(1);
  }
  in1.set_size(&w_emlrtRTEI, &sp, in1.size(0), loop_ub);
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
  b_in1.set_size(&cb_emlrtRTEI, &sp, 3, loop_ub);
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
  in1.set_size(&cb_emlrtRTEI, &sp, 3, in1.size(1));
  loop_ub = b_in1.size(1);
  in1.set_size(&cb_emlrtRTEI, &sp, in1.size(0), b_in1.size(1));
  for (int32_T i{0}; i < loop_ub; i++) {
    in1[3 * i] = b_in1[3 * i];
    in1[3 * i + 1] = b_in1[3 * i + 1];
    in1[3 * i + 2] = b_in1[3 * i + 2];
  }
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)&sp);
}

void dictionary_FISP_withRelaxDecay(
    const emlrtStack *sp, const coder::array<real_T, 2U> &r,
    const coder::array<real_T, 1U> &flip, const coder::array<real_T, 1U> &tr,
    real_T te, const coder::array<real_T, 1U> &phase, real_T ti,
    real_T phasetwist, real_T spins, real_T frames, real_T Inv,
    real_T waitingDuration, coder::array<real32_T, 2U> &Mx,
    coder::array<real32_T, 2U> &My, coder::array<real32_T, 2U> &Mz)
{
  coder::array<real_T, 3U> Z;
  coder::array<real_T, 2U> Bti;
  coder::array<real_T, 2U> M;
  coder::array<real_T, 2U> M0;
  coder::array<real_T, 2U> b_r;
  coder::array<real_T, 2U> phirange;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  real_T E1;
  real_T E2;
  real_T Rz_tmp;
  real_T b_Rz_tmp;
  real_T delta1;
  real_T delta2;
  real_T nt;
  int32_T b_ntilecols;
  int32_T b_outsize_idx_1;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T ibtile;
  int32_T loop_ub;
  int32_T ntilecols;
  int32_T outsize_idx_1;
  boolean_T b;
  boolean_T b_overflow;
  boolean_T overflow;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  nt = flip.size(0);
  if (flip.size(0) > frames) {
    nt = frames;
  }
  // df=0;
  if (!(nt >= 0.0)) {
    emlrtNonNegativeCheckR2012b(nt, &c_emlrtDCI, (emlrtConstCTX)sp);
  }
  E1 = muDoubleScalarFloor(nt);
  if (nt != E1) {
    emlrtIntegerCheckR2012b(nt, &b_emlrtDCI, (emlrtConstCTX)sp);
  }
  i = r.size(0);
  Mx.set_size(&g_emlrtRTEI, sp, static_cast<int32_T>(nt), r.size(0));
  if (nt != E1) {
    emlrtIntegerCheckR2012b(nt, &d_emlrtDCI, (emlrtConstCTX)sp);
  }
  ibtile = static_cast<int32_T>(nt) * r.size(0);
  for (i1 = 0; i1 < ibtile; i1++) {
    Mx[i1] = 0.0F;
  }
  My.set_size(&h_emlrtRTEI, sp, static_cast<int32_T>(nt), r.size(0));
  for (i1 = 0; i1 < ibtile; i1++) {
    My[i1] = 0.0F;
  }
  Mz.set_size(&i_emlrtRTEI, sp, static_cast<int32_T>(nt), r.size(0));
  for (i1 = 0; i1 < ibtile; i1++) {
    Mz[i1] = 0.0F;
  }
  st.site = &emlrtRSI;
  E2 = -phasetwist / 2.0;
  delta2 = phasetwist / 2.0;
  b = !(spins >= 0.0);
  if (b) {
    if (muDoubleScalarIsNaN(spins)) {
      emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
                                    "Coder:toolbox:MustNotBeNaN",
                                    "Coder:toolbox:MustNotBeNaN", 3, 4, 1, "N");
    }
    phirange.set_size(&k_emlrtRTEI, &st, 1, 0);
  } else {
    E1 = muDoubleScalarFloor(spins);
    i1 = static_cast<int32_T>(E1);
    phirange.set_size(&j_emlrtRTEI, &st, 1, i1);
    if (static_cast<int32_T>(E1) >= 1) {
      ibtile = static_cast<int32_T>(E1) - 1;
      phirange[static_cast<int32_T>(E1) - 1] = delta2;
      if (phirange.size(1) >= 2) {
        phirange[0] = E2;
        if (phirange.size(1) >= 3) {
          if (E2 == -delta2) {
            delta2 /= static_cast<real_T>(phirange.size(1)) - 1.0;
            for (int32_T k{2}; k <= ibtile; k++) {
              phirange[k - 1] =
                  static_cast<real_T>(((k << 1) - phirange.size(1)) - 1) *
                  delta2;
            }
            if ((phirange.size(1) & 1) == 1) {
              phirange[phirange.size(1) >> 1] = 0.0;
            }
          } else if (((E2 < 0.0) != (delta2 < 0.0)) &&
                     ((muDoubleScalarAbs(E2) > 8.9884656743115785E+307) ||
                      (muDoubleScalarAbs(delta2) > 8.9884656743115785E+307))) {
            delta1 = E2 / (static_cast<real_T>(phirange.size(1)) - 1.0);
            delta2 /= static_cast<real_T>(phirange.size(1)) - 1.0;
            for (int32_T k{0}; k <= i1 - 3; k++) {
              phirange[k + 1] = (E2 + delta2 * (static_cast<real_T>(k) + 1.0)) -
                                delta1 * (static_cast<real_T>(k) + 1.0);
            }
          } else {
            delta1 =
                (delta2 - E2) / (static_cast<real_T>(phirange.size(1)) - 1.0);
            for (int32_T k{0}; k <= i1 - 3; k++) {
              phirange[k + 1] = E2 + (static_cast<real_T>(k) + 1.0) * delta1;
            }
          }
        }
      }
    }
  }
  Z.set_size(&l_emlrtRTEI, sp, 3, 3, Z.size(2));
  if (b) {
    emlrtNonNegativeCheckR2012b(spins, &f_emlrtDCI, (emlrtConstCTX)sp);
  }
  E1 = static_cast<int32_T>(muDoubleScalarFloor(spins));
  if (spins != E1) {
    emlrtIntegerCheckR2012b(spins, &e_emlrtDCI, (emlrtConstCTX)sp);
  }
  Z.set_size(&l_emlrtRTEI, sp, Z.size(0), Z.size(1),
             static_cast<int32_T>(spins));
  loop_ub = static_cast<int32_T>(spins);
  emlrtForLoopVectorCheckR2021a(1.0, 1.0, spins, mxDOUBLE_CLASS,
                                static_cast<int32_T>(spins), &emlrtRTEI,
                                (emlrtConstCTX)sp);
  for (ibtile = 0; ibtile < loop_ub; ibtile++) {
    if ((ibtile + 1 < 1) || (ibtile + 1 > phirange.size(1))) {
      emlrtDynamicBoundsCheckR2012b(ibtile + 1, 1, phirange.size(1), &emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    Rz_tmp = muDoubleScalarSin(phirange[ibtile]);
    b_Rz_tmp = muDoubleScalarCos(phirange[ibtile]);
    if ((static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U) < 1) ||
        (static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U) >
         Z.size(2))) {
      emlrtDynamicBoundsCheckR2012b(
          static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U), 1,
          Z.size(2), &f_emlrtBCI, (emlrtConstCTX)sp);
    }
    Z[9 * ibtile] = b_Rz_tmp;
    Z[9 * ibtile + 3] = -Rz_tmp;
    Z[9 * ibtile + 6] = 0.0;
    Z[9 * ibtile + 1] = Rz_tmp;
    Z[9 * ibtile + 4] = b_Rz_tmp;
    Z[9 * ibtile + 7] = 0.0;
    Z[9 * ibtile + 2] = 0.0;
    Z[9 * ibtile + 5] = 0.0;
    Z[9 * ibtile + 8] = 1.0;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  if (spins != E1) {
    emlrtIntegerCheckR2012b(spins, &emlrtDCI, (emlrtConstCTX)sp);
  }
  if (r.size(0) - 1 >= 0) {
    outsize_idx_1 = static_cast<int32_T>(spins);
    ntilecols = static_cast<int32_T>(spins);
    overflow = (static_cast<int32_T>(spins) > 2147483646);
    b_outsize_idx_1 = static_cast<int32_T>(spins);
    b_ntilecols = static_cast<int32_T>(spins);
    b_overflow = (static_cast<int32_T>(spins) > 2147483646);
    i2 = static_cast<int32_T>(nt);
  }
  for (int32_T n{0}; n < i; n++) {
    real_T T1;
    real_T T2;
    real_T b_E1;
    real_T b_E2;
    real_T df;
    int32_T jtilecol;
    if (n + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(n + 1, 1, i, &i_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    T1 = r[n] / 1000.0;
    if (n + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(n + 1, 1, i, &j_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    T2 = r[n + r.size(0)] / 1000.0;
    if (n + 1 > i) {
      emlrtDynamicBoundsCheckR2012b(n + 1, 1, i, &k_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    df = r[n + r.size(0) * 2];
    st.site = &b_emlrtRSI;
    b_st.site = &m_emlrtRSI;
    coder::internal::assertValidSizeArg(b_st, spins);
    M0.set_size(&m_emlrtRTEI, &st, 3, outsize_idx_1);
    b_st.site = &n_emlrtRSI;
    if (overflow) {
      c_st.site = &o_emlrtRSI;
      coder::check_forloop_overflow_error(c_st);
    }
    for (jtilecol = 0; jtilecol < ntilecols; jtilecol++) {
      ibtile = jtilecol * 3;
      M0[ibtile] = 0.0;
      M0[ibtile + 1] = 0.0;
      M0[ibtile + 2] = 1.0;
    }
    st.site = &c_emlrtRSI;
    b_st.site = &m_emlrtRSI;
    coder::internal::assertValidSizeArg(b_st, spins);
    M.set_size(&m_emlrtRTEI, &st, 3, b_outsize_idx_1);
    b_st.site = &n_emlrtRSI;
    if (b_overflow) {
      c_st.site = &o_emlrtRSI;
      coder::check_forloop_overflow_error(c_st);
    }
    for (jtilecol = 0; jtilecol < b_ntilecols; jtilecol++) {
      ibtile = jtilecol * 3;
      M[ibtile] = 0.0;
      M[ibtile + 1] = 0.0;
      M[ibtile + 2] = 0.0;
    }
    delta1 = 6.2831853071795862 * df * waitingDuration;
    b_E1 = muDoubleScalarExp(-waitingDuration / T1);
    b_E2 = muDoubleScalarExp(-waitingDuration / T2);
    Rz_tmp = muDoubleScalarSin(delta1);
    b_Rz_tmp = muDoubleScalarCos(delta1);
    for (int32_T iEcho{0}; iEcho < 2; iEcho++) {
      __m128d r1;
      __m128d r2;
      real_T a[9];
      real_T c_E2[9];
      real_T c_Rz_tmp[9];
      int32_T vectorUB;
      // Rudy: run it twice to get to steady state
      // iEcho
      for (int32_T k{0}; k < i2; k++) {
        real_T b_Z[3];
        if (k + 1 == 1) {
          if (Inv == 1.0) {
            st.site = &d_emlrtRSI;
            b_st.site = &d_emlrtRSI;
            throt(b_st, 3.1415926535897931, 0.0, a);
            b_st.site = &t_emlrtRSI;
            coder::internal::blas::mtimes(b_st, a, M0, M);
          }
          //
          // 	Function simulates free precession and decay
          // 	over a time interval T, given relaxation times T1 and T2
          // 	and off-resonance df.  Times in ms, off-resonance in Hz.
          delta1 = 6.2831853071795862 * df * ti;
          //  Resonant precession, radians.
          E1 = muDoubleScalarExp(-ti / T1);
          E2 = muDoubleScalarExp(-ti / T2);
          nt = muDoubleScalarSin(delta1);
          delta2 = muDoubleScalarCos(delta1);
          Bti.set_size(&n_emlrtRTEI, sp, 3, static_cast<int32_T>(spins));
          for (i1 = 0; i1 < loop_ub; i1++) {
            Bti[3 * i1] = 0.0;
            Bti[3 * i1 + 1] = 0.0;
            Bti[3 * i1 + 2] = 1.0 - E1;
          }
          st.site = &e_emlrtRSI;
          a[0] = E2;
          a[3] = 0.0;
          a[6] = 0.0;
          a[1] = 0.0;
          a[4] = E2;
          a[7] = 0.0;
          a[2] = 0.0;
          a[5] = 0.0;
          a[8] = E1;
          c_Rz_tmp[0] = delta2;
          c_Rz_tmp[3] = -nt;
          c_Rz_tmp[6] = 0.0;
          c_Rz_tmp[1] = nt;
          c_Rz_tmp[4] = delta2;
          c_Rz_tmp[7] = 0.0;
          c_Rz_tmp[2] = 0.0;
          c_Rz_tmp[5] = 0.0;
          c_Rz_tmp[8] = 1.0;
          for (i1 = 0; i1 < 3; i1++) {
            E1 = a[i1];
            delta2 = a[i1 + 3];
            delta1 = a[i1 + 6];
            for (ibtile = 0; ibtile < 3; ibtile++) {
              c_E2[i1 + 3 * ibtile] = (E1 * c_Rz_tmp[3 * ibtile] +
                                       delta2 * c_Rz_tmp[3 * ibtile + 1]) +
                                      delta1 * c_Rz_tmp[3 * ibtile + 2];
            }
          }
          b_st.site = &t_emlrtRSI;
          coder::internal::blas::mtimes(b_st, c_E2, M, b_r);
          if ((b_r.size(1) != static_cast<int32_T>(spins)) &&
              ((b_r.size(1) != 1) && (static_cast<int32_T>(spins) != 1))) {
            emlrtDimSizeImpxCheckR2021b(b_r.size(1),
                                        static_cast<int32_T>(spins), &emlrtECI,
                                        (emlrtConstCTX)sp);
          }
          if (b_r.size(1) == Bti.size(1)) {
            M.set_size(&o_emlrtRTEI, sp, 3, b_r.size(1));
            ibtile = 3 * b_r.size(1);
            jtilecol = (ibtile / 2) << 1;
            vectorUB = jtilecol - 2;
            for (i1 = 0; i1 <= vectorUB; i1 += 2) {
              r1 = _mm_loadu_pd(&b_r[i1]);
              r2 = _mm_loadu_pd(&Bti[i1]);
              _mm_storeu_pd(&M[i1], _mm_add_pd(r1, r2));
            }
            for (i1 = jtilecol; i1 < ibtile; i1++) {
              M[i1] = b_r[i1] + Bti[i1];
            }
          } else {
            st.site = &e_emlrtRSI;
            plus(st, M, b_r, Bti);
          }
          //  spoiling
          emlrtForLoopVectorCheckR2021a(1.0, 1.0, spins, mxDOUBLE_CLASS,
                                        static_cast<int32_T>(spins),
                                        &b_emlrtRTEI, (emlrtConstCTX)sp);
          for (ibtile = 0; ibtile < loop_ub; ibtile++) {
            if ((ibtile + 1 < 1) || (ibtile + 1 > Z.size(2))) {
              emlrtDynamicBoundsCheckR2012b(ibtile + 1, 1, Z.size(2),
                                            &b_emlrtBCI, (emlrtConstCTX)sp);
            }
            if ((ibtile + 1 < 1) || (ibtile + 1 > M.size(1))) {
              emlrtDynamicBoundsCheckR2012b(ibtile + 1, 1, M.size(1),
                                            &c_emlrtBCI, (emlrtConstCTX)sp);
            }
            if ((static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U) <
                 1) ||
                (static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U) >
                 M.size(1))) {
              emlrtDynamicBoundsCheckR2012b(
                  static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U), 1,
                  M.size(1), &g_emlrtBCI, (emlrtConstCTX)sp);
            }
            for (i1 = 0; i1 < 3; i1++) {
              b_Z[i1] = (Z[i1 + 9 * ibtile] * M[3 * ibtile] +
                         Z[(i1 + 9 * ibtile) + 3] * M[3 * ibtile + 1]) +
                        Z[(i1 + 9 * ibtile) + 6] * M[3 * ibtile + 2];
            }
            M[3 * ibtile] = b_Z[0];
            M[3 * ibtile + 1] = b_Z[1];
            M[3 * ibtile + 2] = b_Z[2];
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b((emlrtConstCTX)sp);
            }
          }
        }
        //  alpha pulse
        st.site = &f_emlrtRSI;
        if (k + 1 > flip.size(0)) {
          emlrtDynamicBoundsCheckR2012b(k + 1, 1, flip.size(0), &l_emlrtBCI,
                                        &st);
        }
        if (k + 1 > phase.size(0)) {
          emlrtDynamicBoundsCheckR2012b(k + 1, 1, phase.size(0), &m_emlrtBCI,
                                        &st);
        }
        b_st.site = &f_emlrtRSI;
        throt(b_st, flip[k], phase[k], a);
        Bti.set_size(&p_emlrtRTEI, &st, 3, M.size(1));
        ibtile = M.size(0) * M.size(1) - 1;
        for (i1 = 0; i1 <= ibtile; i1++) {
          Bti[i1] = M[i1];
        }
        b_st.site = &t_emlrtRSI;
        coder::internal::blas::mtimes(b_st, a, Bti, M);
        //
        // 	Function simulates free precession and decay
        // 	over a time interval T, given relaxation times T1 and T2
        // 	and off-resonance df.  Times in ms, off-resonance in Hz.
        delta1 = 6.2831853071795862 * df * te;
        //  Resonant precession, radians.
        E1 = muDoubleScalarExp(-te / T1);
        E2 = muDoubleScalarExp(-te / T2);
        nt = muDoubleScalarSin(delta1);
        delta2 = muDoubleScalarCos(delta1);
        Bti.set_size(&q_emlrtRTEI, sp, 3, static_cast<int32_T>(spins));
        for (i1 = 0; i1 < loop_ub; i1++) {
          Bti[3 * i1] = 0.0;
          Bti[3 * i1 + 1] = 0.0;
          Bti[3 * i1 + 2] = 1.0 - E1;
        }
        st.site = &g_emlrtRSI;
        a[0] = E2;
        a[3] = 0.0;
        a[6] = 0.0;
        a[1] = 0.0;
        a[4] = E2;
        a[7] = 0.0;
        a[2] = 0.0;
        a[5] = 0.0;
        a[8] = E1;
        c_Rz_tmp[0] = delta2;
        c_Rz_tmp[3] = -nt;
        c_Rz_tmp[6] = 0.0;
        c_Rz_tmp[1] = nt;
        c_Rz_tmp[4] = delta2;
        c_Rz_tmp[7] = 0.0;
        c_Rz_tmp[2] = 0.0;
        c_Rz_tmp[5] = 0.0;
        c_Rz_tmp[8] = 1.0;
        for (i1 = 0; i1 < 3; i1++) {
          E1 = a[i1];
          delta2 = a[i1 + 3];
          delta1 = a[i1 + 6];
          for (ibtile = 0; ibtile < 3; ibtile++) {
            c_E2[i1 + 3 * ibtile] = (E1 * c_Rz_tmp[3 * ibtile] +
                                     delta2 * c_Rz_tmp[3 * ibtile + 1]) +
                                    delta1 * c_Rz_tmp[3 * ibtile + 2];
          }
        }
        b_st.site = &t_emlrtRSI;
        coder::internal::blas::mtimes(b_st, c_E2, M, b_r);
        if ((b_r.size(1) != static_cast<int32_T>(spins)) &&
            ((b_r.size(1) != 1) && (static_cast<int32_T>(spins) != 1))) {
          emlrtDimSizeImpxCheckR2021b(b_r.size(1), static_cast<int32_T>(spins),
                                      &b_emlrtECI, (emlrtConstCTX)sp);
        }
        if (b_r.size(1) == Bti.size(1)) {
          M.set_size(&r_emlrtRTEI, sp, 3, b_r.size(1));
          ibtile = 3 * b_r.size(1);
          jtilecol = (ibtile / 2) << 1;
          vectorUB = jtilecol - 2;
          for (i1 = 0; i1 <= vectorUB; i1 += 2) {
            r1 = _mm_loadu_pd(&b_r[i1]);
            r2 = _mm_loadu_pd(&Bti[i1]);
            _mm_storeu_pd(&M[i1], _mm_add_pd(r1, r2));
          }
          for (i1 = jtilecol; i1 < ibtile; i1++) {
            M[i1] = b_r[i1] + Bti[i1];
          }
        } else {
          st.site = &g_emlrtRSI;
          plus(st, M, b_r, Bti);
        }
        if (iEcho + 1 == 2) {
          // Rudy: run it twice to get to steady state
          ibtile = M.size(1);
          phirange.set_size(&s_emlrtRTEI, sp, 1, M.size(1));
          for (i1 = 0; i1 < ibtile; i1++) {
            phirange[i1] = M[3 * i1];
          }
          if (k + 1 > Mx.size(0)) {
            emlrtDynamicBoundsCheckR2012b(k + 1, 1, Mx.size(0), &o_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if (n + 1 > Mx.size(1)) {
            emlrtDynamicBoundsCheckR2012b(n + 1, 1, Mx.size(1), &p_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          st.site = &h_emlrtRSI;
          Mx[k + Mx.size(0) * n] =
              static_cast<real32_T>(coder::mean(st, phirange));
          phirange.set_size(&t_emlrtRTEI, sp, 1, M.size(1));
          for (i1 = 0; i1 < ibtile; i1++) {
            phirange[i1] = M[3 * i1 + 1];
          }
          if (k + 1 > My.size(0)) {
            emlrtDynamicBoundsCheckR2012b(k + 1, 1, My.size(0), &q_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if (n + 1 > My.size(1)) {
            emlrtDynamicBoundsCheckR2012b(n + 1, 1, My.size(1), &r_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          st.site = &i_emlrtRSI;
          My[k + My.size(0) * n] =
              static_cast<real32_T>(coder::mean(st, phirange));
          phirange.set_size(&u_emlrtRTEI, sp, 1, M.size(1));
          for (i1 = 0; i1 < ibtile; i1++) {
            phirange[i1] = M[3 * i1 + 2];
          }
          if (k + 1 > Mz.size(0)) {
            emlrtDynamicBoundsCheckR2012b(k + 1, 1, Mz.size(0), &s_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if (n + 1 > Mz.size(1)) {
            emlrtDynamicBoundsCheckR2012b(n + 1, 1, Mz.size(1), &t_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          st.site = &j_emlrtRSI;
          Mz[k + Mz.size(0) * n] =
              static_cast<real32_T>(coder::mean(st, phirange));
        }
        if (k + 1 > tr.size(0)) {
          emlrtDynamicBoundsCheckR2012b(k + 1, 1, tr.size(0), &n_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        delta2 = tr[k] - te;
        //
        // 	Function simulates free precession and decay
        // 	over a time interval T, given relaxation times T1 and T2
        // 	and off-resonance df.  Times in ms, off-resonance in Hz.
        delta1 = 6.2831853071795862 * df * delta2;
        //  Resonant precession, radians.
        E1 = muDoubleScalarExp(-delta2 / T1);
        E2 = muDoubleScalarExp(-delta2 / T2);
        nt = muDoubleScalarSin(delta1);
        delta2 = muDoubleScalarCos(delta1);
        Bti.set_size(&v_emlrtRTEI, sp, 3, static_cast<int32_T>(spins));
        for (i1 = 0; i1 < loop_ub; i1++) {
          Bti[3 * i1] = 0.0;
          Bti[3 * i1 + 1] = 0.0;
          Bti[3 * i1 + 2] = 1.0 - E1;
        }
        st.site = &k_emlrtRSI;
        a[0] = E2;
        a[3] = 0.0;
        a[6] = 0.0;
        a[1] = 0.0;
        a[4] = E2;
        a[7] = 0.0;
        a[2] = 0.0;
        a[5] = 0.0;
        a[8] = E1;
        c_Rz_tmp[0] = delta2;
        c_Rz_tmp[3] = -nt;
        c_Rz_tmp[6] = 0.0;
        c_Rz_tmp[1] = nt;
        c_Rz_tmp[4] = delta2;
        c_Rz_tmp[7] = 0.0;
        c_Rz_tmp[2] = 0.0;
        c_Rz_tmp[5] = 0.0;
        c_Rz_tmp[8] = 1.0;
        for (i1 = 0; i1 < 3; i1++) {
          E1 = a[i1];
          delta2 = a[i1 + 3];
          delta1 = a[i1 + 6];
          for (ibtile = 0; ibtile < 3; ibtile++) {
            c_E2[i1 + 3 * ibtile] = (E1 * c_Rz_tmp[3 * ibtile] +
                                     delta2 * c_Rz_tmp[3 * ibtile + 1]) +
                                    delta1 * c_Rz_tmp[3 * ibtile + 2];
          }
        }
        b_st.site = &t_emlrtRSI;
        coder::internal::blas::mtimes(b_st, c_E2, M, b_r);
        if ((b_r.size(1) != static_cast<int32_T>(spins)) &&
            ((b_r.size(1) != 1) && (static_cast<int32_T>(spins) != 1))) {
          emlrtDimSizeImpxCheckR2021b(b_r.size(1), static_cast<int32_T>(spins),
                                      &c_emlrtECI, (emlrtConstCTX)sp);
        }
        if (b_r.size(1) == Bti.size(1)) {
          M.set_size(&w_emlrtRTEI, sp, 3, b_r.size(1));
          ibtile = 3 * b_r.size(1);
          jtilecol = (ibtile / 2) << 1;
          vectorUB = jtilecol - 2;
          for (i1 = 0; i1 <= vectorUB; i1 += 2) {
            r1 = _mm_loadu_pd(&b_r[i1]);
            r2 = _mm_loadu_pd(&Bti[i1]);
            _mm_storeu_pd(&M[i1], _mm_add_pd(r1, r2));
          }
          for (i1 = jtilecol; i1 < ibtile; i1++) {
            M[i1] = b_r[i1] + Bti[i1];
          }
        } else {
          st.site = &k_emlrtRSI;
          plus(st, M, b_r, Bti);
        }
        emlrtForLoopVectorCheckR2021a(1.0, 1.0, spins, mxDOUBLE_CLASS,
                                      static_cast<int32_T>(spins), &c_emlrtRTEI,
                                      (emlrtConstCTX)sp);
        for (ibtile = 0; ibtile < loop_ub; ibtile++) {
          if ((ibtile + 1 < 1) || (ibtile + 1 > Z.size(2))) {
            emlrtDynamicBoundsCheckR2012b(ibtile + 1, 1, Z.size(2), &d_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if ((ibtile + 1 < 1) || (ibtile + 1 > M.size(1))) {
            emlrtDynamicBoundsCheckR2012b(ibtile + 1, 1, M.size(1), &e_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          if ((static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U) < 1) ||
              (static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U) >
               M.size(1))) {
            emlrtDynamicBoundsCheckR2012b(
                static_cast<int32_T>(static_cast<uint32_T>(ibtile) + 1U), 1,
                M.size(1), &h_emlrtBCI, (emlrtConstCTX)sp);
          }
          for (i1 = 0; i1 < 3; i1++) {
            b_Z[i1] = (Z[i1 + 9 * ibtile] * M[3 * ibtile] +
                       Z[(i1 + 9 * ibtile) + 3] * M[3 * ibtile + 1]) +
                      Z[(i1 + 9 * ibtile) + 6] * M[3 * ibtile + 2];
          }
          M[3 * ibtile] = b_Z[0];
          M[3 * ibtile + 1] = b_Z[1];
          M[3 * ibtile + 2] = b_Z[2];
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b((emlrtConstCTX)sp);
          }
        }
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
      Bti.set_size(&x_emlrtRTEI, sp, 3, static_cast<int32_T>(spins));
      for (i1 = 0; i1 < loop_ub; i1++) {
        Bti[3 * i1] = 0.0;
        Bti[3 * i1 + 1] = 0.0;
        Bti[3 * i1 + 2] = 1.0 - b_E1;
      }
      st.site = &l_emlrtRSI;
      a[0] = b_E2;
      a[3] = 0.0;
      a[6] = 0.0;
      a[1] = 0.0;
      a[4] = b_E2;
      a[7] = 0.0;
      a[2] = 0.0;
      a[5] = 0.0;
      a[8] = b_E1;
      c_Rz_tmp[0] = b_Rz_tmp;
      c_Rz_tmp[3] = -Rz_tmp;
      c_Rz_tmp[6] = 0.0;
      c_Rz_tmp[1] = Rz_tmp;
      c_Rz_tmp[4] = b_Rz_tmp;
      c_Rz_tmp[7] = 0.0;
      c_Rz_tmp[2] = 0.0;
      c_Rz_tmp[5] = 0.0;
      c_Rz_tmp[8] = 1.0;
      for (i1 = 0; i1 < 3; i1++) {
        E1 = a[i1];
        delta2 = a[i1 + 3];
        delta1 = a[i1 + 6];
        for (ibtile = 0; ibtile < 3; ibtile++) {
          c_E2[i1 + 3 * ibtile] =
              (E1 * c_Rz_tmp[3 * ibtile] + delta2 * c_Rz_tmp[3 * ibtile + 1]) +
              delta1 * c_Rz_tmp[3 * ibtile + 2];
        }
      }
      b_st.site = &t_emlrtRSI;
      coder::internal::blas::mtimes(b_st, c_E2, M, M0);
      if ((M0.size(1) != static_cast<int32_T>(spins)) &&
          ((M0.size(1) != 1) && (static_cast<int32_T>(spins) != 1))) {
        emlrtDimSizeImpxCheckR2021b(M0.size(1), static_cast<int32_T>(spins),
                                    &d_emlrtECI, (emlrtConstCTX)sp);
      }
      if (M0.size(1) == Bti.size(1)) {
        ibtile = 3 * M0.size(1);
        M0.set_size(&y_emlrtRTEI, sp, 3, M0.size(1));
        jtilecol = (ibtile / 2) << 1;
        vectorUB = jtilecol - 2;
        for (i1 = 0; i1 <= vectorUB; i1 += 2) {
          r1 = _mm_loadu_pd(&M0[i1]);
          r2 = _mm_loadu_pd(&Bti[i1]);
          _mm_storeu_pd(&M0[i1], _mm_add_pd(r1, r2));
        }
        for (i1 = jtilecol; i1 < ibtile; i1++) {
          M0[i1] = M0[i1] + Bti[i1];
        }
      } else {
        st.site = &l_emlrtRSI;
        plus(st, M0, Bti);
      }
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b((emlrtConstCTX)sp);
      }
    }
    // end
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  //  T1 T2 entry
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

// End of code generation (dictionary_FISP_withRelaxDecay.cpp)
