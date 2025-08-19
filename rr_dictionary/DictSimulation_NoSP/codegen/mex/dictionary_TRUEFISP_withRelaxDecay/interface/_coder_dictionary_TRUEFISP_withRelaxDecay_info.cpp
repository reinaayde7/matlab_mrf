//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// _coder_dictionary_TRUEFISP_withRelaxDecay_info.cpp
//
// Code generation for function 'dictionary_TRUEFISP_withRelaxDecay'
//

// Include files
#include "_coder_dictionary_TRUEFISP_withRelaxDecay_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

// Function Declarations
static const mxArray *c_emlrtMexFcnResolvedFunctionsI();

// Function Definitions
static const mxArray *c_emlrtMexFcnResolvedFunctionsI()
{
  const mxArray *nameCaptureInfo;
  const char_T *data[6]{
      "789ced585f6f935014bf35d59818b54fbef9115cd6e196cdc407374a9db394f1674b2686"
      "51b8b3745ca870cb60dfc027bf86cf7e0a1f8df1c1f829fc088e42a1"
      "9090626850db7b1e7a38f757ceef9e73e0177241e3b0d700003c00917d7a1cf9fb71dc8a"
      "fd2d90b53cde88fded5c0c92f566e6be19fe31f69a6d61e8e328b054",
      "04933b751b19966a61311843e040d7363da84f910bc384a281a0301fb061849839280942"
      "28bc3e1842ed529820e00cdd7487e67c90f4e35741bdcd92fdf850d0"
      "8f560e7fd379db7926737d41a4fb078a64f7e4cd5d8573ec11d4b0d2e31999a243a73087"
      "0217aebb58c5f009520d4bd1b79e6e6e536dd97114ddd0b0615baa13",
      "c8f4cda560a089a9862b0a6b0b9c9cc28ac84b9d69ae2b030f7968aa3e0d3535d840b9fa"
      "cf2bd65ff43cb462040f1d1b2f91ef4e215f84e8f66460c294ef6b45"
      "3eb3902f8bd732df6933a3112eeae7c392f5e57dfaffbb53fffdf3cf469d7ce0f9976fb5"
      "f2c5f6b7f8fc827c659fcf47057cad1cce7a9470a6be989c22f6981f",
      "89ef025642fb4cba0f6e01cfa27d8082b8aefcebfa9e97adaf998bd3fa22e43a92e9ff56"
      "a72f0bf9b2782df3bb4ec647747a597cb1adba4e77c7fc0ec3ddfc6c"
      "31fa1eb5afd303eb64bb4b749ae87484f844a797373f9fe834d1e96cfeb23add79e579fd"
      "634ad2476df63da67df1f42468d3aba3d3eb7e2e32ae58ffbd05f5cf",
      "f00b07c2b10335e8ba19fef38afc7faaeb3f2af239857c59bc9679cf35753658a2ef4be2"
      "8b6dd5f57dc730a5a3eecbd794676994b53bd0ba677bfe0a7d87affb"
      "fb4ece4d22fba7bec7c9b909d1eb5cfe753f37f90df83849bf",
      ""};
  nameCaptureInfo = nullptr;
  emlrtNameCaptureMxArrayR2016a(&data[0], 7616U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties()
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[9]{"Version",
                                 "ResolvedFunctions",
                                 "Checksum",
                                 "EntryPoints",
                                 "CoverageInfo",
                                 "IsPolymorphic",
                                 "PropertyList",
                                 "UUID",
                                 "ClassEntryPointIsHandle"};
  const char_T *epFieldName[8]{
      "Name",     "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "FullPath", "TimeStamp",      "Constructor",     "Visible"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 1, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 10);
  emlrtSetField(xEntryPoints, 0, "Name",
                emlrtMxCreateString("dictionary_TRUEFISP_withRelaxDecay"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(10.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(xEntryPoints, 0, "FullPath",
                emlrtMxCreateString(
                    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-"
                    "main_d240531\\rr_dictionary\\DictSimulation_"
                    "NoSP\\dictionary_TRUEFISP_with"
                    "RelaxDecay.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739449.46298611106));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.1.0.2537033 (R2024a)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("YNZPhgwFsMfTikcfhx4D0B"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

// End of code generation (_coder_dictionary_TRUEFISP_withRelaxDecay_info.cpp)
