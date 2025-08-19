//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// _coder_dictionary_FISP_withRelaxDecay_info.cpp
//
// Code generation for function 'dictionary_FISP_withRelaxDecay'
//

// Include files
#include "_coder_dictionary_FISP_withRelaxDecay_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

// Function Declarations
static const mxArray *c_emlrtMexFcnResolvedFunctionsI();

// Function Definitions
static const mxArray *c_emlrtMexFcnResolvedFunctionsI()
{
  const mxArray *nameCaptureInfo;
  const char_T *data[6]{
      "789ced59cd6ed34010dea014212120275e8033a8a969d52271a0751c0ac47163a745c5c8"
      "38f6b671ea9fe06c8ccd9d33efc25370448803e209780c6aaf13ff48"
      "8b23125c70760e198f3f7b3ecf8cfd69b501b5c34e0d00701b607b790ffb5b71dc88fd35"
      "90b53c5e8bfdf55c3cb30d50cfdc37c33fc65e736c047d84035bb5e0",
      "fc4eddb10c5bb591148c2170e1c4313da847c8996142c9b0a0980ef830b2b814340f4228"
      "3c3e1842ed429c5ac01d4e922734d3c1bc1f3f09f5d617ecc707423f"
      "1a39fc55eb75eb912c744589ed1e287da7236fee2a82eb8ca086944e8f931936740a7728"
      "0ae1f9095211bc6fa986ade85b0f37b799a6ecba8a6e68c8706cd50d",
      "64f6f25034aca9a9866714de11053981719e77061af6a0a9fa2cd4d4e08195aafbcd9275"
      "d77371523746debb4e34ec55f191debb468ce8ce7460c284efcb927c"
      "1744be2c5eca5cc36646d32beee79d05ebcbfbe4fa1b91fff6e947ad4c3ef0f8f3d752f9"
      "62bb2a3e9f906fd1f7f32e81af91c3dbe3de0e275cfe6c71fa1eb3af",
      "b303fb78bb9d3c8750c053f41c801097959feaf79fd5bd41ac1b23688805fc7fd56f93c8"
      "97c54b996bd44c2ce054bf57c4175bd5f59bf718f1547d323db1f8a3"
      "de483a0ff8beb5cf5547bfd7f53ba7eb6c6c749d4d75fa6ff2d175f66af2539dfe7d7d45"
      "3aed539d5eddfc7caad354a7b3f917d5e9d633cfeb1e317d7dd4e4df",
      "22d6974e8e83265b1d9d5ed7fd90f19275df2ca87b869fb9108e5da8c109fe37e1aaf4fc"
      "fb927c2e912f8b9732e754536359a7babe2abed8aaaeeb3b86d97fde"
      "7efa82f16c8db177075afb74cfafd0fa7bddbf77ba5f82ed9f5a87d3fd12aad7b9fcebbe"
      "5ff20bc696b782",
      ""};
  nameCaptureInfo = nullptr;
  emlrtNameCaptureMxArrayR2016a(&data[0], 9072U, &nameCaptureInfo);
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
  xInputs = emlrtCreateLogicalMatrix(1, 11);
  emlrtSetField(xEntryPoints, 0, "Name",
                emlrtMxCreateString("dictionary_FISP_withRelaxDecay"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(11.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(xEntryPoints, 0, "FullPath",
                emlrtMxCreateString(
                    "E:\\POSTDOC_UoM\\08_Project_MRF\\3D_MRF_FISP_Prostate-"
                    "main_d240531\\rr_dictionary\\DictSimulation_"
                    "NoSP\\dictionary_FISP_withRela"
                    "xDecay.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739451.38675925927));
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
                emlrtMxCreateString("Z5ZKVG1cknrxI1be6prCIC"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

// End of code generation (_coder_dictionary_FISP_withRelaxDecay_info.cpp)
