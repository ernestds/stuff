/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_untitled_api.h
 *
 * Code generation for function '_coder_untitled_api'
 *
 */

#ifndef _CODER_UNTITLED_API_H
#define _CODER_UNTITLED_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_untitled_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void untitled(int32_T inputArg1, int32_T inputArg2, real_T *outputArg1,
                     real_T *outputArg2);
extern void untitled_api(const mxArray * const prhs[2], int32_T nlhs, const
  mxArray *plhs[2]);
extern void untitled_atexit(void);
extern void untitled_initialize(void);
extern void untitled_terminate(void);
extern void untitled_xil_terminate(void);

#endif

/* End of code generation (_coder_untitled_api.h) */
