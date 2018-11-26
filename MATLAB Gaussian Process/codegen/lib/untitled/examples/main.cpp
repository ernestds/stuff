/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.cpp
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "untitled.h"
#include "main.h"
#include "untitled_terminate.h"
#include "untitled_initialize.h"

/* Function Declarations */
static int argInit_int32_T();
static void main_untitled();

/* Function Definitions */
static int argInit_int32_T()
{
  return 0;
}

static void main_untitled()
{
  double outputArg1;
  double outputArg2;

  /* Initialize function 'untitled' input arguments. */
  /* Call the entry-point 'untitled'. */
  untitled(argInit_int32_T(), argInit_int32_T(), &outputArg1, &outputArg2);
}

int main(int, const char * const [])
{
  /* Initialize the application.
     You do not need to do this more than one time. */
  untitled_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_untitled();

  /* Terminate the application.
     You do not need to do this more than one time. */
  untitled_terminate();
  return 0;
}

/* End of code generation (main.cpp) */
