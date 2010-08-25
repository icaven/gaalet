#ifndef __GAALET_H
#define __GAALET_H

#ifdef __CUDACC__
   #define GAALET_CUDA_HOST_DEVICE __host__ __device__
#else
   #define GAALET_CUDA_HOST_DEVICE
#endif

#include "multivector.h"
#include "algebra.h"

#include "addition.h"
#include "multiplication.h"
#include "reverse.h"
#include "grade.h"
#include "part.h"
#include "inverse.h"
#include "scalar.h"

#endif
