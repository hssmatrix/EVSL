#ifndef STRUCT_H
#define STRUCT_H

#include <complex.h>

/**
 * @brief sparse matrix format: the coordinate (COO) format, 0-based
 *
 * ir, jc, vv : triples for all nonzeros (of size nnz)
 */ 
typedef struct _cooMat {
  int nrows,  /**< number of rows */
      ncols,  /**< number of columns */
      nnz,    /**< number of non-zeros */
      *ir,    /**< row indices of nonzero entries */
      *jc;    /**< column indices of a nonzero entries */
  double *vv; /**< values */
} cooMat;

/*! 
 * @brief sparse matrix format: the compressed sparse row (CSR) format, 0-based
 * 
 * 3-array variant: ia,ja,a, nnz == ia[nrows]
 */ 
typedef struct _csrMat {
  int nrows,   /**< number of rows */
      ncols,   /**< number of columns */
      *ia,     /**<  row pointers (of size nrows+1) */
      *ja;     /**<  column indices (of size nnz) */
  double  *a;  /**<  numeric values (of size nnz) */
} csrMat;

/*!
 * @brief  parameters for polynomial filter 
 *
 * default values are set by set_pol_def
 */
typedef struct _polparams {
  /** @name input to find_pol */
  /**@{*/ 
  int max_deg;        /**< max allowed degree */
  int min_deg ;       /**< min allowed degree */
  int damping;        /**< 0 = no damping, 1 = Jackson, 2 = Lanczos */
  double thresh_ext;  /**< threshold for accepting polynom. for end intervals */
  double thresh_int;  /**< threshold for interior intervals */
  /**@}*/
  
  /** @name output from find_pol */
  /**@{*/ 
  double *mu;         /**< coefficients. allocation done by set_pol */
  double cc;          /**< center of interval - used by chebAv */
  double dd;          /**< half-width of interval - used by chebAv */
  double gam;         /**< center of delta function used */
  double bar;         /**< p(theta)>=bar indicates a wanted Ritz value */
  /**@}*/
  
  /** @name both input to and output from find_pol */
  /**@{*/ 
  int deg ;           /**< if deg == 0 before calling find_deg then
                       the polynomial degree is  computed
                       internally. Otherwise it is of degree deg.
                       [and  thresh_ext and thresh_int are not used]
                       default value=0, set by call to set_pol_def */
  /**@}*/
} polparams;

/**
 * @brief linear solver function prototype: [complex version]
 *
 * n  is the size  of the system,  br, bz are  the right-hand
 * side (real and  imaginary parts of complex vector),  xr, xz will
 * be the  solution (complex vector),  and "data" contains  all the
 * data  needed  by  the  solver. 
 */
typedef void (*linSolFunc)(int n, double *br, double *bz, double *xr, double *xz, void *data);

/** 
 * @brief function prototype for applying the following operations with LB
 *   y = LB  \ x 
 *   y = LB' \ x
 *   y = LB  * x
 *   y = LB' * x
 */
typedef void (*LBFunc)(double *x, double *y, void *data);

/**
 * @brief matvec function prototype 
 */
typedef void (*MVFunc)(double *x, double *y, void *data);

/*!
 * @brief  parameters for rational filter
 *
 * default values are set by set_rat_def
 */
typedef struct _ratparams {
  /*  */
  int num;            /**< number of the poles */
  int pw;             /**< default multiplicity of each pole */
  int method;         /**< type of poles: 0: Cauchy pole, 1: Mid-point */
  double beta;        /**< LS approximation weight */
  double aa;          /**< left endpoint of the interval */
  double bb;          /**< right endpoint of the interval */
  double bar;         /**< rational function value at boundaries */
  /** internal data */
  int *mulp;          /**< multiplicity of each pole */
  int pow;            /**< total multiplicites of all poles */
  /** The following are output - i.e., set by find_ratf */
  complex double *omega; /**< weights allocation done by find_ratf */
  complex double *zk;    /**< locations of poles done by find_ratf */
  //double cc;          // center of interval
  //double dd;          // half-width of interval
  int use_default_solver; /**< if default solver is used */
  /** function and associated data to solve shifted linear system with A-\sigma B */
  linSolFunc *solshift; /**< arrays of function pointers of length `num' */
  void **solshiftdata;  /**< arrays of (void*) of length `num' */ 
} ratparams;


/*!
 * @brief user-provided Mat-Vec function and data for y = A * x
 *
 */
typedef struct _externalMatvec {
  int n;               /**< dimension of A */
  MVFunc func;         /**< function pointer */
  void *data;          /**< data */
} externalMatvec;

/*!
 * @brief wrapper of all global variables in EVSL
 *
 */
typedef struct _evsldata {
  externalMatvec Amatvec;  /**< external matvec routine and the associated data for A */
  int hasB;                /**< if right-hand matrix B is set, i.e., if it is a 
                                generalized eigenvalue problem */
  int isDefaultLB;         /**< if B is factored by the default solver */
  LBFunc LB_mult, LBT_mult, LB_solv, LBT_solv; /**< functions to perform 
                                                    y = LB * x, y = LB' * x,
                                                    y = LB \ x, and y = LB' \ x */
  void *LB_func_data; /**< associated data with LB */
  double *LB_func_work;  /**< work space for performing LBFunc,
                              for matvec_gen, size of n
                              for solve (A-SIGMA B), size of 2*n */
} evslData;

/* global variable: evslData */
extern evslData evsldata;

#endif
