#ifndef WRAPPERS_H
#define WRAPPERS_H

#ifdef USING_R

#include <Rcpp.h>

#include "arbor.h"

arbor::qsm::QSM as_qsm(Rcpp::DataFrame df);
arbor::qsm::QSF as_qsf(Rcpp::List x);
Rcpp::DataFrame as_dataframe(const arbor::qsm::QSM& qsm);
Rcpp::DataFrame as_dataframe(const PointCloud& cloud);

#endif

#endif
