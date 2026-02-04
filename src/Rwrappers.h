#ifndef WRAPPERS_H
#define WRAPPERS_H


#include <Rcpp.h>
#include "QSM.h"
#include "QSF.h"

QSM as_qsm(Rcpp::DataFrame df);
QSF as_qsf(Rcpp::List x);
Rcpp::DataFrame as_dataframe(const QSM& qsm);

#endif
