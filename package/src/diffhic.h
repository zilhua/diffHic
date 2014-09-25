#ifndef HICCUP_H
#define HICCUP_H

#include <deque>
#include <queue>
#include <map>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "R.h"
#include "Rinternals.h"

const double very_low_value=std::pow(10.0, -8.0);

extern "C" {

SEXP check_input(SEXP, SEXP);

SEXP cluster_2d (SEXP, SEXP, SEXP, SEXP, SEXP); 

SEXP split_clusters (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP collect_background (SEXP, SEXP, SEXP, SEXP);

SEXP count_connect(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP count_marginals (SEXP, SEXP, SEXP);

SEXP count_patch(SEXP, SEXP, SEXP);

SEXP iterative_correction(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP report_hic_pairs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 

SEXP pair_stats (SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

}

#endif
