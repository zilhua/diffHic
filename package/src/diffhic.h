#ifndef HICCUP_H
#define HICCUP_H

#include <deque>
#include <queue>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cmath>

#include "R.h"
#include "Rinternals.h"

extern "C" {

SEXP check_input(SEXP, SEXP);

SEXP cap_input(SEXP, SEXP, SEXP);

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
