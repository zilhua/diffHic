#include "diffhic.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

extern "C" { 

static const R_CallMethodDef all_call_entries[] = {
	CALLDEF(check_input, 2),
	CALLDEF(cap_input, 3),
	CALLDEF(cluster_2d, 5),
	CALLDEF(split_clusters, 6),
	CALLDEF(collect_background, 5),
	CALLDEF(count_connect, 5),
	CALLDEF(count_marginals, 3),
	CALLDEF(count_patch, 3),
	CALLDEF(iterative_correction, 9),
	CALLDEF(report_hic_pairs, 11),
    CALLDEF(pair_stats, 9),
  	{NULL, NULL, 0}
};

void attribute_visible R_init_diffHic(DllInfo *info)
{
	R_registerRoutines(info, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}

}
