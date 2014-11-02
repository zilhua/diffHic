#include "diffhic.h"

SEXP collect_background (SEXP anchor, SEXP target, SEXP counts, SEXP flank, SEXP tlen) try {
	if (!isInteger(anchor) || !isInteger(target)) { throw std::runtime_error("anchor/target vectors must be integer"); }
	const int npair=LENGTH(anchor);
	if (LENGTH(target)!=npair) { throw std::runtime_error("anchor/target vectors must have the same length"); }
	if (!isInteger(counts)) { throw std::runtime_error("count matrix should be integer"); }
	const int nlibs=LENGTH(counts)/npair;
	if (nlibs*npair!=LENGTH(counts)) { throw std::runtime_error("dimensions of count matrix are not consistent with number of pairs"); }

	// Setting up pointers.
	const int * aptr=INTEGER(anchor),
		  * tptr=INTEGER(target);
	const int ** cptrs=(const int**)R_alloc(nlibs, sizeof(int*));
	cptrs[0]=INTEGER(counts);
	for (int i=1; i<nlibs; ++i) { cptrs[i]=cptrs[i-1]+npair; }

	// Determining the flank width and anchor/target lengths.
	if(!isInteger(flank) || LENGTH(flank)!=1) { throw std::runtime_error("flank width must be an integer scalar"); }
	const int flank_width=asInteger(flank);
	if (!isInteger(tlen) || LENGTH(tlen)!=1) { throw std::runtime_error("target length must be an integer scalar"); }
	const int target_len=asInteger(tlen);

	SEXP output=PROTECT(allocMatrix(INTSXP, npair, nlibs));
try {
	// Setting up output constructs.
	int** optrs=(int**)R_alloc(nlibs, sizeof(int*));
	optrs[0]=INTEGER(output);
	for (int i=1; i<nlibs; ++i) { optrs[i]=optrs[i-1]+npair; }

	int left_index=0, right_index=0;
    int left_edge, right_edge, over_len, lib;
	for (int i=0; i<npair; ++i) {
		const int& cura=aptr[i];
		const int& curt=tptr[i];
		if (i) { 
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]=optrs[lib][i-1]; }
		} else {
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]=0; }
		}

		// Figuring out the effective left and right edges. 
		left_edge = curt - flank_width;
		right_edge = curt + flank_width + 1;
		over_len = right_edge - target_len;
		if (left_edge < 0) { right_edge -= left_edge; }
 	   	if (over_len > 0) { left_edge -= over_len; }

		// Computing the running sum at each entry.
		while (aptr[left_index] < cura || (aptr[left_index]==cura && tptr[left_index] < left_edge)) { 
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]-=cptrs[lib][left_index]; }
			++left_index;	
		}
		while (right_index < npair && cura==aptr[right_index] && right_edge > tptr[right_index]) {
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]+=cptrs[lib][right_index]; }
			++right_index;	
		}
	}
} catch (std::exception& e) {
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}
