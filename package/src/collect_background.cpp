#include "diffhic.h"

SEXP collect_background (SEXP anchor, SEXP target, SEXP counts, SEXP flank) try {
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

	// Determining the flank width.
	if(!isInteger(flank) || LENGTH(flank)!=1) { throw std::runtime_error("flank width must be an integer scalar"); }
	const int flank_width=asInteger(flank);

	SEXP output=PROTECT(allocMatrix(INTSXP, npair, nlibs));
try {
	int** optrs=(int**)R_alloc(nlibs, sizeof(int*));
	optrs[0]=INTEGER(output);
	for (int i=1; i<nlibs; ++i) { optrs[i]=optrs[i-1]+npair; }

	int left_index=0, right_index=0, lib;
	for (int i=0; i<npair; ++i) {
		const int& cura=aptr[i];
		const int& curt=tptr[i];
		if (i) { 
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]=optrs[lib][i-1]; }
		} else {
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]=0; }
		}

		// Computing the running sum at each entry.
		while (aptr[left_index] < cura || (aptr[left_index]==cura && tptr[left_index] + flank_width < curt)) { 
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]-=cptrs[lib][left_index]; }
			++left_index;	
		}
		while (right_index < npair && cura==aptr[right_index] && curt + flank_width >= tptr[right_index]) {
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
