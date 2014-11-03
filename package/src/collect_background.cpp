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

	// Determining the flank width and target lengths.
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

/* This is an equivalent function to run through the diagonals. 
   It looks superficially similar to the function above, but it
   handles the edge cases differently.
 */

SEXP collect_diagonal (SEXP diagonal, SEXP position, SEXP counts, SEXP flank, SEXP clen) try {
	if (!isInteger(diagonal) || !isInteger(position)) { throw std::runtime_error("diagonal/position vectors must be integer"); }
	const int npair=LENGTH(diagonal);
	if (LENGTH(position)!=npair) { throw std::runtime_error("diagonal/position vectors must have the same length"); }
	if (!isInteger(counts)) { throw std::runtime_error("count matrix should be integer"); }
	const int nlibs=LENGTH(counts)/npair;
	if (nlibs*npair!=LENGTH(counts)) { throw std::runtime_error("dimensions of count matrix are not consistent with number of pairs"); }

	// Setting up pointers.
	const int * dptr=INTEGER(diagonal),
		  * pptr=INTEGER(position);
	const int ** cptrs=(const int**)R_alloc(nlibs, sizeof(int*));
	cptrs[0]=INTEGER(counts);
	for (int i=1; i<nlibs; ++i) { cptrs[i]=cptrs[i-1]+npair; }

	// Determining the flank width.
	if(!isInteger(flank) || LENGTH(flank)!=1) { throw std::runtime_error("flank width must be an integer scalar"); }
	const int flank_width=asInteger(flank),
		  eff_full_flank = flank_width*2 + 1;
	if (!isInteger(clen) || LENGTH(clen)!=1) { throw std::runtime_error("chromosome length must be an integer scalar"); }
	const int chromo_len=asInteger(clen);

	// Setting up the remedial matrices.
	std::deque<int> remedial_start(eff_full_flank, -1), remedial_end(eff_full_flank, -1);

	SEXP output=PROTECT(allocVector(VECSXP, 2));
try {
	SET_VECTOR_ELT(output, 0, allocMatrix(INTSXP, npair, nlibs));
	SET_VECTOR_ELT(output, 1, allocVector(INTSXP, npair));

	// Setting up output pointers.
	int** optrs=(int**)R_alloc(nlibs, sizeof(int*));
	optrs[0]=INTEGER(VECTOR_ELT(output, 0));
	for (int i=1; i<nlibs; ++i) { optrs[i]=optrs[i-1]+npair; }
	int* nptr=INTEGER(VECTOR_ELT(output, 1));

	int left_index=0, right_index=0;
    int left_edge, right_edge, over_len, lib, eff_chromo_len;
   	int	remedial_diff, remedial_i=npair;

	for (int i=0; i<npair; ++i) {
		const int& curdiag=dptr[i];
		const int& curpos=pptr[i];
		if (i) { 
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]=optrs[lib][i-1]; }
		} else {
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]=0; }
		}
		nptr[i] = eff_full_flank;

		/* Figuring out the effective left and right edges. Note the truncation of the
		 * chromosome length as we get to the outer diagonals. This is because the diagonals
		 * get shorter as we get further away from the center.
		 */
		eff_chromo_len = chromo_len - curdiag;
		left_edge = curpos - flank_width;
		right_edge = curpos + flank_width + 1;
		over_len = right_edge - eff_chromo_len;
		if (left_edge < 0) { right_edge -= left_edge; }
 	   	if (over_len > 0) { left_edge -= over_len; }

		// Computing the running sum at each entry.
		while (dptr[left_index] < curdiag || (dptr[left_index]==curdiag && pptr[left_index] < left_edge)) { 
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]-=cptrs[lib][left_index]; }
			++left_index;
		}
		while (right_index < npair && curdiag==dptr[right_index] && right_edge > pptr[right_index]) {
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i]+=cptrs[lib][right_index]; }
			++right_index;
		}

		// Storing some data, to use in remedial activities.
		remedial_diff = eff_full_flank - eff_chromo_len;
		if (remedial_diff >= 0) {
			if (curpos==0) {
				remedial_start[remedial_diff]=i;
 			} else if (curpos==eff_chromo_len-1) { 
				remedial_end[remedial_diff]=i;
			}
			if (remedial_i==npair && remedial_diff > 0) { remedial_i=i; }
		}
	}

	/* Remedial action is required for some pairs. If the background area does
 	 * not fit within the single diagonal, it is rerouted to the edges of the
 	 * map, i.e., for a diagonal starting at (x, 1) and ending at (n, 1+n-x),
 	 * we start to count elements at (x-1, 1), (x-2, 1), ... along with (n,
 	 * n-x), (n, n-x-1), ... where n is the total number of locks in the
 	 * chromosome.
	 */
	int remedial_left, remedial_right, remedial_temp;
	for (int i=remedial_i; i<npair; ++i) {
		const int& curdiag=dptr[i];
		const int& curpos=pptr[i];
		eff_chromo_len = chromo_len - curdiag;
		remedial_diff = eff_full_flank - eff_chromo_len;

		/* Some finessing is done here to ensure that the number of extra
 		 * counts added are evenly distributed on both left and right
 		 * edges.  Note that eff_chromo_len must be even if remedial_diff
 		 * is odd, as eff_full_flank is always odd; integer solutions guaranteed.
		 */
		remedial_left = remedial_right = remedial_diff/2;
		if (remedial_diff % 2 == 1) {
			if (curpos <= eff_chromo_len/2) {
				++remedial_left;
			} else {
				++remedial_right;
			}
		}
//		Rprintf("%i %i %i %i\n", curdiag, curpos, remedial_left, remedial_right);
		
		// Checking that it doesn't try to get stuff past the diagonal.
		remedial_temp = remedial_left - curdiag;
		if (remedial_temp > 0) {
			nptr[i] -= remedial_temp;
			remedial_left = curdiag;
		}
		remedial_temp = remedial_right - curdiag;
		if (remedial_temp > 0) {
			nptr[i] -= remedial_temp;
			remedial_right = curdiag;
		}

		// Adding the additional counts.
		for (remedial_i=0; remedial_i < remedial_left; ++remedial_i) {
			const int& remix = remedial_start[remedial_diff - remedial_i - 1];
			if (remix < 0) { continue; }
//			Rprintf("\tLeft adding %i %i\n", dptr[remix], pptr[remix]);
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i] += cptrs[lib][remix]; }
		} 
		for (remedial_i=0; remedial_i < remedial_right; ++remedial_i) {
			const int& remix = remedial_end[remedial_diff - remedial_i - 1];
			if (remix < 0) { continue; }
//			Rprintf("\tRight adding %i %i\n", dptr[remix], pptr[remix]);
			for (lib=0; lib<nlibs; ++lib) { optrs[lib][i] += cptrs[lib][remix]; }
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
