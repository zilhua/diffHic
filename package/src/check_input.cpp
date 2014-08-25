#include "diffhic.h"

/* This just provides a function that quickly and efficiently
 * checks the incoming pair counts for (a) anchor >= target,
 * (b) ordering of anchor and targets, and (c) a positive count.
 */

SEXP check_input (SEXP anchor, SEXP target, SEXP count) try { 
	if (!isInteger(anchor)) { throw std::runtime_error("anchor should be an integer vector"); }
	if (!isInteger(target)) { throw std::runtime_error("target should be an integer vector"); }
	if (!isInteger(count)) { throw std::runtime_error("count should be an integer vector"); }
	const int nlen=LENGTH(count);
	if (LENGTH(anchor)!=nlen || LENGTH(target)!=nlen) { throw std::runtime_error("vectors should be of the same length"); }

	const int * aptr=INTEGER(anchor),
		  * tptr=INTEGER(target),
		  * cptr=INTEGER(count);

	if (nlen) {
		if (aptr[0] < tptr[0]) { throw std::runtime_error("anchor should be greater than or equal to target"); }
		for (int i=1; i<nlen; ++i) {
			if (aptr[i] < tptr[i]) { throw std::runtime_error("anchor should be greater than or equal to target"); }
			if (aptr[i] < aptr[i-1] || (aptr[i]==aptr[i-1] && tptr[i] < tptr[i-1])) {
				throw std::runtime_error("pairs should be sorted by anchor and target"); 
			}
			if (cptr[i]<=0) { throw std::runtime_error("counts should be positive"); }			
		}
	}
	return ScalarLogical(1);
} catch (std::exception& e){
	return mkString(e.what());
}