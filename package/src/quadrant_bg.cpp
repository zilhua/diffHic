#include "diffhic.h"

/* Setting up a couple of structures to direct how the additions are performed.
 */

struct basic {
	basic(int w, int t, bool i) : width(w), level(0), intra(i), tlen(t) {}
	virtual void set (int, int)=0;
	virtual ~basic() {};
	virtual bool bump_level () { 
		if (level >= width) { return false; }
		++level;
		return true;
	}
	int row, left, right;
protected:
	int level, width, tlen;
	bool intra;

	void restrain () {
		if (intra) {
			if (right > row) { 
				right=row+1; 
				if (left > right) { left=right; }
			}
		} else if (right > tlen) { right=tlen; } // For intra's, right will hit diagonal; no need to worry about tlen.
	}
	void set_lefty (int a, int t) {
		left=t-width;
		if (left < 0) { left=0; }
		right=(level ? t+1 : t);
		restrain();
	}
	void set_righty (int a, int t) {
		left=(level ? t : t+1); 
		right=t+width+1; 
		restrain();
	}
	void set_top (int a) { row=a+level; }
	void set_bottom (int a) { row=a-level; }
};

struct upperleft : public basic {
	upperleft(int w, int t, bool i) : basic(w, t, i) {}
	~upperleft() {};
	void set(int a, int t) {
		set_top(a);
		set_lefty(a, t);
	}
};

// For each quadrant.

struct upperright : public basic {
	upperright(int w, int t, bool i) : basic(w, t, i) {}
	~upperright() {};
	void set(int a, int t) {
		set_top(a);
		set_righty(a, t);
	}
};

struct lowerleft : public basic {
	lowerleft(int w, int t, bool i) : basic(w, t, i) {}
	~lowerleft() {};
	void set(int a, int t) {
		set_bottom(a);
		set_lefty(a, t);
	}
};

struct lowerright : public basic { 
	lowerright(int w, int t, bool i) : basic(w, t, i) {}
	~lowerright() {};
	void set(int a, int t) {
		set_bottom(a);
		set_righty(a, t);
	}
};

// Same again, for each horizontal/vertical extension.

struct upper : public basic {
	upper(int w, int t, bool i) : basic(w, t, i) { level=1; }
	~upper() {};
	void set(int a, int t) {
		set_top(a);
		left=t;
		right=t+1; // No need to restrain, we're moving away from the diagonal.
	}	
};

struct lower : public basic {
	lower(int w, int t, bool i) : basic(w, t, i) { level=1; }
	~lower() {};
	void set(int a, int t) { 
		set_bottom(a);
		left=t;
		right=t+1;
		restrain(); 
	}	
};

struct left : public basic  {
	left(int w, int t, bool i) : basic(w, t, i) {}
	~left() {};
	bool bump_level() { return false; } // First bump is the last.
	void set(int a, int t) {
		row=a;
		set_lefty(a, t);
	}
};

struct right : public basic { 
	right(int w, int t, bool i) : basic(w, t, i) {}
	~right() {};
	bool bump_level() { return false; }
	void set(int a, int t) {
		row=a;
		set_righty(a, t);
	}
};

/* Main loop */

extern "C" {

SEXP quadrant_bg (SEXP anchor, SEXP target, 
		SEXP abundance_int, SEXP abundance_dec, SEXP mult,
		SEXP width, SEXP alen, SEXP tlen, SEXP issame) try {
	if (!isInteger(anchor) || !isInteger(target)) { throw std::runtime_error("anchor/target vectors must be integer"); }
	const int npair=LENGTH(anchor);
	if (LENGTH(target)!=npair) { throw std::runtime_error("anchor/target vectors must have the same length"); }
	if (!isInteger(abundance_int) || !isInteger(abundance_dec)) { throw std::runtime_error("vector of abundances should be integer"); }
	if (LENGTH(abundance_int)!=npair || LENGTH(abundance_dec)!=npair) { 
		throw std::runtime_error("vector of abundances should be the same length as that of the indices"); }
  	
	// Setting up pointers.
	const int * aptr=INTEGER(anchor), 
	  * tptr=INTEGER(target),
	  * biptr=INTEGER(abundance_int),
	  * bdptr=INTEGER(abundance_dec);

	// Determining the flank width.
	if (!isReal(mult) || LENGTH(mult)!=1) { throw std::runtime_error("multiplier must be a double-precision scalar"); }
	const double multiplier=asReal(mult);
	if (!isInteger(width) || LENGTH(width)!=1) { throw std::runtime_error("flank width must be an integer scalar"); }
	const int flank_width=asInteger(width);

	if (!isInteger(alen) || LENGTH(alen)!=1) { throw std::runtime_error("anchor length must be an integer scalar"); }
	const int alength=asInteger(alen);
	if (!isInteger(tlen) || LENGTH(tlen)!=1) { throw std::runtime_error("anchor length must be an integer scalar"); }
	const int tlength=asInteger(tlen);
	if (!isLogical(issame) || LENGTH(issame)!=1) { throw std::runtime_error("same chromosome specifier must be a logical scalar"); }
	const bool intrachr=asLogical(issame);

	SEXP output=PROTECT(allocVector(REALSXP, npair));
	try {
		double * optr=REAL(output);
		int curpair;
		for (curpair=0; curpair<npair; ++curpair) { optr[curpair]=0; }
		int running_sum_int, running_sum_dec, 
			left_index, right_index,
			left_edge, right_edge, cur_anchor; 

		int* nptr=(int*)R_alloc(npair, sizeof(int));
		double* temp_int=(double*)R_alloc(npair, sizeof(double));
		double* temp_dec=(double*)R_alloc(npair, sizeof(double));
		double temp_val;

		// Iterating over all quadrants.
		upperleft ul(flank_width, tlength, intrachr);
		upperright ur(flank_width, tlength, intrachr);
		lowerleft ll(flank_width, tlength, intrachr);
		lowerright lr(flank_width, tlength, intrachr);		
		upper up(flank_width, tlength, intrachr);
		lower lo(flank_width, tlength, intrachr);
		left le(flank_width, tlength, intrachr);
		right ri(flank_width, tlength, intrachr);
		basic* current=NULL;

		for (int quadtype=0; quadtype<8; ++quadtype) {
			switch(quadtype) { 
				case 0: current=&ul; break;
				case 1: current=&ur; break;
				case 2: current=&ll; break;
				case 3: current=&lr; break;
				case 4: current=&up; break;
				case 5: current=&lo; break;
				case 6: current=&le; break;
				case 7: current=&ri; break;
			}
			for (curpair=0; curpair<npair; ++curpair) { 
				nptr[curpair]=0; 
				temp_int[curpair]=0;
				temp_dec[curpair]=0;
			}

			// Iterating across all flank widths.
			do {
				running_sum_int=0;
				running_sum_dec=0;
				left_index=0;
				right_index=0;
				
				for (curpair=0; curpair<npair; ++curpair) {
					current->set(aptr[curpair], tptr[curpair]);
					cur_anchor=current->row;
					if (cur_anchor >= alength) { break; }
					left_edge=current->left;
					right_edge=current->right;

					// Identifying all bin pairs in the relevant range.
					while (left_index < npair && (aptr[left_index] < cur_anchor || 
							(aptr[left_index]==cur_anchor && tptr[left_index] < left_edge))) {
						running_sum_int -= biptr[left_index];
						running_sum_dec -= bdptr[left_index];
						++left_index;
					}

					while (right_index < npair && (aptr[right_index]<cur_anchor || 
							(aptr[right_index]==cur_anchor && tptr[right_index] < right_edge))) { 
						running_sum_int += biptr[right_index];
						running_sum_dec += bdptr[right_index];
						++right_index;		
					}

					if (cur_anchor >= 0) {
						// Figuring out the actual number of boxes.
						temp_int[curpair] += running_sum_int;
						temp_dec[curpair] += running_sum_dec;
						nptr[curpair] += right_edge - left_edge;
					}
				}
			} while (current->bump_level());

			// Checking if it exceeds the previous maxima.
			for (curpair=0; curpair<npair; ++curpair) {
// 			    if (quadtype==4) { Rprintf("%i %i %i\n", curpair+1, temp_int[curpair], temp_dec[curpair]); }	
				if (nptr[curpair]) {
					temp_val = (temp_int[curpair] + temp_dec[curpair]/multiplier)/nptr[curpair];
					if (optr[curpair] < temp_val) { optr[curpair]=temp_val; }
				}
			}
		}
	} catch (std::exception &e) {
		UNPROTECT(1);
		throw;
	}
	   
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}

}
