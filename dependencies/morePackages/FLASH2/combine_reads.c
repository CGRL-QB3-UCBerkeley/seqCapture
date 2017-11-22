/*
 * combine_reads.c:  This file contains the code implementing the core algorithm
 * to combine reads in FLASH.
 */

/*
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013, 2014 Eric Biggers
 *
 * This file is part of FLASH, a fast tool to merge overlapping paired-end
 * reads.
 *
 * FLASH is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 *
 * FLASH is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FLASH; if not, see http://www.gnu.org/licenses/.
 */

#include "combine_reads.h"
#include "read.h"
#include "util.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#if defined(__GNUC__) && defined(__SSE2__)
#  define WITH_SSE2
#endif

#ifdef WITH_SSE2
#  include <xmmintrin.h>
#  include <emmintrin.h>
#endif

#ifdef WITH_SSE2

/* Sum the values an 8 x 8 bit vector and return a 32-bit result.  */
static inline uint32_t
hsum32_v8(__m128i v)
{
	v = _mm_sad_epu8(v, _mm_set1_epi8(0));
	return (uint32_t)_mm_extract_epi16(v, 0) +
	       (uint32_t)_mm_extract_epi16(v, 4);
}

/* Sum the values an 8 x 16 bit vector and return a 32-bit result.  */
static inline uint32_t
hsum32_v16(__m128i v)
{
	__m128i mask = _mm_set1_epi32(0x0000ffff);
	v = _mm_add_epi32(v & mask, _mm_srli_si128(v, 2) & mask);
	v = _mm_add_epi32(v, _mm_srli_si128(v, 4));
	v = _mm_add_epi32(v, _mm_srli_si128(v, 8));
	return _mm_cvtsi128_si32(v);
}

#endif /* WITH_SSE2 */

/*
 * Compute mismatch statistics between two sequences.
 *
 * @seq_1, @seq_2:
 *	The two sequences to compare (ASCII characters A, C, G, T, N).
 * @qual_1, @qual_2:
 *	Quality scores for the two sequences, based at 0.
 * @haveN
 *	As an optimization, this can be set to %false to indicate that neither
 *	sequence contains an uncalled base (represented as an N character).
 * @len_p
 *	Pointer to the length of the sequence.  This value will be updated to
 *	subtract the number of positions at which an uncalled base (N) exists in
 *	either sequence.
 * @num_mismatches_ret
 *	Location into which to return the number of positions at which the bases
 *	were mismatched.
 * @mismatch_qual_total_ret
 *	Location into which to return the sum of lesser quality scores at
 *	mismatch sites.
 */

static inline void
compute_mismatch_stats(const char * restrict seq_1,
		       const char * restrict seq_2,
		       const char * restrict qual_1,
		       const char * restrict qual_2,
		       bool haveN,
		       int * restrict len_p,
		       unsigned * restrict num_mismatches_ret,
		       unsigned * restrict mismatch_qual_total_ret)
{
	int num_uncalled = 0;
	unsigned num_mismatches = 0;
	unsigned mismatch_qual_total = 0;
	int len = *len_p;

	if (haveN) {
		for (int i = 0; i < len; i++) {
			if (seq_1[i] == 'N' || seq_2[i] == 'N') {
				num_uncalled++;
			} else {
				if (seq_1[i] != seq_2[i])  {
					num_mismatches++;
					mismatch_qual_total += min(qual_1[i], qual_2[i]);
				}
			}
		}
	} else {
		/* This part of the 'if' statement is for optimization purposes
		 * only; its behavior is equivalent to the block above, except
		 * this block assumes there are no N characters in the input,
		 * and therefore no further checks for N's are needed.
		 *
		 * Note: this optimization is only useful if most reads don't
		 * contain N characters.  */

	#ifdef WITH_SSE2

		/* Optional vectorized implementation (about twice as fast as
		 * nonvectorized on x86_64).  */

		while (len >= 16) {

			/* 16 x 8 bit counters for number of mismatches  */
			__m128i num_mismatches_v8 = _mm_set1_epi8(0);

			/* 8 x 16 bit counters for mismatch quality total  */
			__m128i mismatch_qual_total_v16 = _mm_set1_epi16(0);

			/* The counters of num_mismatches_v8 will overflow if
			 * 256 mismatches are detected at the same position
			 * modulo 16 bytes.  So, don't process 4096 or more
			 * bytes before reducing the counters.
			 *
			 * mismatch_qual_total_v16 would overflow even faster,
			 * but we use 16-bit counters for it.  */

			int todo = min(len, 255 * 16) & ~0xf;
			len -= todo;

			do {
				/* Load 16 bases  */
				__m128i s1_v8 = _mm_loadu_si128((const void *)seq_1);
				__m128i s2_v8 = _mm_loadu_si128((const void *)seq_2);

				/* Load 16 quality scores  */
				__m128i q1_v8 = _mm_loadu_si128((const void *)qual_1);
				__m128i q2_v8 = _mm_loadu_si128((const void *)qual_2);

				/* Compare bases with each other and negate the
				 * result.  This will produce 0xff in bytes
				 * where the bases differ and 0x00 in bytes
				 * where the bases were the same.  */
				__m128i cmpresult = ~_mm_cmpeq_epi8(s1_v8, s2_v8);

				/* Tally mismatched bases.  Subtracting 0x00 and
				 * 0xff is equivalent to adding 0 and 1,
				 * respectively.  */
				num_mismatches_v8 = _mm_sub_epi8(num_mismatches_v8,
								 cmpresult);

				/* Tally quality scores for mismatched bases.
				 */

				/* Get minimum of each quality score.  */
				__m128i qmin_v8 = _mm_min_epu8(q1_v8, q2_v8);

				/* Select only quality scores at mismatch sites
				 */
				__m128i qadd_v8 = qmin_v8 & cmpresult;

				/* Double the precision (8 => 16 bits) and tally  */
				__m128i qadd_v16_1 = _mm_unpacklo_epi8(qadd_v8,
								       _mm_set1_epi8(0));

				__m128i qadd_v16_2 = _mm_unpackhi_epi8(qadd_v8,
								       _mm_set1_epi8(0));

				mismatch_qual_total_v16 = _mm_add_epi16(mismatch_qual_total_v16,
									qadd_v16_1);

				mismatch_qual_total_v16 = _mm_add_epi16(mismatch_qual_total_v16,
									qadd_v16_2);

				/* Advance pointers  */
				seq_1 += 16, seq_2 += 16;
				qual_1 += 16, qual_2 += 16;
				todo -= 16;
			} while (todo);

			/* Reduce the counters.  */
			num_mismatches += hsum32_v8(num_mismatches_v8);
			mismatch_qual_total += hsum32_v16(mismatch_qual_total_v16);
		}

	#endif /* WITH_SSE2  */

	#if 0
		/* Verify the values computed by the vectorized implementation.
		 */
		{
			int veclen = *len_p & ~0xf;
			const char *_seq_1 = seq_1 - veclen;
			const char *_seq_2 = seq_2 - veclen;
			const char *_qual_1 = qual_1 - veclen;
			const char *_qual_2 = qual_2 - veclen;
			unsigned _num_mismatches = 0;
			unsigned _mismatch_qual_total = 0;
			for (int i = 0; i < veclen; i++) {
				if (_seq_1[i] != _seq_2[i])  {
					_num_mismatches++;
					_mismatch_qual_total += min(_qual_1[i], _qual_2[i]);
				}
			}
			assert(num_mismatches == _num_mismatches);
			assert(mismatch_qual_total == _mismatch_qual_total);
		}
	#endif

		/* Process any remainder that wasn't processed by the vectorized
		 * implementation.  */
		for (int i = 0; i < len; i++) {
			if (seq_1[i] != seq_2[i])  {
				num_mismatches++;
				mismatch_qual_total += min(qual_1[i], qual_2[i]);
			}
		}
	}

	/* Return results in pointer arguments  */
	*num_mismatches_ret = num_mismatches;
	*mismatch_qual_total_ret = mismatch_qual_total;
	*len_p -= num_uncalled;
}

#define NO_ALIGNMENT INT_MIN



static inline int *
pair_align(const struct read *read_1, const struct read *read_2,
	   int min_overlap, int max_overlap, int min_overlap_outie, float max_mismatch_density,
	   bool allow_outies, bool * was_outie)
{
	bool haveN = memchr(read_1->seq, 'N', read_1->seq_len) ||
		     memchr(read_2->seq, 'N', read_2->seq_len);

	/* Best (smallest) mismatch density that has been found so far in an
	 * overlap. */
	float best_mismatch_density = max_mismatch_density + 1.0f;
	float best_qual_score = 0.0f;
	int best_position = NO_ALIGNMENT;
	bool best_was_outie = false;
	bool doing_outie = false;
	int start;
	int end;
	int read_offset = 0;
	int best_offset = 0;
again:
	/* Require at least min_overlap bases overlap, and require that the
	 * second read is not overlapped such that it is completely contained in
	 * the first read.  */
	start = max(0, read_1->seq_len - read_2->seq_len);
	end = read_1->seq_len - min_overlap + 1;

	if (start == 0) {
		read_offset = (read_2->seq_len - read_1->seq_len);

		if (read_offset < 0) {
			read_offset *= -1;
		}
	}


	for (int i = start; i < end; i++) {
		unsigned num_mismatches;
		unsigned mismatch_qual_total;
		int overlap_len = read_1->seq_len - i;

		/*This modification allows for engulf cases of the read*/
		compute_mismatch_stats(read_1->seq + i,
				       read_2->seq + read_offset,
				       read_1->qual + i,
				       read_2->qual + read_offset,
				       haveN,
				       &overlap_len,
				       &num_mismatches,
				       &mismatch_qual_total);

		/*Logic to make sure read makes minimum requirements*/
		if (((!doing_outie && overlap_len >= min_overlap) || (doing_outie && overlap_len >= min_overlap_outie))) {
			float score_len = (float)min(overlap_len, max_overlap);
			float qual_score = mismatch_qual_total / score_len;
			float mismatch_density = num_mismatches / score_len;
				
			if (mismatch_density <= best_mismatch_density &&
			    (mismatch_density < best_mismatch_density ||
			     qual_score < best_qual_score))
			{
				best_qual_score       = qual_score;
				best_mismatch_density = mismatch_density;
				best_position         = i;
				best_offset           = read_offset;
				best_was_outie        = doing_outie;
			}
		}
		if (read_offset != 0) {
			read_offset -= (i + 1);
			i--;
		}
		
	}

	/*At one point I will remove the goto statments - in non-modified code*/

	if (allow_outies) {
        const struct read *tmp = read_1;
		read_1 = read_2;
		read_2 = tmp;
		allow_outies = false;
		doing_outie = true;
		goto again;
	}

	int *position_and_offset = (int*)malloc(2 * sizeof(int));

	position_and_offset[0] = best_position;
	position_and_offset[1] = best_offset;


	if (best_mismatch_density > max_mismatch_density) {
		position_and_offset[0] = NO_ALIGNMENT;
		position_and_offset[1] = NO_ALIGNMENT;
	}
//		return position_and_offset;

	*was_outie = best_was_outie;
	return position_and_offset;
}

/* Fills in the combined read from the specified alignment.  */
static void
generate_combined_read(const struct read *read_1,
		       const struct read *read_2,
		       struct read *combined_read,
		       int overlap_begin,
		       int read_offset,
	               bool was_outie,
		       bool cap_mismatch_quals)
{
	/* Length of the overlapping part of two reads.  */

	int overlap_len = read_1->seq_len - overlap_begin;
	
	/* Length of the part of the second read not overlapped with the first
	 * read.  */
	int remaining_len = read_2->seq_len - overlap_len;
	int combined_seq_len = read_1->seq_len + remaining_len;

	const char * restrict seq_1 = read_1->seq;
	const char * restrict seq_2 = read_2->seq;
	const char * restrict qual_1 = read_1->qual;
	const char * restrict qual_2 = read_2->qual;
	char * restrict combined_seq;
	char * restrict combined_qual;
	if (was_outie) {
		//Case of outie, Engulf case
		//So, read 2 is engulfed by read 1
		if (read_offset != 0) {
			const struct read *tmp = read_2;
			read_2 = read_1;
			read_1 = tmp;

			seq_1 = read_1->seq;
			seq_2 = read_2->seq;
			qual_1 = read_1->qual;
			qual_2 = read_2->qual;

			overlap_begin = read_offset;
			remaining_len = 0;
			combined_seq_len = read_2->seq_len+read_offset;
			overlap_len = read_2->seq_len;
		} else {
			combined_seq_len = read_1->seq_len - overlap_begin;
			seq_1 += overlap_begin;
			overlap_begin = 0;
			remaining_len = 0;	
		}

    /*'Typical Case*/
    } else {
		//Innie, engulf case
        if (read_offset != 0) {
            //R1 is shorter than R2 (take first part of read R2)

            /*ignore first part of read 2*/
            seq_2 += read_offset;
            qual_2 += read_offset;

            /*don't take any part of the first read*/
            overlap_begin = 0;

            /*Seq_2 minus first part of read 1*/
			combined_seq_len = read_2->seq_len - read_offset;
			remaining_len = read_2->seq_len - read_1->seq_len - read_offset;
           

		} else {
			//nothing fun happens
		}
	}

	/*Allocates the correct size for each condition*/
	
    combined_read->seq = xrealloc(combined_read->seq,
                                              combined_seq_len);
                combined_read->seq_bufsz = combined_seq_len;
		combined_read->qual = xrealloc(combined_read->qual,
					       combined_seq_len);
		combined_read->qual_bufsz = combined_seq_len;

	/* David Modification
 	 * --------|------------
 	 *          ------------|---------
 	 *          Only keeps the overlying parts in outie
 	 */

         combined_seq = combined_read->seq;
	 combined_qual = combined_read->qual;

         combined_read->seq_len = combined_seq_len;
         combined_read->qual_len = combined_seq_len;

	/* Copy the beginning of read 1 (not in the overlapped region).  */
	//if (read_offset == 0) {
		while (overlap_begin--) {
			*combined_seq++ = *seq_1++;
			*combined_qual++ = *qual_1++;
		}
	//} 
	/*else if (read_offset != 0) {
		while (read_offset--) {
			*combined_seq++ = *seq_2++;
			*combined_qual++ = *qual_2++;
		}
		read_offset = -1;
		overlap_len = read_1->seq_len;

	}*/

	
	/* Copy the overlapped region.  */
	while (overlap_len-- > 0) {
			
		if (*seq_1 == *seq_2) {
			/* Same base in both reads.  Take the higher quality
			 * value. */
			*combined_seq = *seq_1;
			*combined_qual = max(*qual_1, *qual_2);
		} else {
			/* Different bases in the two reads; use the higher
			 * quality one.
			 *
			 * The old way of calculating the resulting quality
			 * value (params->cap_mismatch_quals == %true) is to use
			 * the lower quality value, and use a quality value of
			 * at most 2 (+ phred_offset in the final output--- here
			 * the quality values are all scaled to start at 0).
			 * The motivation for this behavior is that the read
			 * combination shows there was sequencing error at the
			 * mismatch location, so the corresponding base call in
			 * the combined read should be given a low quality
			 * score.
			 *
			 * The new way (params->cap_mismatch_quals == %false,
			 * default as of FLASH v1.2.8) is to use the absolute
			 * value of the difference in quality scores, but at
			 * least 2.  This allows a base call with a high quality
			 * score to override a base call with a low quality
			 * score without too much penalty.
			 */

			if (cap_mismatch_quals)
				*combined_qual = min(min(*qual_1, *qual_2), 2);
			else
				*combined_qual = max(abs(*qual_1 - *qual_2), 2);

			if (*qual_1 > *qual_2) {
				*combined_seq = *seq_1;
			} else if (*qual_1 < *qual_2) {
				*combined_seq = *seq_2;
			} else {
				/* Same quality value; take the base from the
				 * first read if the base from the second read
				 * is an 'N'; otherwise take the base from the
				 * second read. */

				if (*seq_2 == 'N')
					*combined_seq = *seq_1;
				else
					*combined_seq = *seq_2;
			}
		}
		combined_seq++;
		combined_qual++;
		seq_1++;
		seq_2++;
		qual_1++;
		qual_2++;
	}
	/* Copy the end of read 2 (not in the overlapped region).  */
	//if (!was_outie && read_offset == 0 && remaining_len > 0) {
		while (remaining_len--) {
			if (read_offset != 0 && was_outie) {
				*combined_seq++ = *seq_2++;
				*combined_qual++ = *qual_2++;
			} else {
				*combined_seq++ = *seq_2++;
				*combined_qual++ = *qual_2++;
			}
		}
	/*This check and sees if there is an offset, but also that the offset wasn't for the innie case*/
	/*} else if (read_offset > 0) {
			while (read_offset--) {
				*combined_seq++ = *seq_2++;
			}
	}*/
}

/* This is the entry point for the core algorithm of FLASH.  The following
 * function attempts to combine @read_1 with @read_2, and writes the result into
 * @combined_read.  COMBINED_AS_INNIE or COMBINED_AS_OUTIE is returned if
 * combination was successful.  COMBINED_AS_OUTIE is only possible if
 * params->allow_outies is set.
 *
 * Note: @read_2 is provided to this function after having been
 * reverse-complemented.  Hence, the code just aligns the reads in the forward
 * orientation, which is equivalent to aligning the original reads in the
 * desired reverse-complement orientation.
 *
 * Please see the help output of FLASH for the description of the min_overlap,
 * max_overlap, and max_mismatch_density parameters.  (--min-overlap,
 * --max-overlap, and --max-mismatch-density on the command line).
 *
 * You may also want to read the original FLASH publication for a description of
 * the algorithm used here:
 *
 *  Title:   FLASH: fast length adjustment of short reads to improve genome assemblies
 *  Authors: Tanja Magoč and Steven L. Salzberg
 *  URL:     http://bioinformatics.oxfordjournals.org/content/27/21/2957.full
 *
 */






/*Function for Junk - removal of primer dimers
 *50 percent of quality scores must be greater the 2 as default*/

double check_for_junk_read(char * qual_1, char * qual_2, int cut_off) {

	double count = 0;
	double junk = 0;
	int i_1 = 0, i_2 = 0;


	/*Could vectorize this in the future*/
	while (qual_1[i_1] != '\0' && qual_2[i_2] != '\0') {
		if (qual_1[i_1] != '\0') {
			if (qual_1[i_1]  <= cut_off) {
				junk++;
			}
			count++;
			i_1++;
		}
	
		if (qual_2[i_2] != '\0') {
			if (qual_2[i_2] <= cut_off) {
				junk++;
			}
			count++;
			i_2++;
		}
	}
	/*Both doubles so do not have to worry about interger division*/
	return junk/count;
}


enum combine_status
combine_reads(const struct read *read_1, const struct read *read_2,
	      struct read *combined_read,
	      const struct combine_params *params)
{
	int overlap_begin, read_offset, *overlap_and_offset;
	enum combine_status status;
	bool was_outie = false;

	/*check read isn't junk;*/
	if (params->discard_reads) {
		double percent_junk = check_for_junk_read(read_1->qual, read_2->qual, params->qual_score_cutoff);
		if (percent_junk*100 > params->percent_cutoff) {
			return JUNK_READ;
		}	
	}
	/* Do the alignment.  */

	overlap_and_offset = pair_align(read_1, read_2,
				   params->min_overlap,
				   params->max_overlap,
				   params->min_overlap_outie,
				   params->max_mismatch_density,
				   params->allow_outies,
				   &was_outie);
	overlap_begin = overlap_and_offset[0];
	read_offset = overlap_and_offset[1];

	/*
	 * If overlap_begin == NO_ALIGNMENT, then no sufficient overlap between
	 * the reads was found.
	 *
	 * If !@was_outie, then the pair forms an "innie" overlap, and
	 * overlap_begin is the 0-based position in read_1 at which read_2
	 * begins.  (Shown below with read 2 already reverse complemented!)
	 *
	 *		0	  overlap_begin
	 *	        |         |
	 *	Read 1: ------------------>
	 *	Read 2:           ---------------------->
	 *
	 * If @was_outie, then the pair forms an "outie" overlap, and
	 * overlap_begin is the 0-based position in read_2 at which read_1
	 * begins. (Shown below with read 2 already reverse complemented!)
	 *
	 *	        0         overlap_begin
	 *	        |         |
	 *	Read 2: ------------------>
	 *	Read 1:           ---------------------->
	 */

	if (overlap_begin == NO_ALIGNMENT)
		return NOT_COMBINED;
	if (!was_outie) {
		status = COMBINED_AS_INNIE;
	} else {
		const struct read *tmp;

		/* Simplify generation of the combined read by turning the outie
		 * case into the innie case.  */

		tmp = read_1;
		read_1 = read_2;
		read_2 = tmp;

		status = COMBINED_AS_OUTIE;
		/*
		 * Now it's just:
		 *
		 *		0	  overlap_begin
		 *	        |         |
		 *	Read 1: ------------------>
		 *	Read 2:           ---------------------->
		 *
		 * The same as the "innie" case.
		 */
	}

    if (params->skip_overlap) {
        status = NOT_COMBINED;
    } else {
	    /* Fill in the combined read.  */
	    generate_combined_read(read_1, read_2, combined_read,
			       overlap_begin, read_offset, was_outie, params->cap_mismatch_quals);
    }
	return status;
}
