/*
 * Copyright (C) 2026 Netherlands Forensic Institute
 *
 * This file is part of FDSTools, data analysis tools for Massively
 * Parallel Sequencing of forensic DNA markers.
 *
 * FDSTools is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FDSTools is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FDSTools.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * This file is based on original code from FLASH (Fast Length
 * Adjustment of SHort reads), a fast tool to merge overlapping
 * paired-end reads, version 1.2.11. The original, unmodified code
 * was made available under the below license.
 *
 * Copyright (C) 2012 Tanja Magoc
 * Copyright (C) 2012, 2013, 2014 Eric Biggers
 *
 * FLASH is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * FLASH is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))


#if defined(__SSE2__) || defined(_M_AMD64) || defined(_M_X64) || (defined(_M_IX86_FP) && _M_IX86_FP == 2)
/**************************************************************************************************
  SSE2-enabled Implementation
**************************************************************************************************/
#include <xmmintrin.h>
#include <emmintrin.h>


/*
 * Sum 16 8-bit values and return a 32-bit result.
 */
static inline uint32_t hsum32_v8(__m128i v) {
    v = _mm_sad_epu8(v, _mm_setzero_si128());
    return (uint32_t)_mm_extract_epi16(v, 0) +
           (uint32_t)_mm_extract_epi16(v, 4);
}


/*
 * Sum 8 16-bit values and return a 32-bit result.
 */
static inline uint32_t hsum32_v16(__m128i v) {
    __m128i mask = _mm_set1_epi32(0x0000ffff);
    v = _mm_add_epi32(_mm_and_si128(v, mask), _mm_and_si128(_mm_srli_si128(v, 2), mask));
    v = _mm_add_epi32(v, _mm_srli_si128(v, 4));
    v = _mm_add_epi32(v, _mm_srli_si128(v, 8));
    return _mm_cvtsi128_si32(v);
}


/*
 * Compute mismatch statistics between two sequences.
 *
 * @seq_1, @seq_2:
 *  The two sequences to compare (ASCII characters A, C, G, T, N).
 * @qual_1, @qual_2:
 *  Quality scores for the two sequences, based at 0.
 * @haveN
 *  As an optimization, this can be set to %false to indicate that neither
 *  sequence contains an uncalled base (represented as an N character).
 * @len_p
 *  Pointer to the length of the overlap.  This value will be updated to
 *  subtract the number of positions at which an uncalled base (N) exists in
 *  either sequence.
 * @num_mismatches_ret
 *  Location into which to return the number of positions at which the bases
 *  were mismatched.
 * @mismatch_qual_total_ret
 *  Location into which to return the sum of lesser quality scores at
 *  mismatch sites.
 */
static inline void _compute_mismatch_stats(
        const char * restrict seq_1,
        const char * restrict seq_2,
        const char * restrict qual_1,
        const char * restrict qual_2,
        int haveN,
        int * restrict len_p,
        unsigned * restrict num_mismatches_ret,
        unsigned * restrict mismatch_qual_total_ret) {
    int num_uncalled = 0;
    unsigned num_mismatches = 0;
    unsigned mismatch_qual_total = 0;
    int len = *len_p;

    if (haveN) {
        for (int i = 0; i < len; i++) {
            if (seq_1[i] == 'N' || seq_2[i] == 'N') {
                num_uncalled++;
            } else {
                if (seq_1[i] != seq_2[i]) {
                    num_mismatches++;
                    mismatch_qual_total += min(qual_1[i], qual_2[i]);
                }
            }
        }
    } else {
        while (len >= 16) {

            // 16 x 8 bit counters for number of mismatches
            __m128i num_mismatches_v8 = _mm_setzero_si128();

            // 8 x 16 bit counters for mismatch quality total
            __m128i mismatch_qual_total_v16 = _mm_setzero_si128();

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

                // Load 16 bases and quality scores.
                __m128i s1_v8 = _mm_loadu_si128((const void *)seq_1);
                __m128i s2_v8 = _mm_loadu_si128((const void *)seq_2);
                __m128i q1_v8 = _mm_loadu_si128((const void *)qual_1);
                __m128i q2_v8 = _mm_loadu_si128((const void *)qual_2);

                /* Compare bases with each other and negate the
                 * result.  This will produce 0xff in bytes
                 * where the bases differ and 0x00 in bytes
                 * where the bases were the same.  */
                __m128i cmpresult = _mm_xor_si128(_mm_cmpeq_epi8(s1_v8, s2_v8), _mm_cmpeq_epi32(s1_v8, s1_v8));

                /* Tally mismatched bases.  Subtracting 0x00 and
                 * 0xff is equivalent to adding 0 and 1,
                 * respectively.  */
                num_mismatches_v8 = _mm_sub_epi8(num_mismatches_v8, cmpresult);

                // Tally quality scores for mismatched bases.

                // Get minimum of each quality score.
                __m128i qmin_v8 = _mm_min_epu8(q1_v8, q2_v8);

                // Select only quality scores at mismatch sites.
                __m128i qadd_v8 = _mm_and_si128(qmin_v8, cmpresult);

                // Double the precision (8 => 16 bits) and tally.
                __m128i qadd_v16_1 = _mm_unpacklo_epi8(qadd_v8, _mm_setzero_si128());
                __m128i qadd_v16_2 = _mm_unpackhi_epi8(qadd_v8, _mm_setzero_si128());

                mismatch_qual_total_v16 = _mm_add_epi16(mismatch_qual_total_v16, qadd_v16_1);
                mismatch_qual_total_v16 = _mm_add_epi16(mismatch_qual_total_v16, qadd_v16_2);

                // Advance pointers.
                seq_1 += 16, seq_2 += 16;
                qual_1 += 16, qual_2 += 16;
                todo -= 16;
            } while (todo);

            // Reduce the counters.
            num_mismatches += hsum32_v8(num_mismatches_v8);
            mismatch_qual_total += hsum32_v16(mismatch_qual_total_v16);
        }

        // Process any remainder that wasn't processed by the vectorized implementation.
        for (int i = 0; i < len; i++) {
            if (seq_1[i] != seq_2[i]) {
                num_mismatches++;
                mismatch_qual_total += min(qual_1[i], qual_2[i]);
            }
        }
    }

    // Return results in pointer arguments.
    *num_mismatches_ret = num_mismatches;
    *mismatch_qual_total_ret = mismatch_qual_total;
    *len_p -= num_uncalled;
}


#else

/**************************************************************************************************
  Plain Implementation
**************************************************************************************************/

static inline void _compute_mismatch_stats(
        const char * restrict seq_1,
        const char * restrict seq_2,
        const char * restrict qual_1,
        const char * restrict qual_2,
        int haveN,
        int * restrict len_p,
        unsigned * restrict num_mismatches_ret,
        unsigned * restrict mismatch_qual_total_ret) {
    int num_uncalled = 0;
    unsigned num_mismatches = 0;
    unsigned mismatch_qual_total = 0;
    int len = *len_p;

    if (haveN) {
        for (int i = 0; i < len; i++) {
            if (seq_1[i] == 'N' || seq_2[i] == 'N') {
                num_uncalled++;
            } else {
                if (seq_1[i] != seq_2[i]) {
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
        for (int i = 0; i < len; i++) {
            if (seq_1[i] != seq_2[i]) {
                num_mismatches++;
                mismatch_qual_total += min(qual_1[i], qual_2[i]);
            }
        }
    }

    // Return results in pointer arguments.
    *num_mismatches_ret = num_mismatches;
    *mismatch_qual_total_ret = mismatch_qual_total;
    *len_p -= num_uncalled;
}

#endif

/*************************************************************************************************/

static inline PyObject *_align_pair(
        const char *read_1_seq, const char *read_1_qual, int read_1_len,
        const char *read_2_seq, const char *read_2_qual, int read_2_len,
        int min_overlap, float max_mismatch_density, int allow_outies) {
    int haveN = memchr(read_1_seq, 'N', read_1_len) || memchr(read_2_seq, 'N', read_2_len);

    // Best (smallest) mismatch density that has been found so far in an overlap.
    float best_mismatch_density = max_mismatch_density + 1.0f;
    float best_qual_score = 0.0f;
    int best_position = INT_MIN;
    int best_was_outie = 0;
    int doing_outie = 0;
    int start;
    int end;

again:
    /* Require at least min_overlap bases overlap, and require that the
     * second read is not overlapped such that it is completely contained in
     * the first read.  */
    start = max(0, read_1_len - read_2_len);
    end = read_1_len - min_overlap;
    for (int i = start; i <= end; i++) {
        unsigned num_mismatches;
        unsigned mismatch_qual_total;
        int overlap_len = read_1_len - i;

        _compute_mismatch_stats(
                    read_1_seq + i,
                    read_2_seq,
                    read_1_qual + i,
                    read_2_qual,
                    haveN,
                    &overlap_len,
                    &num_mismatches,
                    &mismatch_qual_total);

        if (overlap_len >= min_overlap) {
            float score_len = (float) overlap_len;
            float qual_score = mismatch_qual_total / score_len;
            float mismatch_density = num_mismatches / score_len;

            if (mismatch_density <= best_mismatch_density &&
                (mismatch_density < best_mismatch_density ||
                 qual_score < best_qual_score)) {
                best_qual_score       = qual_score;
                best_mismatch_density = mismatch_density;
                best_position         = i;
                best_was_outie        = doing_outie;
            }
        }
    }

    if (allow_outies) {
        const char *tmp = read_1_seq;
        read_1_seq = read_2_seq;
        read_2_seq = tmp;
        tmp = read_1_qual;
        read_1_qual = read_2_qual;
        read_2_qual = tmp;
        int tmplen = read_1_len;
        read_1_len = read_2_len;
        read_2_len = tmplen;
        allow_outies = 0;
        doing_outie = 1;
        goto again;
    }

    if (best_mismatch_density > max_mismatch_density) {
        return Py_BuildValue("(OOO)", Py_None, Py_None, Py_None);
    }

    return Py_BuildValue("(iiO)", best_position, best_position - read_1_len,
        best_was_outie? Py_True : Py_False);
}


#define METHOD_DOC "align_pair(read1, read2, min_overlap, max_mismatch_density, allow_outies)\n"\
                   "\n"\
                   "Find overlap in the given read pair.\n"\
                   "\n"\
                   "Input arguments are:\n"\
                   "  read1: tuple(str, str, str)  The header, sequence and quality values of R1.\n"\
                   "  read2: tuple(str, str, str)  The header, sequence and quality values of R2.\n"\
                   "  min_overlap: int             The minimum required number of overlapping bases.\n"\
                   "  max_mismatch_density: float  The maximum mismatch density allowed.\n"\
                   "  allow_outies: bool           Whether the 'outie' configuration is allowed.\n"\
                   "\n"\
                   "Return (overlap_pos: int, overlap_len: int, was_outie: bool), or (None, None, None)\n"\
                   "if the reads could not be combined."


static char *argumentNames[] = {
        "read1", "read2", "min_overlap", "max_mismatch_density", "allow_outies", NULL};


static PyObject *align_pair(PyObject *self, PyObject *args, PyObject *kwargs) {
    void *null; // Sink for ignored arguments.
    const char *read_1_seq, *read_2_seq, *read_1_qual, *read_2_qual;
    Py_ssize_t read_1_len, read_2_len;
    int min_overlap, allow_outies;
    float max_mismatch_density;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(ss#s)(ss#s)ifp", argumentNames,
            &null, &read_1_seq, &read_1_len, &read_1_qual,
            &null, &read_2_seq, &read_2_len, &read_2_qual,
            &min_overlap, &max_mismatch_density, &allow_outies)) {
        return NULL;
    }

    return _align_pair(read_1_seq, read_1_qual, read_1_len,
                       read_2_seq, read_2_qual, read_2_len,
                       min_overlap, max_mismatch_density, allow_outies);
}


static PyMethodDef AlignerMethods[] = {
    {"align_pair", (PyCFunction)align_pair, METH_VARARGS | METH_KEYWORDS, METHOD_DOC},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "pair_align",    /* m_name */
    NULL,            /* m_doc */
    -1,              /* m_size */
    AlignerMethods,  /* m_methods */
    NULL,            /* m_slots */
    NULL,            /* m_traverse */
    NULL,            /* m_clear */
    NULL             /* m_free */
};


PyMODINIT_FUNC PyInit_pair_align(void) {
    return PyModule_Create(&moduledef);
}
