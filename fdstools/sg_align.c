/*
 * Copyright (C) 2020 Jerry Hoogenboom
 *
 * This file is based on original code from the TSSV project: Targeted
 * characterisation of short structural variation, version 0.4.0, which was
 * originally made available under the below license.
 *
 * Copyright (c) 2012-2016 Jeroen F.J. Laros <j.f.j.laros@lumc.nl>
 * Copyright (c) 2016 Jerry Hoogenboom <j.hoogenboom@nfi.minvenj.nl>
 * Copyright (c) 2012 Jaap W.F. van der Heijden
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include <Python.h>

#define METHOD_DOC "align(seq1, seq2, indel_socre)\n"\
                   "\n"\
                   "Return the minimum number of mismatches and the corresponding position by\n"\
                   "aligning seq2 against seq1. Return (len(seq2), 0) if len(seq1) < len(seq2).\n"\
                   "\n"\
                   "An insertion or deletion of one base is counted as indel_score mismatches."


/*
Return the minimum of a and b.
*/
static __inline unsigned char _min(const unsigned char a, const unsigned char b) {
    if (a < b) {
        return a;
    }
    return b;
}


/*
Saturating addition of a and b.
*/
static __inline unsigned char _sadd(const unsigned char a, const unsigned char b) {
    if (a > 255 - b) {
        return 255;
    }
    return a + b;
}


#ifdef DEBUG
#include <stdio.h>
#endif


#if defined(_MSC_VER) || defined(__SSE2__)
/**************************************************************************************************
  SSE2-enabled Implementation
**************************************************************************************************/
#include <xmmintrin.h>
#include <emmintrin.h>

#ifdef _MSC_VER
#define _cpuid(cpuInfo, function_id) __cpuidex(cpuInfo, function_id, 0)
#else
#include <cpuid.h>
static __inline void _cpuid(int cpuInfo[4], int function_id) {
    __cpuid_count(function_id, 0, cpuInfo[0], cpuInfo[1], cpuInfo[2], cpuInfo[3]);
}
#endif

#ifdef DEBUG
#define capitalsonly(x) ((x>= 0x41 && x<=0x5A)? x : 0x2D)
static void _print128seq(__m128i var, char* name) {
    char *val = (char*) &var;
    printf("%s: %c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c\n", name,
           capitalsonly(val[0]), capitalsonly(val[1]), capitalsonly(val[2]), capitalsonly(val[3]),
           capitalsonly(val[4]), capitalsonly(val[5]), capitalsonly(val[6]), capitalsonly(val[7]),
           capitalsonly(val[8]), capitalsonly(val[9]), capitalsonly(val[10]), capitalsonly(val[11]),
           capitalsonly(val[12]), capitalsonly(val[13]), capitalsonly(val[14]), capitalsonly(val[15]));
}
static void _print128num(__m128i var, char* name) {
    char *val = (char*) &var;
    printf("%s: %3u %3u %3u %3u %3u %3u %3u %3u %3u %3u %3u %3u %3u %3u %3u %3u\n", name,
           val[0] & 0x0FF, val[1] & 0x0FF, val[2] & 0x0FF, val[3] & 0x0FF,
           val[4] & 0x0FF, val[5] & 0x0FF, val[6] & 0x0FF, val[7] & 0x0FF,
           val[8] & 0x0FF, val[9] & 0x0FF, val[10] & 0x0FF, val[11] & 0x0FF,
           val[12] & 0x0FF, val[13] & 0x0FF, val[14] & 0x0FF, val[15] & 0x0FF);
}
#endif


/*
Allocates and returns a pointer to the alignment matrix.
The actual matrix starts at the first 16-byte aligned byte and is stored diagonally.
The seq1len argument MUST NOT be less than the seq2len argument (this is not checked).
*/
static unsigned char *_sse2_make_matrix(const unsigned int seq1len, const unsigned int seq2len,
                                        const unsigned char indel_score) {
    const unsigned int width = (seq2len+31) & ~0x0F,
                       height = seq1len + seq2len + 1;
    unsigned char *mem = malloc(width * height + 16),
                  *matrix = (unsigned char*)(((Py_uintptr_t)mem + 15) & ~(Py_uintptr_t)0x0F),
                  *cell,
                  score;
    unsigned int i, j;

    // Set the first column to 0.
    for (i = 0, cell = matrix; i <= seq1len; i++, cell += width) {
        *cell = 0;
    }

    // Set the first row to 0, 1, 2, 3, 4, ... times the indel_score
    for (i = 0, cell = matrix, score = 0; i <= seq2len; i++, cell += width + 1) {
        *cell = score;
        score = _sadd(score, indel_score);
        for (j = 1; j <= seq2len - i; j++) {
            *(cell + j) = 255;  // This protects the second row.
        }
    }

    return mem;
}


/*
Reverse a sequence. Copies seq to seqr, from right to left. Len is the length of the sequence,
excluding the terminating NUL byte.
*/
static void _revseq(const char *seq, char *seqr, const unsigned int len){
    const char *p = seq;
    char *q = seqr + len;
    *(q--) = 0;
    while (q >= seqr) {
        *(q--) = *(p++);
    }
}


/*
Fill the alignment matrix as created with _sse2_make_matrix().
*/
static void _sse2_align(unsigned char *mem, const unsigned int seq1len, const unsigned int seq2len,
                        const char *seq1, const char *seq2, const unsigned char indel_score) {
    unsigned int x = 1,
                 y = 1,
                 width = (seq2len+31) & ~0x0F,
                 end = seq1len + _min(16, seq2len),
                 limit;
    unsigned char *matrix = (unsigned char*)(((Py_uintptr_t)mem + 15) & ~(Py_uintptr_t)0x0F),
                  *d = matrix,
                  *l = matrix + width,
                  *i = l + width + 1;

    const __m128i ones = _mm_set1_epi8(1),
        indel_scores = _mm_set1_epi8(indel_score);
    __m128i md = _mm_load_si128((__m128i*)d),
            ml = _mm_load_si128((__m128i*)l),
            mu = _mm_loadu_si128((__m128i*)(l + 1)),
            mi,
            mx,
            my;

    // Get copy of seq2 and reverse of seq1, making sure
    // that we can read 16 bytes (of garbage) past the end.
    char *seq1r = malloc(seq1len + 16),
         *seq2f = malloc(seq2len + 16);
    strcpy(seq2f, seq2);
    _revseq(seq1, seq1r, seq1len);
    mx = _mm_loadu_si128((__m128i*)seq2f);
    my = _mm_loadu_si128((__m128i*)(seq1r + seq1len - 1));

#ifdef DEBUG
    printf("Matrix size %ux%u (%ux%u)\n", seq2len + 1, seq1len + 1, width, seq1len + seq2len + 1);
    _print128seq(mx, "start x");
    _print128seq(my, "start y");
#endif

    while (y < end) {
        mi = _mm_min_epu8(_mm_adds_epu8(_mm_min_epu8(ml, mu), indel_scores),
             _mm_adds_epu8(md, _mm_add_epi8(_mm_cmpeq_epi8(mx, my), ones)));
        _mm_storeu_si128((__m128i*)i, mi);

#ifdef DEBUG
        printf("(%2u,%2u) d=%u l=%u i=%u end=%u\n", x, y, d-matrix, l-matrix, i-matrix, end);
        _print128num(md, "md");
        _print128num(ml, "ml");
        _print128num(mu, "mu");
        _print128num(mi, "mi");
#endif

        if (++y < end) {
            // Move down one row.
            l += width;
            i += width;
            md = ml;
            mu = mi;
            ml = _mm_load_si128((__m128i*)l);

            // Move to the next base in seq1r.
            my = _mm_slli_si128(my, 1);

            limit = y - x;
            if (limit < seq1len) {
                my = _mm_insert_epi16(my, *(short*)(seq1r + seq1len - limit - 1), 0);
            }

#ifdef DEBUG
            _print128seq(my, "shift y");
#endif

        }
        else if ((x += 16) <= seq2len) {
            // Move right 16 columns.
            y = x;
            end += _min(16, seq2len + 1 - x);
            d += 16 * width + 16;
            l = d + width;
            i = l + width + 1;
            md = _mm_load_si128((__m128i*)d);
            ml = _mm_load_si128((__m128i*)l);
            mu = _mm_loadu_si128((__m128i*)(l + 1));
            mx = _mm_loadu_si128((__m128i*)(seq2f + x - 1));
            my = _mm_loadu_si128((__m128i*)(seq1r + seq1len - 1));

#ifdef DEBUG
            _print128seq(mx, "jump  x");
            _print128seq(my, "jump  y");
#endif

        }
        else break;
    }

    free(seq2f);
    free(seq1r);

#ifdef DEBUG
    // Print matrix.
    for (y = 0; y < seq1len + seq2len + 1; y++) {
        for (x = 0; x < width; x++) {
            printf("%4i ", *(mem++));
        }
        printf("\n");
    }
    printf("         ");
    i = matrix;
    for (x = 0; x <= seq2len - 1; x++) {
        printf("%4c ", seq2[x]);
    }
    for (y = 0; y <= seq1len; y++) {
        i = matrix + y * width;
        printf("\n");
        if (y == 0) printf("     ");
        else printf("%4c ", seq1[y-1]);
        for (x = 0; x <= seq2len; x++, i += width + 1) {
            printf("%4i ", *i);
        }
    }
    printf("\n");
#endif
}


/*
Return the smallest possible sequence distance, along with the ending position inside seq1.

Oprates on an alignment matrix as created with _sse2_make_matrix().
*/
static PyObject *_sse2_find_min(unsigned char *mem, const unsigned int seq1len, const unsigned int seq2len) {
    const unsigned int width = (seq2len+31) & ~0x0F;
    unsigned char *matrix = (unsigned char*)(((Py_uintptr_t)mem + 15) & ~(Py_uintptr_t)0x0F);
    unsigned int distance = seq2len;
    unsigned int position = 0;
    unsigned int i;

    unsigned char *cell = matrix + seq2len * width + seq2len;
    for (i = 0; i <= seq1len; i++, cell += width) {
        if (*cell < distance) {
            distance = *cell;
            position = i;
        }
    }

    return Py_BuildValue("II", distance, position);
}


/*
Align two sequences, finding one in the other.

Input arguments are two Python strings and a Python integer.
The second string will be searched inside the first.
The integer specifies the indel penalty score.

Return two python integers: the sequence distance, and the position of the match.
*/
static PyObject *sse2_align(PyObject *self, PyObject *args) {
    // Parse Python arguments.
    const char *seq1, *seq2;
    const unsigned int seq1len, seq2len;
    const unsigned char indel_score;
    unsigned char *matrix;
    PyObject *result;
    if (!PyArg_ParseTuple(args, "s#s#b", &seq1, &seq1len, &seq2, &seq2len, &indel_score)) {
        return NULL;
    }

    // This algorithm is designed to find seq2 in seq1.
    // Give the most pessimistic possible score if seq1 is shorter.
    if (seq1len < seq2len) {
        return Py_BuildValue("II", seq2len, 0);
    }

    matrix = _sse2_make_matrix(seq1len, seq2len, indel_score);
    _sse2_align(matrix, seq1len, seq2len, seq1, seq2, indel_score);
    result = _sse2_find_min(matrix, seq1len, seq2len);
    free(matrix);
    return result;
}


#endif

/**************************************************************************************************
  Plain Implementation
**************************************************************************************************/

/*
Allocates and returns a pointer to the alignment matrix.
*/
static unsigned char *_make_matrix(const unsigned int rows, const unsigned int columns,
                                   const unsigned char indel_score) {
    unsigned char *matrix = malloc(rows * columns * sizeof(char));
    unsigned char score = 0;
    unsigned int i;
    for (i = 1; i < rows; i++) {
        *(matrix + i*columns) = 0;
    }
    for (i = 0; i < columns; i++) {
        *(matrix + i) = score;
        score = _sadd(score, indel_score);
    }
    return matrix;
}


/*
Fill the alignment matrix as created with _make_matrix().
*/
static void _align(unsigned char *matrix, const unsigned int rows, const unsigned int columns,
                   const char *seq1, const char *seq2, const unsigned char indel_score) {
    unsigned int r, c;
    unsigned char *d = matrix,
                  *l = d + 1,
                  *u = d + columns,
                  *i = l + columns;
    for (r = 1; r < rows; r++, d++, l++, u++, i++) {
        for (c = 1; c < columns; c++, d++, l++, u++, i++) {
            *i = _min(
                _sadd(_min(*l, *u), indel_score),
                _sadd(*d, seq1[r - 1] != seq2[c - 1]));
        }
    }

#ifdef DEBUG
    // Print matrix.
    printf("         ");
    for (c = 0; c < columns - 1; c++) {
        printf("%4c ", seq2[c]);
    }
    for (r = 0; r < rows; r++) {
        printf("\n");
        if (r == 0) printf("     ");
        else printf("%4c ", seq1[r - 1]);
        for (c = 0; c < columns; c++) {
            printf("%4i ", *(matrix + r*columns + c));
        }
    }
    printf("\n");
#endif

}


/*
Return the smallest possible sequence distance, along with the ending position inside seq1.

Oprates on an alignment matrix as created with _make_matrix().
*/
static PyObject *_find_min(unsigned char *matrix, const unsigned int rows, const unsigned int columns) {
    unsigned int distance = columns - 1;
    unsigned int position = 0;
    unsigned int r;

    for (r = 0; r < rows; r++) {
        if (*(matrix + r*columns + columns-1) < distance) {
            distance = *(matrix + r*columns + columns-1);
            position = r;
        }
    }

    return Py_BuildValue("II", distance, position);
}


/*
Align two sequences, finding one in the other.

Input arguments are two Python strings and a Python integer.
The second string will be searched inside the first.
The integer specifies the indel penalty score.

Return two python integers: the sequence distance, and the position of the match.
*/
static PyObject *align(PyObject *self, PyObject *args) {
    // Parse Python arguments.
    const char *seq1, *seq2;
    const unsigned int seq1len, seq2len;
    const unsigned char indel_score;
    unsigned int rows, columns;
    unsigned char *matrix;
    PyObject *result;
    if (!PyArg_ParseTuple(args, "s#s#b", &seq1, &seq1len, &seq2, &seq2len, &indel_score)) {
        return NULL;
    }

    // This algorithm is designed to find seq2 in seq1.
    // Give the most pessimistic possible score if seq1 is shorter.
    if (seq1len < seq2len) {
        return Py_BuildValue("II", seq2len, 0);
    }

    // Allocate memory for the matrix.
    rows = seq1len + 1;
    columns = seq2len + 1;
    matrix = _make_matrix(rows, columns, indel_score);

    // Perform the alignment.
    _align(matrix, rows, columns, seq1, seq2, indel_score);
    result = _find_min(matrix, rows, columns);

    // Free the matrix memory and return the result.
    free(matrix);
    return result;
}


/*************************************************************************************************/


static PyMethodDef AlignerMethods[] = {
    {"align",  align, METH_VARARGS, METHOD_DOC},
    {"_simple_align",  align, METH_VARARGS, METHOD_DOC},
#if defined(_MSC_VER) || defined(__SSE2__)
    {"_sse2_align",  sse2_align, METH_VARARGS, METHOD_DOC},
#endif
    {NULL, NULL, 0, NULL}  /* Sentinel */
};


PyMODINIT_FUNC initsg_align(void) {
#if defined(_MSC_VER) || defined(__SSE2__)
    int cpuInfo[4];
    _cpuid(cpuInfo, 1);
#ifdef DEBUG
    printf("CPUID: %02x %02x %02x %02x\n", cpuInfo[0], cpuInfo[1], cpuInfo[2], cpuInfo[3]);
#endif
    if (cpuInfo[3] & (1 << 26)) {
        AlignerMethods[0].ml_meth = sse2_align;
    }
#endif

    (void) Py_InitModule("sg_align", AlignerMethods);
}
