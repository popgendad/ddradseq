/* file: local_align.c
 * description: Calculates the local sequence alignment by Smith-Waterman algorithm
 * author: Adapted from https://github.com/attractivechaos/klib/blob/master/ksw.c
 *         by Heng Li by Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include "ddradseq.h"

/* Define alphabet size */
#define ALPHA_SIZE 5

/* Globally scoped variables */
const ALIGN_RESULT g_defr = { 0, -1, -1, -1, -1, -1, -1 };

/* Function prototypes */
static ALIGN_QUERY* align_init(int qlen, const char *query, const char *mat, FILE *lf);

ALIGN_RESULT smith_waterman(ALIGN_QUERY *q, int tlen, const char *target, int _gapo,
                            int _gape, int xtra, FILE *lf);

static void revseq(int l, char *s);


ALIGN_RESULT local_align(int qlen, char *query, int tlen, char *target,
                         const char *mat, int gapo, int gape, int xtra, FILE *lf)
{
	ALIGN_QUERY *q;
	ALIGN_RESULT r;
	ALIGN_RESULT rr;

	q = align_init(qlen, query, mat, lf);
	r = smith_waterman(q, tlen, target, gapo, gape, xtra, lf);
	free(q);
	if (((xtra & KSW_XSTART) == 0 || (xtra & KSW_XSUBO)) && r.score < (xtra & 0xffff))
		return r;
	revseq(r.query_end + 1, query);

	/* +1 because qe/te points to the exact end */
	/* not the position after the end */
	revseq(r.target_end + 1, target);
	q = align_init(r.query_end + 1, query, mat, lf);
	rr = smith_waterman(q, tlen, target, gapo, gape, KSW_XSTOP | r.score, lf);
	revseq(r.query_end + 1, query);
	revseq(r.target_end + 1, target);
	free(q);
	if (r.score == rr.score)
	{
		r.target_begin = r.target_end - rr.target_end;
		r.query_begin = r.query_end - rr.query_end;
	}
	return r;
}

static ALIGN_QUERY *align_init(int qlen, const char *query, const char *mat, FILE *lf)
{
	int slen = 0;
	int a = 0;
	int tmp = 0;
	int p = 0;
	ALIGN_QUERY *q = NULL;

	/* Number of values per __m128i */
	p = 16;

	/* Segmented length */
	slen = (qlen + p - 1) / p;

	/* Allocate memory for query profile */
	if ((q = malloc(sizeof(ALIGN_QUERY) + 256 + 16 * slen * (ALPHA_SIZE + 4))) == NULL)
	{
		logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return NULL;
	}

	/* Align memory */
	q->qp = (__m128i*)(((size_t)q + sizeof(ALIGN_QUERY) + 15u) >> 4 << 4);
	q->H0 = q->qp + slen * ALPHA_SIZE;
	q->H1 = q->H0 + slen;
	q->E  = q->H1 + slen;
	q->Hmax = q->E + slen;
	q->slen = slen;
	q->qlen = qlen;

	/* Compute shift */
	tmp = ALPHA_SIZE * ALPHA_SIZE;

	/* Find the minimum and maximum score */
	for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; a++)
	{
		if (mat[a] < (char)q->shift)
			q->shift = mat[a];
		if (mat[a] > (char)q->mdiff)
			q->mdiff = mat[a];
	}
	q->max = q->mdiff;

	/* NB: q->shift is uint8_t */
	q->shift = 256 - q->shift;

	/* Difference between the min and max scores */
	q->mdiff += q->shift;

	char *t = (char*)q->qp;
	for (a = 0; a < ALPHA_SIZE; a++)
	{
		int i = 0;
		int k = 0;
		int nlen = slen * p;
		const char *ma = mat + a * ALPHA_SIZE;

		for (i = 0; i < slen; i++)
			for (k = i; k < nlen; k += slen)
				*t++ = (k >= qlen ? 0 : ma[(unsigned char)query[k]]) + q->shift;
	}
	return q;
}

ALIGN_RESULT smith_waterman(ALIGN_QUERY *q, int tlen, const char *target,
                            int _gapo, int _gape, int xtra, FILE *lf)
{
	int slen = 0;
	int i = 0;
	int m_b = 0;
	int n_b = 0;
	int te = -1;
	int gmax = 0;
	int minsc = 0;
	int endsc = 0;
	uint64_t *b = NULL;
	__m128i zero;
	__m128i gapoe;
	__m128i gape;
	__m128i shift;
	__m128i *H0;
	__m128i *H1;
	__m128i *E;
	__m128i *Hmax;
	ALIGN_RESULT r;

#define __max_16(ret, xx) do { \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
		(xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
		(ret) = _mm_extract_epi16((xx), 0) & 0x00ff; \
	} while (0)

	/* Initialization */
	r = g_defr;
	minsc = xtra & KSW_XSUBO ? xtra & 0xffff : 0x10000;
	endsc = xtra & KSW_XSTOP ? xtra & 0xffff : 0x10000;
	m_b = n_b = 0;
	b = 0;
	zero = _mm_set1_epi32(0);
	gapoe = _mm_set1_epi8(_gapo + _gape);
	gape = _mm_set1_epi8(_gape);
	shift = _mm_set1_epi8(q->shift);
	H0 = q->H0;
	H1 = q->H1;
	E = q->E;
	Hmax = q->Hmax;
	slen = q->slen;
	for (i = 0; i < slen; i++)
	{
		_mm_store_si128(E + i, zero);
		_mm_store_si128(H0 + i, zero);
		_mm_store_si128(Hmax + i, zero);
	}

	/* Core loop */
	for (i = 0; i < tlen; i++)
	{
		int j = 0;
		int k = 0;
		int cmp = 0;
		int imax = 0;
		__m128i e;
		__m128i h;
		__m128i f = zero;
		__m128i max = zero;
		__m128i *S = q->qp + target[i] * slen;

		/* h={2,5,8,11,14,17,-1,-1} in the above example */
		h = _mm_load_si128(H0 + slen - 1);

		/* h=H(i-1,-1); << instead of >> because x64 is little-endian */
		h = _mm_slli_si128(h, 1);

		for (j = 0; LIKELY(j < slen); j++)
		{
			/* SW cells are computed in the following order:
			 *	 H(i,j)	  = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			 *	 E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
			 *	 F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
			 */

			/* Compute H'(i,j); note that at the beginning, h=H'(i-1,j-1) */
			h = _mm_adds_epu8(h, _mm_load_si128(S + j));

			/* h=H'(i-1,j-1)+S(i,j) */
			h = _mm_subs_epu8(h, shift);

			/* e=E'(i,j) */
			e = _mm_load_si128(E + j);
			h = _mm_max_epu8(h, e);

			/* h=H'(i,j) */
			h = _mm_max_epu8(h, f);

			/* Set max */
			max = _mm_max_epu8(max, h);

			/* Save to H'(i,j) */
			_mm_store_si128(H1 + j, h);

			/* Now compute E'(i+1,j) */
			/* h=H'(i,j)-gapo */
			h = _mm_subs_epu8(h, gapoe);

			/* e=E'(i,j)-gape */
			e = _mm_subs_epu8(e, gape);

			/* e=E'(i+1,j) */
			e = _mm_max_epu8(e, h);

			/* Save to E'(i+1,j) */
			_mm_store_si128(E + j, e);

			/* Now compute F'(i,j+1) */
			f = _mm_subs_epu8(f, gape);
			f = _mm_max_epu8(f, h);

			/* get H'(i-1,j) and prepare for the next j */
			/* h=H'(i-1,j) */
			h = _mm_load_si128(H0 + j);
		}

		/* NB: we do not need to set E(i,j) as we disallow */
		/* adjecent insertion and then deletion */
		/* this block mimics SWPS3; NB: H(i,j) updated in the */
		/* lazy-F loop cannot exceed max */
		for (k = 0; LIKELY(k < 16); k++)
		{
			f = _mm_slli_si128(f, 1);
			for (j = 0; LIKELY(j < slen); j++)
			{
				h = _mm_load_si128(H1 + j);

				/* h=H'(i,j) */
				h = _mm_max_epu8(h, f);
				_mm_store_si128(H1 + j, h);
				h = _mm_subs_epu8(h, gapoe);
				f = _mm_subs_epu8(f, gape);
				cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h), zero));
				if (UNLIKELY(cmp == 0xffff))
					goto end_loop16;
			}
		}
end_loop16:
		/* imax is the maximum number in max */
		__max_16(imax, max);

		/* Write the b array; this condition adds branching unfornately */
		if (imax >= minsc)
		{
			/* Then append */
			if (n_b == 0 || (int)b[n_b-1] + 1 != i)
			{
				if (n_b == m_b)
				{
					m_b = m_b ? m_b << 1 : 8;
					if ((b = realloc(b, 8 * m_b)) == NULL)
					{
						logerror(lf, "%s:%d Memory reallocation failure.\n", __func__, __LINE__);
						return r;
					}
				}
				b[n_b++] = (uint64_t)imax << 32 | i;
			}
			else if ((int)(b[n_b-1] >> 32) < imax)
			{
				/* Modify the last */
				b[n_b-1] = (uint64_t)imax << 32 | i;
			}
		}
		if (imax > gmax)
		{
			gmax = imax;

			/* te is the end position on the target */
			te = i;

			/* Keep the H1 vector */
			for (j = 0; LIKELY(j < slen); j++)
				_mm_store_si128 (Hmax + j, _mm_load_si128 (H1 + j));

			if (gmax + q->shift >= 255 || gmax >= endsc)
				break;
		}
		S = H1;
		H1 = H0;

		/* Swap H0 and H1 */
		H0 = S;
	}

	r.score = gmax + q->shift < 255 ? gmax : 255;
	r.target_end = te;

	/* Get a->query_end, the end of query match */
	/* find the 2nd best score */
	if (r.score != 255)
	{
		int max = -1;
		int low = 0;
		int high = 0;
		int qlen = slen * 16;
		char *t = (char*)Hmax;

		for (i = 0; i < qlen; i++, t++)
		{
			if ((int)*t > max)
			{
				max = *t;
				r.query_end = i / 16 + i % 16 * slen;
			}
		}
		if (b)
		{
			i = (r.score + q->max - 1) / q->max;
			low = te - i;
			high = te + i;
			for (i = 0; i < n_b; i++)
			{
				int e = (int)b[i];
				if ((e < low || e > high) && (int)(b[i] >> 32) > r.score2)
				{
					r.score2 = b[i] >> 32;
					r.target_end2 = e;
				}
			}
		}
	}
	free (b);

	return r;
}

static void revseq(int l, char *s)
{
	int i = 0;
	int t = 0;

	for (i = 0; i < l >> 1; i++)
	{
		t = s[i];
		s[i] = s[l - 1 - i];
		s[l - 1 - i] = t;
	}
}
