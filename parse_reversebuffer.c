/* file: parse_reversebuffer.c
 * description: Parses reverse fastQ entries in the buffer
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include "khash.h"
#include "ddradseq.h"

int parse_reversebuffer(const CMD *cp, char *buff, const size_t nl, const khash_t(pool_hash) *h,
                        const khash_t(mates) *m)
{
	bool *skip = NULL;
	char *q = buff;
	char *copy = NULL;
	char *idline = NULL;
	char *mkey = NULL;
	char *pstart = NULL;
	char *pend = NULL;
	char *flowcell = NULL;
	char *barcode_sequence = NULL;
	char *index_sequence = NULL;
	char *dna_sequence = NULL;
	char *qual_sequence = NULL;
	int ret = 0;
	size_t add_bytes = 0;
	size_t l = 0;
	size_t ll = 0;
	ptrdiff_t plen = 0;
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khint_t mk = 0;
	khash_t(barcode) *b = NULL;
	khash_t(pool) *p = NULL;
	BARCODE *bc = NULL;
	POOL *pl = NULL;
	FILE *lf = cp->lf;

	/* Indicator variable whether to skip processing a line */
	skip = calloc(nl, sizeof(bool));
	if (!skip)
	{
		logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 1;
	}

	/* Iterate through lines in the buffer */
	for (l = 0; l < nl; l++)
	{
		ll = strlen(q);
		if (!skip[l])
		{
			switch (l % 4)
			{
				case 0:
					/* Make a copy of the Illumina identifier line */
					copy = strdup(q);
					if (!copy)
					{
						logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					idline = strdup(q);
					if (!idline)
					{
						logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}

					/* Parse Illumina identifier line and construct fastQ hash key */
					pstart = strchr(idline, ':');
					pend = strchr(idline, ' ');
					plen = pend - pstart;
					mkey = strndup(pstart+1, plen - 1);
					if (UNLIKELY(!mkey))
					{
						logerror(lf, "%s:%d fastQ header parsing error.\n", __func__, __LINE__);
						return 1;
					}
					/* Parse flowcell identifier */
					pstart = strchr(pstart+1, ':');
					pend = strchr(pstart+1, ':');
					plen = pend - pstart;
					flowcell = strndup(pstart+1, plen - 1);
					if (!flowcell)
					{
						logerror(lf, "%s:%d Illumina ID parsing failure.\n", __func__, __LINE__);
						return 1;
					}
					/* Parse index sequence */
					pstart = strrchr(idline, ':');
					index_sequence = strdup(pstart+1);
					if (!index_sequence)
					{
						logerror(lf, "%s:%d Illumina ID parsing failure.\n", __func__, __LINE__);
						return 1;
					}

					/* Lookup flow cell identifier */
					i = kh_get(pool_hash, h, flowcell);

					/* Flow cell identifier is not present in database */
					if (i == kh_end(h))
					{
						logwarn(lf, "Hash lookup failure using key %s.\n", flowcell);
						logwarn(lf, "Skipping sequence: %s\n", idline);
						skip[l+1] = true;
						skip[l+2] = true;
						skip[l+3] = true;
						free(idline);
						free(copy);
						free(flowcell);
						break;
					}
					else
						p = kh_value(h, i);

					/* Lookup pool identifier */
					j = kh_get(pool, p, index_sequence);
					pl = kh_value(p, j);
					b = pl->b;

					/* Retrieve barcode sequence of mate */
					mk = kh_get(mates, m, mkey);
					if (mk == kh_end(m))
					{
						logwarn(lf, "Hash lookup failure using key %s.\n", mkey);
						logwarn(lf, "Skipping sequence: %s\n", idline);
						skip[l+1] = true;
						skip[l+2] = true;
						skip[l+3] = true;
						free(mkey);
						free(idline);
						free(copy);
						break;
					}
					barcode_sequence = kh_value(m, mk);
					free(mkey);

					/* Get the barcode entry of read's mate */
					k = kh_get(barcode, b, barcode_sequence);
					if (k != kh_end(b))
						bc = kh_value(b, k);
					else
					{
						skip[l+1] = true;
						skip[l+2] = true;
						skip[l+3] = true;
						free(idline);
					}

					/* Free memory */
					free(flowcell);
					free(index_sequence);
					free(copy);
					break;
				case 1:
					/* Sequence line */
					dna_sequence = strdup(q);
					if (!dna_sequence)
					{
						logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					break;
				case 2:
					/* Quality identifier line */
					break;
				case 3:
					/* Quality sequence line */
					qual_sequence = strdup(q);
					if (!qual_sequence)
					{
						logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					add_bytes = strlen(idline) + strlen(dna_sequence) +
								strlen(qual_sequence) + 5u;
					char *t = NULL;
					if ((bc->curr_bytes + add_bytes) >= BUFLEN)
					{
						ret = flush_buffer(REVERSE, bc, lf);
						if (ret)
						{
							logerror(lf, "%s:%d Problem writing to file.\n", __func__, __LINE__);
							return 1;
						}
					}
					bc->curr_bytes += add_bytes;
					t = malloc(add_bytes + 1u);
					if (!t)
					{
						logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					sprintf(t, "%s\n%s\n+\n%s\n", idline, dna_sequence, qual_sequence);
					strcat(bc->buffer, t);

					/* Free alloc'd memory for fastQ entry */
					free(t);
					free(idline);
					free(dna_sequence);
					free(qual_sequence);
					break;
			}
		}
		q += ll + 1u;
	}

	/* Free memory from the heap */
	free(skip);

	return 0;
}
