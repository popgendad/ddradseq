/* file pair_mates.c
 * description: Pair mates in two fastQ files
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#include "ddradseq.h"
#include "khash.h"

extern int errno;

int pair_mates(const char *filename, const khash_t(fastq) *h, const char *ffor,
               const char *frev, FILE *lf)
{
	char **buf = NULL;
	char *idline = NULL;
	char *mkey = NULL;
	char *pstart = NULL;
	char *pend = NULL;
	char *errstr = NULL;
	int i = 0;
	size_t l = 0;
	size_t lc = 0;
	size_t pos = 0;
	size_t strl = 0;
	ptrdiff_t plen = 0;
	khint_t k = 0;
	gzFile in;
	gzFile fout;
	gzFile rout;
	FASTQ *e = NULL;

	/* Allocate memory for buffer from heap */
	buf = malloc(BSIZE * sizeof(char*));
	if (UNLIKELY(!buf))
	{
		logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 1;
	}
	for (i = 0; i < BSIZE; i++)
	{
		buf[i] = malloc(MAX_LINE_LENGTH);
		if (UNLIKELY(!buf[i]))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
	}

	/* Open the fastQ input stream */
	in = gzopen(filename, "rb");
	if (!in)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Unable to open input file \'%s\': %s.\n", __func__, __LINE__,
		         filename, errstr);
		return 1;
	}

	/* Open the output fastQ file streams */
	fout = gzopen(ffor, "wb");
	if (!fout)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Unable to forward output file \'%s\': %s.\n", __func__,
		         __LINE__, ffor, errstr);
		return 1;
	}

	rout = gzopen(frev, "wb");
	if (!rout)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Unable to reverse output file \'%s\': %s.\n", __func__,
		         __LINE__, frev, errstr);
		return 1;
	}

	/* Enter data from the fastQ input file into the database */
	while (1)
	{
		/* Fill up the buffer */
		for (lc = 0; lc < BSIZE; lc++)
		{
			/* Clear the buffer */
			memset(buf[lc], 0, MAX_LINE_LENGTH);

			/* Get line from the fastQ input stream */
			if (gzgets(in, buf[lc], MAX_LINE_LENGTH) == Z_NULL)
				break;
		}

		/* Iterate through lines in the buffer */
		for (l = 0; l < lc; l++)
		{
			/* We are at the end of one fastQ entry */
			if (l % 4 == 3)
			{
				/* Parse entry identifier */
				pos = strcspn(buf[l-3], "\n");
				buf[l-3][pos] = '\0';
				strl = strlen(&buf[l-3][1]);

				/* Construct fastQ hash key */
				idline = malloc(strl + 1u);
				if (UNLIKELY(!idline))
				{
					logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
					return 1;
				}
				strcpy(idline, &buf[l-3][1]);

				/* Parse Illumina identifier line and construct hash key */
				pstart = strchr(idline, ':');
				pend = strchr(idline, ' ');
				plen = pend - pstart;
				mkey = strndup(pstart+1, plen - 1);
				if (UNLIKELY(!mkey))
				{
					logerror(lf, "%s:%d fastQ header parsing error.\n", __func__, __LINE__);
					return 1;
				}
				k = kh_get(fastq, h, mkey);
				if (k != kh_end(h))
					e = kh_value(h, k);
				free(mkey);
				free(idline);

				if (e != NULL)
				{
					/* Parse DNA sequence */
					pos = strcspn(buf[l-2], "\n");
					buf[l-2][pos] = '\0';

					/* Parse quality sequence */
					pos = strcspn(buf[l], "\n");
					buf[l][pos] = '\0';

					/* Need to construct output file streams */
					gzprintf(fout, "@%s\n%s\n+\n%s\n", e->id, e->seq, e->qual);
					gzprintf(rout, "@%s\n%s\n+\n%s\n", &buf[l-3][1], &buf[l-2][0], &buf[l][0]);
				}
			}
		}

		/* If we are at the end of the file */
		if (lc < BSIZE)
			break;
	}

	/* Free memory for buffer to heap */
	for (i = 0; i < BSIZE; i++)
		free(buf[i]);
	free(buf);

	/* Close input stream */
	gzclose(in);
	gzclose(fout);
	gzclose(rout);

	return 0;
}
