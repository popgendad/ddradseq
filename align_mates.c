/* file align_mates.c
 * description: Align mates in two fastQ files and trim 3' end of reverse sequences
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include "ddradseq.h"

#define NBASES 4

/*Globally scoped variables */
const char seq_nt4_table[256] = {
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
const char alpha[5] = "ACGTN";
extern int errno;

int align_mates(const CMD *cp, const char *forin, const char *revin, const char *forout, const char *revout)
{
	char **fbuf = NULL;
	char **rbuf = NULL;
	char *errstr = NULL;
	char mat[25];
	int i = 0;
	int j = 0;
	int k = 0;
	int xtra = KSW_XSTART;
	const int sa = 1;
	const int sb = 3;
	const int gap_open = cp->gapo;
	const int gap_extend = cp->gape;
	const int min_score = cp->score;
	unsigned int count = 0;
	size_t l = 0;
	size_t lc = 0;
	FILE *lf = cp->lf;
	gzFile fin;
	gzFile rin;
	gzFile fout;
	gzFile rout;

	/* Allocate buffer memory from the heap */
	fbuf = malloc(BSIZE * sizeof(char*));
	if (UNLIKELY(fbuf == NULL))
	{
		logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 1;
	}
	rbuf = malloc(BSIZE * sizeof(char*));
	if (UNLIKELY(rbuf == NULL))
	{
		logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 1;
	}
	for (i = 0; i < BSIZE; i++)
	{
		fbuf[i] = malloc(MAX_LINE_LENGTH);
		if (UNLIKELY(!fbuf[i]))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
		rbuf[i] = malloc(MAX_LINE_LENGTH);
		if (UNLIKELY(!rbuf[i]))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
	}

	/* Open input forward fastQ file stream */
	fin = gzopen(forin, "rb");
	if (!fin)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Failed to open input forward fastQ file \'%s\': %s.\n",
		         __func__, __LINE__, forin, errstr);
		return 1;
	}

	/* Open input reverse fastQ file stream */
	rin = gzopen(revin, "rb");
	if (!rin)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Failed to open input reverse fastQ file \'%s\': %s.\n",
		         __func__, __LINE__, revin, errstr);
		return 1;
	}

	/* Open output forward fastQ file stream */
	fout = gzopen(forout, "wb");
	if (!fout)
	{
		logerror(lf, "%s:%d Failed to open forward output fastQ file \'%s\': %s.\n",
		         __func__, __LINE__, forout, errstr);
		return 1;
	}

	/* Open output reverse fastQ file stream */
	rout = gzopen(revout, "wb");
	if (!rout)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Failed to open reverse output fastQ file \'%s\': %s.\n",
		         __func__, __LINE__, revout, errstr);
		return 1;
	}

	/* Initialize the scoring matrix */
	for (i = k = 0; i < NBASES; i++)
	{
		for (j = 0; j < NBASES; j++)
			mat[k++] = (i == j) ? sa : -sb;

		/* Ambiguous base */
		mat[k++] = 0;
	}
	for (j = 0; j <= NBASES; j++)
		mat[k++] = 0;

	while (1)
	{
		/* Fill up the forward buffer */
		for (lc = 0; lc < BSIZE; lc++)
		{
			/* Clear the buffer */
			memset(fbuf[lc], 0, MAX_LINE_LENGTH);

			/* Get line from the fastQ input stream */
			if (gzgets(fin, fbuf[lc], MAX_LINE_LENGTH) == Z_NULL)
				break;
		}

		/* Fill up the reverse buffer */
		for (lc = 0; lc < BSIZE; lc++)
		{
			/* Clear the buffer */
			memset(rbuf[lc], 0, MAX_LINE_LENGTH);

			/* Get line from the fastQ input stream */
			if (gzgets(rin, rbuf[lc], MAX_LINE_LENGTH) == Z_NULL)
				break;
		}

		/* Iterate through lines in the buffers */
		for (l = 0; l < lc; l++)
		{
			/* We are at the end of one fastQ entry */
			if (l % 4 == 3)
			{
				ALIGN_RESULT r;
				char *target = NULL;
				char *query = NULL;
				target = strdup(&fbuf[l-2][0]);
				query = revcom(&rbuf[l-2][0], lf);
				if (!target || !query)
					return 1;
				int tlen = (int)strlen(target);
				target[tlen] = '\0';
				int qlen = (int)strlen(query);

				/* Transform sequences */
				for (i = 0; i < qlen; i++)
					query[i] = seq_nt4_table[(unsigned char)query[i]];
				for (i = 0; i < tlen; i++)
					target[i] = seq_nt4_table[(unsigned char)target[i]];

				/* Do the alignment */
				r = local_align(qlen, query, tlen, target, mat, gap_open, gap_extend, xtra, lf);
				free(target);
				free(query);

				/* Actually trim the sequence */
				if (r.score >= min_score)
				{
					/* Test trimming criterion */
					if (r.target_begin == 0 && r.query_begin > 0)
					{
						int new_end_pos = qlen - r.query_begin;
						char *seq = &rbuf[l-2][0];
						char *qual = &rbuf[l][0];
						seq[new_end_pos] = '\n';
						seq[new_end_pos+1] = '\0';
						qual[new_end_pos] = '\n';
						qual[new_end_pos+1] = '\0';
						count++;
					}
				}

				/* Write sequences to file */
				gzprintf(fout, "%s%s+\n%s", &fbuf[l-3][0], &fbuf[l-2][0], &fbuf[l][0]);
				gzprintf(rout, "%s%s+\n%s", &rbuf[l-3][0], &rbuf[l-2][0], &rbuf[l][0]);
			}
		}

		/* If we are at the end of the file */
		if (lc < BSIZE)
			break;
	}

	/* Print informational message to logfile */
	loginfo(lf, "%u sequences trimmed.\n", count);

	/* Free memory from the heap */
	for (i = 0; i < BSIZE; i++)
	{
		free(fbuf[i]);
		free(rbuf[i]);
	}
	free(fbuf);
	free(rbuf);

	/* Close all file streams */
	gzclose(fin);
	gzclose(rin);
	gzclose(fout);
	gzclose(rout);

	return 0;
}
