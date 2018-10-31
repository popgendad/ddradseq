/* file: parse_fastq.c
 * description: Parses a fastQ file by index sequence
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <errno.h>
#include "khash.h"
#include "ddradseq.h"

extern int errno;

int parse_fastq(const CMD *cp, const int orient, const char *filename, khash_t(pool_hash) *h,
                khash_t(mates) *m)
{
	char *r = NULL;
	char *q = NULL;
	char *errstr = NULL;
	char buffer[BUFLEN];
	int ret = 0;
	size_t numlines = 0;
	size_t bytes_read = 0;
	size_t buff_rem = 0;
	khint_t i = 0;
	khint_t j = 0;
	khint_t k = 0;
	khash_t(barcode) *b = NULL;
	khash_t(pool) *p = NULL;
	BARCODE *bc = NULL;
	POOL *pl = NULL;
	FILE *lf = cp->lf;
	gzFile fin;

	/* Print informational message to log */
	loginfo(lf, "Parsing fastQ file \'%s\'.\n", filename);

	/* Open input file */
	fin = gzopen(filename, "rb");
	if (!fin)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Unable to open file \'%s\': %s.\n", __func__, __LINE__,
		         filename, errstr);
		return 1;
	}

	/* Initialize buffer */
	memset(buffer, 0, sizeof(buffer));
	r = &buffer[0];

	/* Iterate through blocks from input fastQ file */
	while (1)
	{
		/* Read block from file into input buffer */
		bytes_read = gzread(fin, &buffer[buff_rem], BUFLEN - buff_rem - 1);
		if (bytes_read < 0)
		{
			logerror(lf, "%s:%d Failed to read data from file \'%s\': %s.\n",
			         __func__, __LINE__, filename);
			return 1;
		}
		/* Set null terminating character on input buffer */
		buffer[bytes_read + buff_rem] = '\0';

		/* Count lines in input buffer */
		q = &buffer[0];
		numlines = count_lines(q);
		r = clean_buffer(q, &numlines);
		if (orient == FORWARD)
			ret = parse_forwardbuffer(cp, q, numlines, h, m);
		else
			ret = parse_reversebuffer(cp, q, numlines, h, m);
		if (ret)
			return 1;
		buff_rem = reset_buffer(q, r);

		/* Check if we are at the end of file */
		if (gzeof(fin))
			break;
	}

	/* Flush remaining data in buffers */
	for (i = kh_begin(h); i != kh_end(h); i++)
	{
		if (kh_exist(h, i))
		{
			p = kh_value(h, i);
			for (j = kh_begin(p); j != kh_end(p); j++)
			{
				if (kh_exist(p, j))
				{
					pl = kh_value(p, j);
					b = pl->b;
					for (k = kh_begin(b); k != kh_end(b); k++)
					{
						if (kh_exist(b, k))
						{
							bc = kh_value(b, k);
							if (bc->curr_bytes > 0)
							{
								ret = flush_buffer(orient, bc, lf);
								if (ret)
								{
									logerror(lf, "%s:%d Problem writing buffer to file.\n",
									         __func__, __LINE__);
									return 1;
								}
							}
						}
					}
				}
			}
		}
	}

	/* Close input file */
	gzclose(fin);

	/* Print informational message to log */
	loginfo(lf, "Successfully parsed fastQ file \'%s\'.\n", filename);

	return 0;
}
