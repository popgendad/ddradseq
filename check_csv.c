/* file: check_csv.c
 * description: Check the integrity of the input CSV file
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
#include "ddradseq.h"

static int compare(const void *a, const void *b);

int check_csv(const CMD *cp)
{
	char buf[MAX_LINE_LENGTH];
	char **list = NULL;
	char *p = NULL;
	char *q = NULL;
	unsigned int i = 0;
	unsigned int nlines = 0;
	ptrdiff_t diffp;
	ptrdiff_t diffq;
	FILE *lf = cp->lf;
	gzFile in;

	in = gzopen(cp->csvfile, "rb");
	if (!in)
	{
		logerror(lf, "%s:%d Could not read CSV database file %s into memory.\n",
			     __func__, __LINE__, cp->csvfile);
		return 1;;
	}

	/* Read CSV and count the number of entries */
	while (gzgets(in, buf, MAX_LINE_LENGTH) != Z_NULL)
		nlines++;

	/* Allocate memory for list from heap */
	list = malloc(nlines * sizeof(char*));
	if (UNLIKELY(!list))
	{
		logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 1;
	}
	for (i = 0; i < nlines; i++)
	{
		list[i] = malloc(MAX_LINE_LENGTH);
		if (UNLIKELY(!list[i]))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
	}

	/* Rewind input stream */
	gzrewind(in);

	/* Populate the list of lines */
	for (i = 0; i < nlines; i++)
		gzgets(in, list[i], MAX_LINE_LENGTH);

	/* Sort the list */
	qsort(list, nlines, sizeof(char*), compare);

	/* Iterate through CSV lines */
	for (i = 1; i < nlines; i++)
	{
		/* Look for duplicate lines */
		if (string_equal(list[i], list[i-1]))
		{
			logerror(lf, "%s:%d: CSV database file contains identical lines.\n",
			         __func__, __LINE__);
			return 1;
		}

		/* Delimit last field */
		p = strrchr(list[i-1], ',');
		q = strrchr(list[i], ',');
		diffp = p - list[i-1];
		diffq = q - list[i];

		/* If substrings are the same length */
		if (diffp == diffq)
		{
			if (strncmp(list[i-1], list[i], diffp) == 0)
			{
				if (string_equal(p+1, q+1) == 0)
				{
					logerror(lf, "%s:%d: Different individuals have the same key "
					         "pattern.\n", __func__, __LINE__);
					return 1;
				}
			}
		}
	}

	/* Close input file stream */
	gzclose(in);

	/* Deallocate heap memory */
	for (i = 0; i < nlines; i++)
		free(list[i]);
	free(list);

	return 0;
}

static int compare(const void *a, const void *b)
{
	return strcmp(*(char **) a, *(char **) b);
}