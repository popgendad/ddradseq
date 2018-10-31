/* file: read_csv.c
 * description: Read the input CSV file into hash database
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <zlib.h>
#include "khash.h"
#include "ddradseq.h"

khash_t(pool_hash) *read_csv(const CMD *cp)
{
	const char *csvfile = cp->csvfile;  /* Pointer to CSV database file name */
	const char *outpath = cp->outdir;	/* Pointer to parent of output directories */
	char buf[MAX_LINE_LENGTH];          /* File input buffer */
	char seps[] = ",";			        /* CSV entry separator character */
	char *tok = NULL;			        /* Holds parsed CSV tokens */
	char *r = NULL;				        /* Residual pointer for strtok_r */
	char *tmp = NULL;			        /* Temporary pointer */
	bool trail = false;                 /* Boolean indicator of trailing slash */
	int a = 0;					        /* Return value for database entry */
	size_t strl = 0;			        /* Generic string length holder */
	size_t pathl = 0;			        /* Length of path string */
	gzFile in;					        /* Input file stream */
	khint_t i = 0;                      /* Generic hash iterator */
	khint_t j = 0;                      /* Generic hash iterator */
	khint_t k = 0;                      /* Generic hash iterator */
	khash_t(barcode) *b = NULL;         /* Pointer to barcode hash table */
	khash_t(pool) *p = NULL;            /* Pointer to pool hash table */
	khash_t(pool_hash) *h = NULL;       /* Pointer to flow hash table */
	BARCODE *bc = NULL;                 /* Pointer barcode data structure */
	POOL *pl = NULL;                    /* Pointer to pool data structure */
	FILE *lf = cp->lf;                  /* Pointer to log file stream */

	/* Print informational message to log */
	loginfo(lf, "Parsing CSV database file \'%s\'.\n", csvfile);

	/* Check for trailing slash on outpath */
	pathl = strlen(outpath);
	if (outpath[pathl - 1u] == '/')
		trail = true;

	/* Open input database text file stream */
	in = gzopen(csvfile, "rb");
	if (!in)
	{
		logerror(lf, "%s:%d Could not read CSV database file %s into memory.\n",
			     __func__, __LINE__, csvfile);
		return NULL;
	}

	/* Initialize top-level hash */
	h = kh_init(pool_hash);

	/* Read CSV and populate the database */
	while (gzgets(in, buf, MAX_LINE_LENGTH) != Z_NULL)
	{
		/* Re-initialize measure of outfile path string length */
		pathl = strlen(outpath);

		/* Get the flowcell entry */
		if ((tok = strtok_r(buf, seps, &r)) == NULL)
		{
			logerror(lf, "%s:%d Parsing CSV file failed.\n", __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for flowcell string */
		strl = strlen(tok);
		tmp = malloc(strl + 1u);
		if (UNLIKELY(!tmp))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return NULL;
		}
		if (!cp->across)
			pathl += strl + 1u;
		strcpy(tmp, tok);

		/* Put flowcell string as key in top-level hash */
		i = kh_put(pool_hash, h, tmp, &a);

		/* If this flowcell is a new entry-- */
		/* initialize a second-level hash and add to value */
		if (a)
		{
			p = kh_init(pool);
			kh_value(h, i) = p;
		}
		else
			free(tmp);

		/* Get pool sequence */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			logerror(lf, "%s:%d Parsing CSV file failed.\n", __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for pool sequence string */
		strl = strlen(tok);
		tmp = malloc(strl + 1u);
		if (UNLIKELY(!tmp))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return NULL;
		}
		strcpy(tmp, tok);

		/* Put pool sequence string in second-level hash */
		p = kh_value(h, i);
		j = kh_put(pool, p, tmp, &a);

		/* If this pool sequence is a new entry-- */
		/* initialize a third-level hash and POOL data structure */
		/* Add the new hash to POOL and add POOL to second-level hash */
		if (a)
		{
			pl = malloc(sizeof(POOL));
			if (UNLIKELY(!pl))
			{
				logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
				return NULL;
			}
			b = kh_init(barcode);
			pl->b = b;
			kh_value(p, j) = pl;
		}
		else
			free(tmp);

		/* Get pool value */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			logerror(lf, "%s:%d Parsing CSV file failed.\n", __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for pool identifier */
		strl = strlen(tok);
		tmp = malloc(strl + 1u);
		if (UNLIKELY(!tmp))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return NULL;
		}
		pathl += strl + 1u;
		strcpy(tmp, tok);

		/* Put pool identifier string and directory path */
		/* in POOL data structure */
		pl = kh_value(p, j);
		if (a)
			pl->poolID = tmp;
		else
			free(tmp);
		tmp = malloc(pathl + 1u);
		if (UNLIKELY(!tmp))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return NULL;
		}
		if (trail)
		{
			if (cp->across)
				sprintf(tmp, "%s%s", outpath, pl->poolID);
			else
				sprintf(tmp, "%s%s/%s", outpath, kh_key(h, i), pl->poolID);
		}
		else
		{
			if (cp->across)
				sprintf(tmp, "%s/%s", outpath, pl->poolID);
			else
				sprintf(tmp, "%s/%s/%s", outpath, kh_key(h, i), pl->poolID);
		}
		if (a)
			pl->poolpath = tmp;
		else
			free(tmp);

		/* Get barcode sequence */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			logerror(lf, "%s:%d Parsing CSV file failed.\n", __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for barcode sequence string */
		strl = strlen(tok);
		tmp = malloc(strl + 1u);
		if (UNLIKELY(!tmp))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return NULL;
		}
		strcpy(tmp, tok);

		/* Put barcode sequence in third-level hash */
		b = pl->b;
		k = kh_put(barcode, b, tmp, &a);
		if (b->size == 1)
			pl->barcode_length = strl;
		else
		{
			if (pl->barcode_length != strl)
			{
				logerror(lf, "%s:%d Unequal barcode lengths in CSV file %s.\n",
					     __func__, __LINE__, csvfile);
				return NULL;
			}
		}
		/* If this barcode is a new entry-- */
		/* initialize a fourth-level hash and add to value */
		if (a)
		{
			bc = malloc(sizeof(BARCODE));
			if (UNLIKELY(!bc))
			{
				logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
				return NULL;
			}
			bc->buffer = malloc(BUFLEN);
			if (UNLIKELY(!bc->buffer))
			{
				logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
				return NULL;
			}
			bc->buffer[0] = '\0';
			bc->curr_bytes = 0;
			kh_value(b, k) = bc;
		}
		else
			free(tmp);

		/* Get barcode value */
		if ((tok = strtok_r(NULL, seps, &r)) == NULL)
		{
			logerror(lf, "%s:%d Parsing CSV file failed.\n", __func__, __LINE__);
			return NULL;
		}

		/* Alloc memory for barcode value */
		strl = strcspn(tok, " \n");
		tmp = malloc(strl + 1u);
		if (UNLIKELY(!tmp))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return NULL;
		}
		pathl += strl + 1u;
		strncpy(tmp, tok, strl);
		tmp[strl] = '\0';

		/* Add barcode value to BARCODE data structure */
		bc = kh_value(b, k);
		bc->smplID = tmp;
		pathl += 21u;
		tmp = malloc(pathl + 1u);
		if (UNLIKELY(!tmp))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return NULL;
		}
		sprintf(tmp, "%s/parse/smpl_%s.R1.fq.gz", pl->poolpath, bc->smplID);
		bc->outfile = tmp;
	}

	/* Close input CSV file stream */
	gzclose(in);

	/* Print informational message to log */
	loginfo(lf, "Successfully parsed CSV database file \'%s\'.\n", csvfile);

	return h;
}
