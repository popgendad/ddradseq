/* file: parse_main.c
 * description: Entry point for the parse modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

int parse_main(const CMD *cp)
{
	char **filelist = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int nfiles = 0;
	khash_t(pool_hash) *h = NULL;
	khash_t(mates) *m = NULL;
	FILE *lf = cp->lf;

	/* Check the integrity of the CSV input database file */
	ret = check_csv(cp);
	if (ret)
	{
		logerror(lf, "%s:%d Problem with the format of the CSV database file.\n",
			     __func__, __LINE__);
		return 1;
	}

	/* Read CSV database into memory */
	h = read_csv(cp);
	if (!h)
	{
		logerror(lf, "%s:%d Failed to read CSV database into memory.\n",
			     __func__, __LINE__);
		return 1;
	}

	/* Check for write permissions on parent of output directory */
	ret = create_dirtree(cp, h);
	if (ret)
		return 1;

	/* Initialize hash for mate pair information */
	m = kh_init(mates);
	if (!m)
		return 1;

	/* Get list of all files */
	nfiles = traverse_dirtree(cp, __func__, &filelist);
	if (!filelist)
		return 1;

	if (nfiles < 1)
	{
		logerror(lf, "%s:%d No input fastQ files found.\n", __func__, __LINE__);
		return 1;
	}

	for (i = 0; i < nfiles; i += 2)
	{
		char *ffor = NULL;
		char *frev = NULL;
		size_t spn = 0;

		/* Construct output file names */
		ffor = strdup(filelist[i]);
		if (UNLIKELY(!ffor))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
		frev = strdup(filelist[i+1]);
		if (UNLIKELY(!frev))
		{
			logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
			return 1;
		}
		spn = strcspn(ffor, ".");
		ret = strncmp(ffor, frev, spn);
		if (ret)
		{
			logerror(lf, "%s:%d Files \'%s\' and \'%s\' do not appear to be mate-"
			         "pairs.\n", __func__, __LINE__, ffor, frev);
			return 1;
		}

		/* Print informational update to log file */
		loginfo(lf, "Deciphering mate-pair information for \'%s\' and \'%s\'.\n", ffor, frev);

		/* Read the forward fastQ input file */
		ret = parse_fastq(cp, FORWARD, ffor, h, m);
		if (ret)
			return 1;

		/* Read the reverse fastQ input file */
		ret = parse_fastq(cp, REVERSE, frev, h, m);
		if (ret)
			return 1;
		free(ffor);
		free(frev);
	}

	/* Deallocate memory from the heap */
	for (i = 0; i < nfiles; i++)
		free(filelist[i]);
	free(filelist);
	free_db(h);
	free_matedb(m);

	/* Print informational message to log */
	loginfo(lf, "Parse step of pipeline is complete.\n");

	return 0;
}
