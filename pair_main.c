/* file: pair_main.c
 * description: Entry point for the pair modality
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

int pair_main(const CMD *cp)
{
	char *pch = NULL;
	char **filelist = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int nfiles = 0;
	FILE *lf = cp->lf;

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
		khash_t(fastq) *h = NULL;
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
		pch = strstr(ffor, "parse");
		if (!pch)
			return 1;
		strncpy(pch, "pairs", DNAME_LENGTH);
		pch = strstr(frev, "parse");
		if (!pch)
			return 1;
		strncpy(pch, "pairs", DNAME_LENGTH);

		/* Double-check that files are mates */
		spn = strcspn(ffor, ".");
		ret = strncmp(ffor, frev, spn);
		if (ret)
		{
			logerror(lf, "%s:%d Files \'%s\' and \'%s\' do not appear to be mate-"
				     "pairs.\n", __func__, __LINE__, ffor, frev);
			return 1;
		}

		/* Read forward fastQ file into hash table */
		h = fastq_to_db(filelist[i], lf);
		if (!h)
			return 1;

		/* Print informational update to log file */
		loginfo(lf, "Attempting to pair files \'%s\' and \'%s\'.\n", ffor, frev);

		/* Align mated pairs and write to output file*/
		ret = pair_mates(filelist[i + 1], h, ffor, frev, lf);
		if (ret)
			return 1;

		/* Free allocated memory */
		free(ffor);
		free(frev);
		free_pairdb(h);
	}

	/* Print informational message to logfile */
	if (string_equal(cp->mode, "pair"))
		loginfo(lf, "Done pairing all fastQ files in \'%s\'.\n", cp->parent_indir);
	else
		loginfo(lf, "Done pairing all fastQ files in \'%s\'.\n", cp->outdir);

	/* Deallocate memory */
	for (i = 0; i < nfiles; i++)
		free(filelist[i]);
	free(filelist);

	return 0;
}
