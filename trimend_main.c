/* file: trimend_main.c
 * description: Entry point for the trimend function
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ddradseq.h"

int trimend_main(const CMD *cp)
{
	char *pch = NULL;
	char **filelist = NULL;
	int ret = 0;
	unsigned int i = 0;
	unsigned int nfiles = 0;
	FILE *lf = cp->lf;

	/* Print informational message to log file */
	loginfo(lf, "Beginning to trim 3\' end of reverse sequences in \'%s\'.\n", cp->outdir);

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
		pch = strstr(ffor, "pairs");
		if (!pch)
			return 1;
		strncpy(pch, "final", DNAME_LENGTH);
		pch = strstr(frev, "pairs");
		if (!pch)
			return 1;
		strncpy(pch, "final", DNAME_LENGTH);

		/* Double-check that files are mates */
		spn = strcspn(ffor, ".");
		ret = strncmp(ffor, frev, spn);
		if (ret)
		{
			logerror(lf, "%s:%d Files \'%s\' and \'%s\' do not appear to be mate-"
				     "pairs.\n", __func__, __LINE__, ffor, frev);
			return 1;
		}

		/* Print informational update to log file */
		loginfo(lf, "Attempting to align sequences in \'%s\' and \'%s\'.\n", ffor, frev);

		/* Align mated pairs and write to output file*/
		ret = align_mates(cp, filelist[i], filelist[i+1], ffor, frev);
		if (ret)
			return 1;

		/* Free allocated memory */
		free(ffor);
		free(frev);
	}

	/* Print informational message to log file */
	if (string_equal(cp->mode, "trimend"))
		loginfo(lf, "Done trimming 3\' end of reverse sequences in \'%s\'.\n", cp->parent_indir);
	else
		loginfo(lf, "Done trimming 3\' end of reverse sequences in \'%s\'.\n", cp->outdir);

	/* Deallocate memory */
	for (i = 0; i < nfiles; i++)
		free(filelist[i]);
	free(filelist);

	return 0;
}
