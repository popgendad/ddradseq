/* file: traverse_dirtree.c
 * description: Produces a sorted list of all fastQ files in the input directory tree
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <ftw.h>
#include <fnmatch.h>
#include <libgen.h>
#include "ddradseq.h"

#ifndef USE_FDS
#define USE_FDS 15
#endif

/* Globally scoped variables */
/* nftw can only work with globally scoped variables */
unsigned int n;
char *pattern;
char **f;
extern int errno;

/* Function prototypes */
static int count_fastqfiles(const char *filepath, const struct stat *info,
                            const int typeflag, struct FTW *pathinfo);
static int get_fastqfiles(const char *filepath, const struct stat *info,
                          const int typeflag, struct FTW *pathinfo);
static int count_parsefiles(const char *filepath, const struct stat *info,
	                        const int typeflag, struct FTW *pathinfo);
static int get_parsefiles(const char *filepath, const struct stat *info,
                          const int typeflag, struct FTW *pathinfo);
static int count_pairfiles(const char *filepath, const struct stat *info,
                           const int typeflag, struct FTW *pathinfo);
static int get_pairfiles(const char *filepath, const struct stat *info,
                         const int typeflag, struct FTW *pathinfo);
static int compare(const void *a, const void *b);

unsigned int traverse_dirtree(const CMD *cp, const char *caller, char ***flist)
{
	char *errstr = NULL;
	char *dirpath = NULL;
	int r = 0;
	int i = 0;
	unsigned int nfiles = 0;
	FILE *lf = cp->lf;

	if (string_equal(caller, "pair_main") || string_equal(caller, "trimend_main"))
		dirpath = cp->outdir;
	else
		dirpath = cp->parent_indir;

	/* Check validity of directory path */
	if (dirpath == NULL || *dirpath == '\0')
		return 0;

	if (cp->glob)
		pattern = strdup(cp->glob);

	/* Count number of files in directory tree */
	n = 0;
	if (string_equal(caller, "pair_main"))
		r = nftw(dirpath, count_parsefiles, USE_FDS, FTW_PHYS);
	else if (string_equal(caller, "trimend_main"))
		r = nftw(dirpath, count_pairfiles, USE_FDS, FTW_PHYS);
	else
		r = nftw(dirpath, count_fastqfiles, USE_FDS, FTW_PHYS);

	/* Check for errors */
	if (r < 0)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Directory traversal on %s failed: %s.\n", __func__, __LINE__,
		         dirpath, errstr);
		return 0;
	}
	else if (r > 0)
	{
		logerror(lf, "%s:%d Directory traversal callback failed.\n", __func__, __LINE__);
		return 0;
	}

	/* Get list of filenames */
	f = malloc(n * sizeof(char*));
	if (UNLIKELY(!f))
	{
		logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
		return 0;
	}
	n = 0;
	if (string_equal(caller, "pair_main"))
		r = nftw(dirpath, get_parsefiles, USE_FDS, FTW_PHYS);
	else if (string_equal(caller, "trimend_main"))
		r = nftw(dirpath, get_pairfiles, USE_FDS, FTW_PHYS);
	else
		r = nftw(dirpath, get_fastqfiles, USE_FDS, FTW_PHYS);

	/* Check for errors */
	if (r < 0)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Directory traversal failed: %s.\n", __func__, __LINE__,
		         errstr);
		return 0;
	}
	else if (r > 0)
	{
		logerror(lf, "%s:%d Directory traversal callback failed.\n", __func__, __LINE__);
		return 0;
	}

	/* Sort file list */
	qsort(f, n, sizeof(const char *), compare);

	/* Assign number of files */
	nfiles = n;

	/* Assign address of file list */
	*flist = malloc(n * sizeof(char*));
	for (i = 0; i < n; i++)
	{
		size_t len = strlen(f[i]) + 1;
		(*flist)[i] = malloc(len);
		memcpy((*flist)[i], f[i], len);
		free(f[i]);
	}
	free(f);
	free(pattern);

	return nfiles;
}

static int count_fastqfiles(const char *filepath, const struct stat *info,
                            const int typeflag, struct FTW *pathinfo)
{
	if (typeflag == FTW_F)
	{
		int r;
		char *fullname = strdup(filepath);
		char *fname = basename(fullname);
		r = fnmatch(pattern, fname, 0);
		if (r == 0)
			n++;
		free(fullname);
	}
	return 0;
}

static int get_fastqfiles(const char *filepath, const struct stat *info,
                          const int typeflag, struct FTW *pathinfo)
{
	if (typeflag == FTW_F)
	{
		int r;
		char *fullname = strdup(filepath);
		char *fname = basename(fullname);
		r = fnmatch(pattern, fname, 0);
		if (r == 0)
			f[n++] = strdup(filepath);
		free(fullname);
	}
	return 0;
}


static int count_parsefiles(const char *filepath, const struct stat *info,
                            const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "parse");
		if (p != NULL && q != NULL)
			n++;
	}
	return 0;
}

static int get_parsefiles(const char *filepath, const struct stat *info,
                          const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;
	size_t l = 0;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "parse");
		if (p != NULL && q != NULL)
		{
			l = strlen(filepath);
			f[n] = malloc(l + 1u);
			if (UNLIKELY(!f[n]))
			{
				error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
				return 1;
			}
			strcpy(f[n], filepath);
			n++;
		}
	}
	return 0;
}

static int count_pairfiles(const char *filepath, const struct stat *info,
                           const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "pairs");
		if (p != NULL && q != NULL)
			n++;
	}
	return 0;
}

static int get_pairfiles(const char *filepath, const struct stat *info,
                         const int typeflag, struct FTW *pathinfo)
{
	char *p = NULL;
	char *q = NULL;
	size_t l = 0;

	if (typeflag == FTW_F)
	{
		p = strstr(filepath, ".fq.gz");
		q = strstr(filepath, "pairs");
		if (p != NULL && q != NULL)
		{
			l = strlen(filepath);
			f[n] = malloc(l + 1u);
			if (UNLIKELY(!f[n]))
			{
				error("%s:%d Memory allocation failure.\n", __func__, __LINE__);
				return 1;
			}
			strcpy(f[n], filepath);
			n++;
		}
	}
	return 0;
}

static int compare(const void *a, const void *b)
{
	return strcmp(*(const char **) a, *(const char **) b);
}
