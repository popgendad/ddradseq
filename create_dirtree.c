/* file: create_dirtree.c
 * description: Creates and checks output directory tree
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <errno.h>
#include <dirent.h>
#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "khash.h"
#include "ddradseq.h"

int create_dirtree(const CMD *cp, const khash_t(pool_hash) *h)
{
	char *pooldir = NULL;
	char *flowdir = NULL;
	char *parsedir = NULL;
	char *pairdir = NULL;
	char *trimdir = NULL;
	char *errstr = NULL;
	char filepath[512];
	int writable = 0;
	int status = 0;
	size_t strl = 0;
	khint_t i = 0;
	khint_t j = 0;
	khash_t(pool) *p = NULL;
	POOL *pl = NULL;
	struct dirent *next_file;
	FILE *lf = cp->lf;
	DIR *d;

	/* Check if parent output directory is writable */
	writable = access(cp->parent_outdir, W_OK);
	if (writable < 0)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Cannot write to directory \'%s\': %s.\n", __func__,
		         __LINE__, cp->parent_outdir, errstr);
		return 1;
	}
	else
	{
		/* Test if output directory already exists */
		d = opendir(cp->outdir);

		/* If the directory doesn't already exist, create it */
		if (d)
			closedir(d);
		else if (errno == ENOENT)
		{
			status = mkdir(cp->outdir, S_IRWXU | S_IRGRP | S_IXGRP |
									   S_IROTH | S_IXOTH);
			if (status < 0)
			{
				errstr = strerror(errno);
				logerror(lf, "%s:%d Failed to create output directory \'%s\': %s.\n",
				         __func__, __LINE__, cp->outdir, errstr);
				return 1;
			}
		}
		else
		{
			errstr = strerror(errno);
			logerror(lf, "%s:%d Failed to create output directory \'%s\': %s.\n",
						__func__, __LINE__, cp->outdir, errstr);
			return 1;
		}

		/* Check and create pool subdirectories */
		for (i = kh_begin(h); i != kh_end(h); i++)
		{
			if (kh_exist(h, i))
			{
				if (!cp->across)
				{
					strl = strlen(cp->outdir) + strlen(kh_key(h, i));
					flowdir = malloc(strl + 1u);
					if (UNLIKELY(!flowdir))
					{
						logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
						return 1;
					}
					strcpy(flowdir, cp->outdir);
					strcat(flowdir, kh_key(h, i));
					d = opendir(flowdir);
					if (d)
						closedir(d);
					else if (errno == ENOENT)
					{
						status = mkdir(flowdir, S_IRWXU | S_IRGRP | S_IXGRP |
												S_IROTH | S_IXOTH);
						if (status < 0)
						{
							errstr = strerror(errno);
							logerror(lf, "%s:%d Failed to create flowcell-level output "
							         "directory \'%s\': %s.\n", __func__, __LINE__,
									 flowdir, errstr);
							return 1;
						}
					}
					else
					{
						errstr = strerror(errno);
						logerror(lf, "%s:%d Failed to create flowcell-level directory \'%s\': %s.\n",
									__func__, __LINE__, flowdir, errstr);
						return 1;
					}
					free(flowdir);
				}
				p = kh_value(h, i);
				for (j = kh_begin(p); j != kh_end(p); j++)
				{
					if (kh_exist(p, j))
					{
						pl = kh_value(p, j);
						pooldir = pl->poolpath;
						d = opendir(pooldir);

						/* If the subsubdirectory doesn't already exist */
						/* create it */
						if (d)
							closedir(d);
						else if (errno == ENOENT)
						{
							status = mkdir(pooldir, S_IRWXU | S_IRGRP |
													S_IXGRP | S_IROTH | S_IXOTH);
							if (status < 0)
							{
								char *errstr = strerror(errno);
								logerror(lf, "%s:%d Failed to create pool-level "
								         "directory \'%s\': %s.\n", __func__,
										 __LINE__, pooldir, errstr);
								return 1;
							}
						}
						strl = strlen(pooldir);
						parsedir = malloc(strl + 7u);
						if (UNLIKELY(!parsedir))
						{
							logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
							return 1;
						}
						strcpy(parsedir, pooldir);
						strcat(parsedir, "/parse");
						d = opendir(parsedir);

						/* If the subsubdirectory doesn't already exist */
						/* create it */
						if (d)
						{
							/* If directory already exists-- delete all files */
							while ((next_file = readdir(d)) != NULL)
							{
								sprintf(filepath, "%s/%s", parsedir, next_file->d_name);
								remove(filepath);
							}
							closedir(d);
						}
						else if (errno == ENOENT)
						{
							status = mkdir(parsedir, S_IRWXU | S_IRGRP |
													 S_IXGRP | S_IROTH | S_IXOTH);
							if (status < 0)
							{
								char *errstr = strerror(errno);
								logerror(lf, "%s:%d Failed to create parse directory "
								         "\'%s\': %s.\n", __func__, __LINE__, parsedir,
										 errstr);
								return 1;
							}
						}
						else
						{
							errstr = strerror(errno);
							logerror(lf, "%s:%d Failed to create parse directory \'%s\': %s.\n",
									 __func__, __LINE__, parsedir, errstr);
							return 1;
						}
						free(parsedir);
						strl = strlen(pooldir);
						pairdir = malloc(strl + 7u);
						if (UNLIKELY(!pairdir))
						{
							logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
							return 1;
						}
						strcpy(pairdir, pooldir);
						strcat(pairdir, "/pairs");
						d = opendir(pairdir);

						/* If the subsubdirectory doesn't already exist */
						/* create it */
						if (d)
						{
							/* If directory already exists-- delete all files */
							while ((next_file = readdir(d)) != NULL)
							{
								sprintf(filepath, "%s/%s", pairdir, next_file->d_name);
								remove(filepath);
							}
							closedir(d);
						}
						else if (errno == ENOENT)
						{
							status = mkdir(pairdir, S_IRWXU | S_IRGRP |
													S_IXGRP | S_IROTH | S_IXOTH);
							if (status < 0)
							{
								char *errstr = strerror(errno);
								logerror(lf, "%s:%d Failed to create pairs directory "
								         "\'%s\': %s.\n", __func__, __LINE__, pairdir,
										 errstr);
								return 1;
							}
						}
						else
						{
							errstr = strerror(errno);
							logerror(lf, "%s:%d Failed to create pairs directory \'%s\': %s.\n",
										__func__, __LINE__, pairdir, errstr);
							return 1;
						}
						free(pairdir);
						strl = strlen(pooldir);
						trimdir = malloc(strl + 7u);
						if (UNLIKELY(!trimdir))
						{
							logerror(lf, "%s:%d Memory allocation failure.\n", __func__, __LINE__);
							return 1;
						}
						strcpy(trimdir, pooldir);
						strcat(trimdir, "/final");
						d = opendir(trimdir);

						/* If the subsubdirectory doesn't already exist */
						/* create it */
						if (d)
						{
							/* If directory already exists-- delete all files */
							while ((next_file = readdir(d)) != NULL)
							{
								sprintf(filepath, "%s/%s", trimdir, next_file->d_name);
								remove(filepath);
							}
							closedir(d);
						}
						else if (errno == ENOENT)
						{
							status = mkdir(trimdir, S_IRWXU | S_IRGRP |
													S_IXGRP | S_IROTH | S_IXOTH);
							if (status < 0)
							{
								char *errstr = strerror(errno);
								logerror(lf, "%s:%d Failed to create final directory "
								         "\'%s\': %s.\n", __func__, __LINE__, trimdir,
										 errstr);
								return 1;
							}
						}
						else
						{
							errstr = strerror(errno);
							logerror(lf, "%s:%d Failed to create final directory \'%s\': %s.\n",
										__func__, __LINE__, trimdir, errstr);
							return 1;
						}
						free(trimdir);
					}
				}
			}
		}
	}
	return 0;
}
