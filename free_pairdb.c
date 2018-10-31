/* file: free_pairdb.c
 * description: Deallocates memory used by pairs database
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include "khash.h"
#include "ddradseq.h"

int free_pairdb(khash_t(fastq) *h)
{
	const char *key;
	FASTQ *e = NULL;

	if (h == NULL)
		return 1;
	kh_foreach(h, key, e, free(e->id); free(e->seq); free(e->qual); free(e); free((void*)key););
	kh_destroy(fastq, h);
	return 0;
}
