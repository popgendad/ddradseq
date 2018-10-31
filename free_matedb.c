/* file: free_matedb.c
 * description: Deallocates memory used by mate pair database
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdlib.h>
#include "khash.h"
#include "ddradseq.h"

int free_matedb(khash_t(mates) *m)
{
	const char *key;
	char *v = NULL;

	if (m == NULL)
		return 1;
	kh_foreach(m, key, v, free(v); free((void*)key););
	kh_destroy(mates, m);
	return 0;
}
