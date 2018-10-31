/* file: destroy_cmdline.c
 * description: Destroy command line parameter data structure
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include "ddradseq.h"

int destroy_cmdline(CMD *cp)
{
	fclose(cp->lf);
	free(cp->parent_indir);
	free(cp->parent_outdir);
	free(cp->outdir);
	free(cp->mode);
	free(cp->glob);
	free(cp->csvfile);
	free(cp);
	return 0;
}
