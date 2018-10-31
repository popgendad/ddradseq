/* file: ddradseq.c
 * description: Entry point for the ddradseq program
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include "ddradseq.h"

/* Function prototypes */
extern int parse_main(const CMD*);
extern int trimend_main(const CMD*);
extern int pair_main(const CMD*);

int main(int argc, char *argv[])
{
	int ret = 0;
	CMD *cp = NULL;

	/* Parse the command line options */
	cp = get_cmdline(argc, argv);
	if (!cp)
		return 1;

	/* Initialize the log files */
	ret = log_init(cp);
	if (ret)
		return 1;

	/* Run the parse pipeline stage */
	if (string_equal(cp->mode, "parse") || string_equal(cp->mode, "all"))
	{
		ret = parse_main(cp);
		if (ret)
			return 1;
	}

	/* Run the pair pipeline stage */
	if (string_equal(cp->mode, "pair") || string_equal(cp->mode, "all"))
	{
		ret = pair_main(cp);
		if (ret)
			return 1;
	}

	/* Run the trimend pipeline stage */
	if (string_equal(cp->mode, "trimend") || string_equal(cp->mode, "all"))
	{
		ret = trimend_main(cp);
		if (ret)
			return 1;
	}

	/* Free memory for command line data structure from heap */
	destroy_cmdline(cp);

	return 0;
}
