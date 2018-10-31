/* file: get_cmdline.c
 * description: Populates command line data structure from program arguments
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <argp.h>
#include <errno.h>
#include "ddradseq.h"

extern int errno;
const char *argp_program_version = "ddradseq v1.4";
const char *argp_program_bug_address = "<dgarriga@lummei.net>";
static struct argp_option options[] =
{
  {"across",  'a', 0,      0, "Pool sequences across flow cells [default: false]"},
  {"mode",    'm', "STR",  0, "Run mode of ddradseq program [default: all]"},
  {"out",     'o', "DIR",  0, "Parent directory to write output"},
  {"csv",     'c', "FILE", 0, "CSV file with index and barcode"},
  {"dist",    'd', "INT",  0, "Edit distance for barcode matching [default: 1]"},
  {"score",   's', "INT",  0, "Alignment score to consider mates properly paired [default: 100]"},
  {"gapo",    'g', "INT",  0, "Penalty for opening a gap [default: 5]"},
  {"gape",    'e', "INT",  0, "Penalty for extending open gap [default: 1]"},
  {"pattern", 'p', "STR",  0, "Input fastQ file glob pattern to match [default: \"*.fastq.gz\""},
  {"threads", 't', "INT",  OPTION_HIDDEN, "Number of threads available for concurrency [default: 1]"},
  {0}
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	CMD *cp = state->input;

	switch (key)
	{
		case 'a':
			cp->across = true;
			break;
		case 'm':
			cp->mode = strdup(arg);
			break;
		case 'd':
			cp->dist = atoi(arg);
			break;
		case 'o':
			cp->parent_outdir = strdup(arg);
			break;
		case 's':
			cp->score = atoi(arg);
			break;
		case 'g':
			cp->gapo = atoi(arg);
			break;
		case 'e':
			cp->gape = atoi(arg);
			break;
		case 't':
			cp->nthreads = atoi(arg);
			if (cp->nthreads > 1)
				cp->mt_mode = true;
			break;
		case 'p':
			cp->glob = strdup(arg);
			break;
		case 'c':
			cp->csvfile = strdup(arg);
			break;
		case ARGP_KEY_ARG:
			if (state->arg_num >= 1)
				argp_usage(state);
			cp->parent_indir = strdup(arg);
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 1)
				argp_usage(state);
			break;
		default:
			return ARGP_ERR_UNKNOWN;
	}
	return 0;
}

static char args_doc[] = "INPUT_DIRECTORY";

static char doc[] =
"Parses fastQ files by flow cell, barcode, and/or index.\v"
"Valid run-time modes are \'parse\', \'pair\', and \'trimend\'. See https://github.com/lummeianalytics/ddradseq for documentation";

static struct argp argp = {options, parse_opt, args_doc, doc};

CMD *get_cmdline(int argc, char *argv[])
{
	char *datec = NULL;
	size_t strl = 0;
	time_t rawtime;
	struct tm *timeinfo;
	CMD *cp = NULL;

	/* Allocate memory for command line option structure */
	cp = malloc(sizeof(CMD));
	if (UNLIKELY(!cp))
	{
		perror("Memory allocation failure");
		return NULL;
	}

	/* Set argument defaults */
	cp->across = false;
	cp->mt_mode = false;
	cp->parent_indir = NULL;
	cp->parent_outdir = NULL;
	cp->outdir = NULL;
	cp->csvfile = NULL;
	cp->mode = NULL;
	cp->dist = 1;
	cp->score = 100;
	cp->gapo = 5;
	cp->gape = 1;
	cp->glob = NULL;
	cp->nthreads = 1;
	cp->lf = NULL;

	argp_parse(&argp, argc, argv, 0, 0, cp);

	/* Sanity check */
	if (!cp->mode)
		cp->mode = strdup("all");
	else if (!string_equal(cp->mode, "parse") && !string_equal(cp->mode, "pair")  &&
	    !string_equal(cp->mode, "trimend"))
	{
		fprintf(stderr, "ERROR: %s is not a valid mode.\n", cp->mode);
		return NULL;
	}
	if (!cp->csvfile && (string_equal(cp->mode, "parse") || string_equal(cp->mode, "all")))
	{
		fputs("ERROR: \'--csv\' switch is mandatory when running parse mode.\n", stderr);
		return NULL;
	}
	if (!cp->parent_outdir && (string_equal(cp->mode, "parse") || string_equal(cp->mode, "all")))
	{
		fputs("ERROR: \'--out\' switch is mandatory when running parse mode.\n", stderr);
		return NULL;
	}

	if (!cp->glob && (string_equal(cp->mode, "parse") || string_equal(cp->mode, "all")))
		cp->glob = strdup("*.fastq.gz");

	datec = malloc(DATELEN + 3u);
	if (UNLIKELY(!datec))
	{
		perror("Memory allocation failure");
		return NULL;
	}
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(datec, DATELEN + 3u, "/ddradseq-%F/", timeinfo);
	strl = strlen(cp->parent_outdir);
	if (cp->parent_outdir[strl - 1] == '/')
	{
		strl += 22u;
		cp->outdir = malloc(strl);
		if (UNLIKELY(!cp->outdir))
		{
			perror("Memory allocation failure");
			return NULL;
		}
		strcpy(cp->outdir, cp->parent_outdir);
		strcat(cp->outdir, datec + 1);
	}
	else
	{
		strl += 23u;
		cp->outdir = malloc(strl);
		if (UNLIKELY(!cp->outdir))
		{
			perror("Memory allocation failure");
			return NULL;
		}
		strcpy(cp->outdir, cp->parent_outdir);
		strcat(cp->outdir, datec);
	}
	free(datec);
	return cp;
}
