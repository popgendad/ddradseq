/* file: write_log.c
 * description: Functions for writing to ddradseq logfile and reporting errors
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: January 2017
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdarg.h>
#include "ddradseq.h"

void logerror(FILE *lf, const char *format, ...)
{
	char timestr[80];
	va_list ap;
	va_list copy;

	/* Update time string */
	get_timestr(&timestr[0]);

	fprintf(stderr, "[ddradseq: %s] ERROR -- ", timestr);
	fprintf(lf, "[ddradseq: %s] ERROR -- ", timestr);
	va_start(ap, format);
	va_copy(copy, ap);
	vfprintf(stderr, format, ap);
	vfprintf(lf, format, copy);
	va_end(ap);
	fflush(lf);
}

void loginfo(FILE *lf, const char *format, ...)
{
	char timestr[80];
	va_list ap;

	/* Update time string */
	get_timestr(&timestr[0]);

	fprintf(lf, "[ddradseq: %s] INFO -- ", timestr);
	va_start(ap, format);
	vfprintf(lf, format, ap);
	va_end(ap);
}

void logwarn(FILE *lf, const char *format, ...)
{
	char timestr[80];
	va_list ap;

	/* Update time string */
	get_timestr(&timestr[0]);

	fprintf(lf, "[ddradseq: %s] WARNING -- ", timestr);
	va_start(ap, format);
	vfprintf(lf, format, ap);
	va_end(ap);
}

void error(const char *format, ...)
{
	char timestr[80];
	va_list ap;

	/* Update time string */
	get_timestr(&timestr[0]);

	fprintf(stderr, "[ddradseq: %s] ERROR -- ", timestr);
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
}
