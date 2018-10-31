/* file: get_timestr.c
 * description: Fills in current time for logfile reporting
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <time.h>

#define LEN 80

int get_timestr(char *s)
{
	time_t rawtime;
	struct tm *timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(s, LEN, "%c", timeinfo);
	return 0;
}
