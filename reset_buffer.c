/* file: reset_buffer.c
 * description: Resets the buffer for next read block
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: October 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <string.h>

size_t reset_buffer(char *buff, const char *r)
{
	size_t br = strlen(r);
	memmove(buff, r, br);
	return br;
}
