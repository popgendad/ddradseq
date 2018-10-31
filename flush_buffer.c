/* file: flush_buffer.c
 * description: Dumps a full buffer to file
 * author: Daniel Garrigan Lummei Analytics LLC
 * updated: November 2016
 * email: dgarriga@lummei.net
 * copyright: MIT license
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "khash.h"
#include "ddradseq.h"

#define MAX_ATTEMPTS 100

extern int errno;

int flush_buffer(int orient, BARCODE *bc, FILE *lf)
{
	char *filename = strdup(bc->outfile);
	char *buffer = bc->buffer;
	char *pch = NULL;
	char *errstr = NULL;
	int ret = 0;
	int fd;
	int num_attempts = 0;
	size_t len = bc->curr_bytes;
	struct flock fl = {F_WRLCK, SEEK_SET, 0, 0, 0};
	struct flock fl2;
	mode_t mode;
	gzFile gzf;

	fl.l_pid = getpid();
	memset(&fl2, 0, sizeof(struct flock));

	/* Set permissions if new output file needs to be created */
	mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH;

	/* Convert forward output file name to reverse */
	if (orient == REVERSE)
	{
		pch = strstr(filename, ".R1.fq.gz");
		strncpy(pch, ".R2", 3);
	}

	/* Get output file descriptor */
	fd = open(filename, O_WRONLY | O_CREAT | O_APPEND, mode);
	if (fd < 0)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Unable to open output file \'%s\': %s.\n", __func__,
		         __LINE__, filename, errstr);
		return 1;
	}

	/* Test if output file has lock in 30 second intervals */
	/* Will timeout after MAX_ATTEMPTS attempts to get a lock */
	fcntl(fd, F_GETLK, &fl2);
	num_attempts++;
	while (fl2.l_type != F_UNLCK)
	{
		if (num_attempts > MAX_ATTEMPTS)
		{
			logerror(lf, "%s:%d File \'%s\' is still locked after %d attempts... exiting.\n", __func__,
			         __LINE__, filename, num_attempts);
			return 1;
		}
		else
		{
			sleep(30);
			fcntl(fd, F_GETLK, &fl2);
			num_attempts++;
		}
	}

	/* Set lock on output file */
	if (fcntl(fd, F_SETLKW, &fl) == -1)
	{
		errstr = strerror(errno);
		logerror(lf, "%s:%d Failed to set lock on file \'%s\': %s.\n", __func__,
		         __LINE__, filename, errstr);
		return 1;
	}

	gzf = gzdopen(fd, "ab");
	ret = gzwrite(gzf, buffer, len);
	if (ret <= 0)
	{
		logerror(lf, "%s:%d Problem writing to output file \'%s\': %s.\n", __func__,
		         __LINE__, filename);
		return 1;
	}
	gzclose(gzf);

	/* Reset buffer */
	bc->curr_bytes = 0;
	memset(bc->buffer, 0, BUFLEN);
	bc->buffer[0] = '\0';

	/* Unlock output file */
	fl.l_type = F_UNLCK;
	fcntl(fd, F_SETLK, &fl);

	/* Close output file */
	close(fd);

	/* Free allocated memory */
	free(filename);

	return 0;
}
