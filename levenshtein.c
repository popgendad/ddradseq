/* file: levenshtein.c
 * description: Calculates the Levenshtein distance between two strings
 * author: Daniel Garrigan
 * date: October 2016
 * copyright: MIT license
 * note: Code is from https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C
*/

#include <string.h>
#include "ddradseq.h"

#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))

int levenshtein(const char *s1, const char *s2)
{
	unsigned int x = 0;
	unsigned int y = 0;
	unsigned int s1len = strlen(s1);
	unsigned int s2len = strlen(s2);
	unsigned int matrix[s2len+1][s1len+1];

	matrix[0][0] = 0;
	for (x = 1; x <= s2len; x++)
		matrix[x][0] = matrix[x-1][0] + 1;
	for (y = 1; y <= s1len; y++)
		matrix[0][y] = matrix[0][y-1] + 1;
	for (x = 1; x <= s2len; x++)
		for (y = 1; y <= s1len; y++)
			matrix[x][y] = MIN3(matrix[x-1][y] + 1, matrix[x][y-1] + 1,
						   matrix[x-1][y-1] + (s1[y-1] == s2[x-1] ? 0 : 1));
	return matrix[s2len][s1len];
}
