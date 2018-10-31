/** @file ddradseq.h
 *  @brief Header for the ddradseq program
 *  @author Daniel Garrigan Lummei Analytics LLC
 *  @date November 2016
 *
 *  Contact: dgarriga@lummei.net
 *  Copyright: MIT license
 */

#ifndef DDRADSEQ_H
#define DDRADEQ_H

#if defined __GNUC__ && !defined __GNUC_STDC_INLINE__ && !defined __GNUC_GNU_INLINE__
#define __GNUC_GNU_INLINE__ 1
#endif

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdbool.h>
#include <emmintrin.h>
#include "khash.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/** \def BUFLEN
    \brief File I/O buffer size
*/

#define BUFLEN 0x20000

/** @def MAX_LINE_LENGTH
 *  @brief Maximum line length to read from input file.
 */

#define MAX_LINE_LENGTH 400

/** @def BSIZE
 *  @brief Number of lines in individual parse buffers.
 */

#define BSIZE 4000

/** @def DNAME_LENGTH
 *  @brief Length of terminal output directory name.
 */

#define DNAME_LENGTH 5

/** @def FORWARD
 *  @brief Identifier for forward-oriented reads.
 */

#define FORWARD 1

/** @def REVERSE
 *  @brief Identifier for reverse-oriented reads.
 */

#define REVERSE 2

/** @def DATELEN
 *  @brief Length of data format YYYY-DD-MM.
 */

#define DATELEN 20

/** @def KSW_XBYTE
 *  @brief
 */

#define KSW_XBYTE  0x10000

/** @def KSW_XSTOP
 *  @brief
 */

#define KSW_XSTOP  0x20000

/** @def KSW_XSUBO
 *  @brief
 */

#define KSW_XSUBO  0x40000

/** @def KSW_XSTART
 *  @brief
 */

#define KSW_XSTART 0x80000


/******************************************************
 * Data structure defintions
 ******************************************************/

/** @var typedef struct cmdparam_t CMD
 *  @brief Data structure to hold user-provided command line parameters.
 */

typedef struct cmdparam_t
{
	bool across;          /**< Flag to pool sequences across flow cells. */
	bool mt_mode;         /**< Flag to indicate multi-threaded mode. */
	char *parent_indir;   /**< String holding the full path and name of the parent input directory. */
	char *parent_outdir;  /**< String holding the full path to the parent output directory. */
	char *outdir;         /**< String holding the full path to the output directory. */
	char *csvfile;        /**< String holding the full path of the CSV database input file. */
	char *mode;           /**< String holding the run-time mode of the program. */
	char *glob;           /**< String holding the input fastQ file glob expression. */
	int dist;             /**< The allowable edit distance for a barcode match. */
	int score;            /**< The alignment score to consider mates properly paired. */
	int gapo;             /**< The penalty for opening an alignment gap. */
	int gape;             /**< The penalty for extending an open alignment gap. */
	int nthreads;         /**< The number of threads to use for parallel computation. */
	FILE *lf;             /**< Pointer to the log file output stream. */
} CMD;

/** @var typedef struct fastq_t FASTQ
 *  @brief Data structure to hold a single fastQ entry.
 */

typedef struct fastq_t
{
	char *id;       /**< String holding the Illumina identifier line for a fastQ entry. */
	char *seq;      /**< String holding the DNA sequence line for a fastQ entry. */
	char *qual;     /**< String holding the DNA sequence quality line for a fastQ entry. */
} FASTQ;


/** @var typedef struct ksqr_t ALIGN_RESULT
 *  @brief Data structure to hold the results of a local sequence alignment.
 */

typedef struct ksqr_t
{
	int score;          /**< The score of the highest scoring local alignment. */
	int target_begin;   /**< The position on the target sequence where the highest scoring local alignment begins. */
	int target_end;     /**< The position on the target sequence where the highest scoring local alignment ends. */
	int query_begin;    /**< The position on the query sequence where the highest scoring local alignment begins. */
	int query_end;      /**< The position on the query sequence where the highest scoring local alignment ends. */
	int score2;         /**< The score of the second highest scoring local alignment. */
	int target_end2;    /**< The position on the target sequence where the second highest scoring local alignment ends. */
} ALIGN_RESULT;


/** @var typedef struct kswq_t ALIGN_QUERY
 *  @brief Data structure to hold align query parameters.
 */

typedef struct kswq_t
{
    int qlen;
    int slen;
    unsigned char shift;
    unsigned char mdiff;
    unsigned char max;
    __m128i *qp;
    __m128i *H0;
    __m128i *H1;
    __m128i *E;
    __m128i *Hmax;
} ALIGN_QUERY;


/** @var typedef struct barcode_t BARCODE
 *  @brief Barcode-level data structure.
 */

typedef struct barcode_t
{
	char *smplID;       /**< The sample identifier from the CSV database file. */
	char *outfile;      /**< The full path to the output file associated with a biological sample. */
	char *buffer;       /**< The output buffer associated with a biological sample. */
	size_t curr_bytes;  /**< The number of bytes currently in the output buffer associated with a biological sample. */
} BARCODE;

/** @def KHASH_MAP_INIT_STR(barcode, BARCODE*)
 *  @brief Defines the third-level hash
 */

KHASH_MAP_INIT_STR(barcode, BARCODE*)

/** @var typedef struct pool_t POOL
 *  @brief Pool-level data structure.
 */

typedef struct pool_t
{
	char *poolID;            /**< The sample pool identifier from the CSV database file. */
	char *poolpath;          /**< The full path to the output directory associated with a sample pool. */
	size_t barcode_length;   /**< The length of the pool identifier barcode. */
	khash_t(barcode) *b;     /**< Pointer to the hash of samples associated with this pool. */
} POOL;

/** @def KHASH_MAP_INIT_STR(pool, POOL*)
 *  @brief Defines the second-level hash
 */

KHASH_MAP_INIT_STR(pool, POOL*)

/** @def KHASH_MAP_INIT_STR(pool_hash, khash_t(pool)*)
 *  @brief Defines the top-level hash
 */

KHASH_MAP_INIT_STR(pool_hash, khash_t(pool)*)

/** @def KHASH_MAP_INIT_STR(fastq, FASTQ*)
 *  @brief Defines the hash to hold fastQ entries
 */

KHASH_MAP_INIT_STR(fastq, FASTQ*)

/** @def KHASH_MAP_INIT_STR(mates, char*)
 *  @brief Defines the hash to hold mate information
 */

KHASH_MAP_INIT_STR(mates, char*)


/******************************************************
 * Function prototypes
 ******************************************************/

/******************************************************
 * Parsing functions
 ******************************************************/

/** @fn int parse_fastq(const CMD *cp, const int orient, const char *filename, khash_t(pool_hash) *h, khash_t(mates) *m)
 *  @brief Parses a fastQ file by index sequence.
 *  @param cp Pointer to command line data structure (read-only).
 *  @param orient Orientation of reads in fastQ file (read-only).
 *  @param filename Pointer to string holding fastQ input file name (read-only).
 *  @param h Pointer to pool_hash hash table with parsing database.
 *  @param m Pointer to mate information hash table.
 *  @return Zero on success and non-zero on failure.
 */

extern int parse_fastq(const CMD *cp, const int orient, const char *filename, khash_t(pool_hash) *h, khash_t(mates) *m);


/** @fn int parse_forwardbuffer(const CMD *cp, char *buff, const size_t nl, khash_t(pool_hash) *h, khash_t(mates) *m)
 *  @brief Parses forward fastQ entries in the buffer.
 *  @param cp Pointer to command line data structure (read-only).
 *  @param buff Pointer to string holding the buffer.
 *  @param nl Number of lines in the buffer (read-only).
 *  @param h Pointer to pool_hash hash table with parsing database (read-only).
 *  @param m Pointer to mate information hash table.
 *  @return Zero on success and non-zero on failure.
 */

extern int parse_forwardbuffer(const CMD *cp, char *buff, const size_t nl, const khash_t(pool_hash) *h, khash_t(mates) *m);


/** @fn int parse_reversebuffer(const CMD *cp, char *buff, const size_t nl, khash_t(pool_hash) *h, khash_t(mates) *m)
 *  @brief Parses reverse fastQ entries in the buffer.
 *  @param cp Pointer to command line data structure (read-only).
 *  @param buff Pointer to string holding the buffer.
 *  @param nl Number of lines in the buffer (read-only).
 *  @param h Pointer to pool_hash hash table with parsing database (read-only).
 *  @param m Pointer to mate information hash table (read-only).
 *  @return Zero on success and non-zero on failure.
 */

extern int parse_reversebuffer(const CMD *cp, char *buff, const size_t nl, const khash_t(pool_hash) *h, const khash_t(mates) *m);


/******************************************************
 * Sequence pairing functions
 ******************************************************/

/** @fn int pair_mates(const char *filename, const khash_t(fastq) *h, const char *ffor, const char *frev, FILE *lf)
 *  @brief Pairs mates in two fastQ files.
 *  @param filename Pointer to string for input forward fastQ (read-only).
 *  @param h Pointer to hash table to hold forward sequences (read only).
 *  @param ffor Pointer to string with forward output file name (read only).
 *  @param frev Pointer to string with reverse output file name (read only).
 *  @param lf Pointer to log file stream.
 *  @return Zero on success and non-zero on failure.
 */

extern int pair_mates(const char *filename, const khash_t(fastq) *h, const char *ffor, const char *frev, FILE *lf);


/******************************************************
 * Trimend functions
 ******************************************************/

/** @fn int align_mates(const CMD *cp, const char *fin, const char *rin, const char *fout, const char *rout)
 *  @brief Align mates in two fastQ files and trim 3' end of reverse sequences.
 *  @param cp Pointer to command line data structure (read-only).
 *  @param fin Pointer to string with forward input file name (read-only).
 *  @param rin Pointer to string with reverse input file name (read-only).
 *  @param fout Pointer to string with forward output file name (read-only).
 *  @param rout Pointer to string with reverse output file name (read-only).
 *  @return Zero on success and non-zero on failure.
 */

extern int align_mates(const CMD *cp, const char *fin, const char *rin, const char *fout, const char *rout);


/******************************************************
 * UI functions
 ******************************************************/

/** @fn CMD *get_cmdline(int argc, char *argv[])
 *  @brief Reads command line parameters into CMD data structure.
 *  @param argc Number of command line arguments.
 *  @param argv Pointer to array of strings holding command line arguments.
 *  @return Pointer to command line data structure on success or NULL on failure.
 */

extern CMD *get_cmdline(int argc, char *argv[]);


/******************************************************
 * I/O management functions
 ******************************************************/

/** @fn khash_t(pool_hash) *read_csv(const CMD *cp)
 *  @brief Reads CSV database file into parsing hash database.
 *  @param cp Pointer to command line data structure (read-only).
 *  @return Pointer to pool_hash hash table on success or NULL on failure.
 */

extern khash_t(pool_hash) *read_csv(const CMD *cp);


/** @fn int check_csv(const CMD *cp)
 *  @brief Check the integrity of the input CSV file.
 *  @param cp Pointer to the command line parameter data structure (read-only).
 *  @return Zero on pass, non-zero on fail integrity test.
 */

extern int check_csv(const CMD *cp);


/** @fn khash_t(fastq)* fastq_to_db(const char *filename, FILE *lf)
 *  @brief Populates a fastQ database from fastQ input file.
 *  @param filename Pointer to string holding input fastQ file name (read-only).
 *  @param lf Pointer to log file stream.
 *  @return Pointer to fastQ hash table on success or NULL on failure.
 */

extern khash_t(fastq) *fastq_to_db(const char *filename, FILE *lf);


/******************************************************
 * File system functions
 ******************************************************/

/** @fn int create_dirtree(const CMD *cp, const khash_t(pool_hash) *h)
 *  @brief Creates and checks output directory tree.
 *  @param cp Pointer to command line data structure (read-only).
 *  @param h Pointer to pool_hash hash table with parsing database (read-only).
 *  @return Zero on success and non-zero on failure.
 */

extern int create_dirtree(const CMD *cp, const khash_t(pool_hash) *h);


/** @fn unsigned int traverse_dirtree(const CMD *cp, const char *caller, char ***flist)
 *  @ brief Produces a sorted list of all fastQ files in the input directory tree.
 *  @param [in] cp Pointer to command line data structure (read-only).
 *  @param [in] caller Pointer to string identifying calling function (read-only).
 *  @param [out] flist Pointer to array of input file names.
 *  @return The number of input files found in the input directory tree.
 */

extern unsigned int traverse_dirtree(const CMD *cp, const char *caller, char ***flist);


/******************************************************
 * Buffer management functions
 ******************************************************/

/** @fn char *clean_buffer(char *buff, size_t *nl)
 *  @brief Limits the buffer to hold only entire fastQ entries.
 *  @param buff Pointer to the string holding the buffer.
 *  @param nl Number of new line characters in the buffer.
 *  @return Pointer to new end of buffer.
 */

extern char *clean_buffer(char *buff, size_t *nl);


/** @fn size_t reset_buffer(char *buff, const char *r)
 *  @brief Resets the buffer for next read block.
 *  @param buff Pointer to the string holding the buffer.
 *  @param r Pointer to new end of buffer (read-only).
 *  @return The number of leftover characters in the old buffer.
 */

extern size_t reset_buffer(char *buff, const char *r);


/** @fn size_t reset_buffer(char *buff, const char *r)
 *  @brief Counts newline characters in buffer.
 *  @param buff Pointer to the string holding the buffer (read-only).
 *  @return The number of newline characters in the buffer.
 */

extern size_t count_lines(const char *buff);


/** @fn int flush_buffer(int orient, BARCODE *bc)
 *  @brief Dumps a full buffer to file.
 *  @param orient Orientation of reads in the buffer.
 *  @param bc Pointer to BARCODE data structure.
 *  @param lf Pointer to log file stream.
 *  @return Zero on success and non-zero on failure.
 */

extern int flush_buffer(int orient, BARCODE *bc, FILE *lf);


/******************************************************
 * Memory management functions
 ******************************************************/

/** @fn int destroy_cmdline(CMD *cp)
 *  @brief Destroy command line parameter data structure.
 *  @param cp Pointer to command line data structure.
 *  @return Zero on success and non-zero on failure.
 */

extern int destroy_cmdline(CMD *cp);


/** @fn int free_db(khash_t(pool_hash) *h)
 *  @brief Deallocates memory used by CSV database.
 *  @param h Pointer to pool_hash hash table with parsing database.
 *  @return Zero on success and non-zero on failure.
 */

extern int free_db(khash_t(pool_hash) *h);


/** @fn int free_pairdb (khash_t(fastq) *h)
 *  @brief Deallocates memory used by forward fastQ database.
 *  @param h Pointer to hash table to hold forward sequences.0
 *  @return Zero on success and non-zero on failure.
 */

extern int free_pairdb(khash_t(fastq) *h);


/** @fn int free_matedb(khash_t(mates) *m)
 *  @brief Deallocates memory used by mate pair database.
 *  @param m Pointer to mate information hash table.
 *  @return Zero on success and non-zero on failure.
 */

extern int free_matedb(khash_t(mates) *m);


/******************************************************
 * Alignment functions
 ******************************************************/

/** @fn ALIGN_RESULT local_align(int qlen, char *query, int tlen, char *target, const char *mat, int gapo, int gape, int xtra, FILE *lf)
 *  @brief Calculates the local sequence alignment by Smith-Waterman algorithm.
 *  @param qlen Length of query sequence.
 *  @param query Pointer string holding query sequence.
 *  @param tlen Length of the target sequence.
 *  @param target Pointer to string holding target sequence.
 *  @param mat Scoring matrix in a one-dimension array (read-only).
 *  @param gapo Gap penalty.
 *  @param gape Gap extension penalty.
 *  @param xtra Status variable.
 *  @param lf Pointer to log file stream.
 *  @return ALIGN_RESULT data structure on success
 */

extern ALIGN_RESULT local_align(int qlen, char *query, int tlen, char *target, const char *mat, int gapo, int gape, int xtra, FILE *lf);


/** @fn char *revcom(const char *s, FILE *lf)
 *  @brief Reverse complement a DNA string with full IUPAC alphabet.
 *  @param s Pointer to string to be reverse-complemented (read-only).
 *  @param lf Pointer to log file stream.
 *  @return Pointer to newly allocated reverse-complemented string on success or NULL on failure.
 */

extern char *revcom(const char *s, FILE *lf);


/******************************************************
 * Edit distance functions
 ******************************************************/

/** @fn int levenshtein(const char *s1, const char *s2)
 *  @brief Calculates the Levenshtein distance between two strings.
 *  @param s1 Pointer to the first string (read-only).
 *  @param s2 Pointer to the second string (read-only).
 *  @return Edit distance between the two strings.
 */

extern int levenshtein(const char *s1, const char *s2);


/******************************************************
 * Log file functions
 ******************************************************/

/** @fn int log_init(CMD *cp)
 *  @brief Initialize the ddradseq log file.
 *  @param cp Pointer to command line data structure (read-only).
 *  @return Zero on success and non-zero on failure.
 */

extern int log_init(CMD *cp);


/** @fn void loginfo(FILE *lf, const char *format)
 *  @brief Write informational message to log file.
 */

extern void loginfo(FILE *lf, const char *format, ...);


/** @fn void logwarn(FILE *lf, const char *format)
 *  @brief Write warning message to log file.
 */

extern void logwarn(FILE *lf, const char *format, ...);


/** @fn void logerror(FILE *lf, const char *format)
 *  @brief Report error message to both the logfile and standard error.
 */

extern void logerror(FILE *lf, const char *format, ...);


/** @fn int get_timestr(char *s);
 *  @ brief Fill in current time for log file reporting.
 *  @param s Pointer to time string.
 *  @return Zero on success and non-zero on failure.
 */

extern int get_timestr(char *s);

/******************************************************
 * Error reporting functions
 ******************************************************/

/** @fn void error(const char *format)
 *  @brief Report error to stderr
 */

extern void error(const char *format, ...);


/******************************************************
 * Inline utility functions
 ******************************************************/

/** @fn inline int string_equal(const char *a, const char *b)
 *  @brief Compare two strings for equality
 *  @param a Pointer to first string (read-only).
 *  @param b Pointer to second string (read-only).
 *  @return Non-zero if strings are equal and zero if the strings are not equal.
 */

static inline int string_equal(const char *a, const char *b)
{
	return (strcmp(a, b) == 0);
}

#endif
