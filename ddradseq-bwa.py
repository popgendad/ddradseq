#!/usr/bin/env python2.7
#----------------------------------------------------------
# File: ddradseq-bwa.py
#
# Author: Lummei Analytics LLC
# Last updated: October 2016
#
# Description: A control script for running the read
# assembly step of the ddRadSeq pipeline
#----------------------------------------------------------
import os
import sys
import time
import glob
import re
import logging
import argparse
import textwrap
import subprocess
import multiprocessing
# Globally scoped initialization of logging class
logger = logging.getLogger("ddradseq-bwa")

def main(args):
    # Parse command line arguments
    parameters = getCommandLine(args)
    # Start pipeline logger
    initializeLog()
    # Check for write permissions on output directory
    checkPermissions(parameters.ddradseqDir)
    # Check for executable worker programs in user's path
    checkWorkers()
    # Check for files in existing output directories, if any
    checkExistingDirs(parameters.ddradseqDir)
    # Run the bwa stage of the pipeline
    runBWA(
        parameters.ddradseqDir,
        parameters.referenceFile,
        parameters.numThreads,
        parameters.threadsBWA)
    # Finish pipeline
    logger.info("ddradseq-bwa.py stage completed.")

"""
------------------------------------------------------------
execBWA()
------------------------------------------------------------
Executes the bwa application in parallel
Takes the directory of the output directory and the input
fastA reference genome file as arguments
Return value is trivial
Exits the script if bwa encounters run-time error
(e.g., non-zero return value to shell)
------------------------------------------------------------
"""
def execBWA(threadID, nthreadsBWA, ddradseqDir, referenceFile,
            fileStart, fileEnd, filenamesForwardSort, filenamesReverseSort):
    for i in range(fileStart, fileEnd):
        base = os.path.basename(filenamesForwardSort[i]);
        fileParse = re.split('[_.]', base)
        outfileBWA = "bam/smpl_" + fileParse[1] + ".bam"
        fulloutfileBWA = os.path.join(ddradseqDir, outfileBWA)
        cmdBWA = "bwa mem -t {:d} {} {} {} | samtools view -bS -T {} -o {} -".format(nthreadsBWA, referenceFile,
                                                                                     filenamesForwardSort[i], filenamesReverseSort[i], referenceFile, fulloutfileBWA)
        logger.info(
            "Thread {:d}: Running command: {}".format(threadID, cmdBWA))
        p = subprocess.Popen(
            cmdBWA,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        # Catch any run time errors
        exitCode = p.returncode
        if exitCode:
            logger.error(
                "thread {:d}: bwa command encountered run-time error".format(threadID))
            logger.error(stderr)
        else:
            logger.info(
                "thread {:d}: bwa command executed successfully".format(threadID))
    return 0

"""
------------------------------------------------------------
runBWA()
------------------------------------------------------------
Function to run bwa to map reads
Takes the directory of the output directory and the input
fastA reference genome file as arguments
Return value is trivial
Exits the script if bwa encounters run-time error
(e.g., non-zero return value to shell)
------------------------------------------------------------
"""
def runBWA(ddradseqDir, referenceFile, numThreads, nthreadsBWA):
    logger.info("Running the bwa stage of the pipeline")
    logger.info(
        "Starting to map reads to reference genome in {}.".format(referenceFile))
    # Get the stage start time
    startTime = time.time()
    # Construct sorted arrays of mate pair input file names
    fullPathForward = os.path.join(ddradseqDir, 'final/smpl_*.R1.fq.gz')
    filenamesForward = glob.glob(fullPathForward)
    if not filenamesForward:
        logger.error("bwa error: no input R1 fastQ files found")
        sys.exit('FATAL ERROR: no R1 fastQ files for bwa input found')
    fullPathReverse = os.path.join(ddradseqDir, 'final/smpl_*.R2.fq.gz')
    filenamesReverse = glob.glob(fullPathReverse)
    if not filenamesReverse:
        logger.error("bwa error: no input R2 fastQ files found")
        sys.exit('FATAL ERROR: no R2 fastQ files for bwa input found')
    filenamesForwardSort = sorted(filenamesForward)
    filenamesReverseSort = sorted(filenamesReverse)
    # Make a bwa subdirectory in outdir
    outSubDir = os.path.join(ddradseqDir, "bam")
    try:
        os.makedirs(outSubDir)
    except OSError:
        pass
    # Iterate through input fastQ files
    numFiles = len(filenamesForwardSort)
    fileStart = [0] * numThreads
    fileEnd = [0] * numThreads
    jobs = []
    filesPerThread = (numFiles + numThreads - 1) / numThreads
    for tid in range(numThreads):
        fileStart[tid] = tid * filesPerThread
        fileEnd[tid] = (tid + 1) * filesPerThread
    fileEnd[numThreads - 1] = numFiles
    for tid in range(numThreads):
        proc = multiprocessing.Process(
            target=execBWA,
            args=(tid,
                  nthreadsBWA,
                  ddradseqDir,
                  referenceFile,
                  fileStart[tid],
                  fileEnd[tid],
                  filenamesForwardSort,
                  filenamesReverseSort))
        jobs.append(proc)
    # Run jobs
    for pr in jobs:
        pr.start()
    # Wait for join to finish
    for pr in jobs:
        pr.join()
    # Print information about stage run time
    elapsedTime = time.time() - startTime
    min, sec = divmod(int(elapsedTime), 60)
    hour, min = divmod(min, 60)
    logger.info(
        "Total elapsed bwa time: {0:02d}:{1:02d}:{2:02d}".format(hour, min, sec))
    return 0
"""
------------------------------------------------------------
checkWorkers()
------------------------------------------------------------
Function to check for executable worker programs in the
user's path
Exits the script if any of the three worker programs are
not found
Return value is trivial
------------------------------------------------------------
"""
def checkWorkers():
    # Check for both bwa and samtools programs
    try:
        logger.info("Checking whether the bwa executable is in user PATH")
        pathBWA = subprocess.check_output('which bwa', shell=True)
    except:
        logger.error("bwa executable not found in user PATH")
        sys.exit('FATAL ERROR: bwa executable not found')
    logger.info("bwa executable found at: {}".format(pathBWA))
    try:
        logger.info(
            "Checking whether the samtools executable is in user PATH")
        pathSAMtools = subprocess.check_output(
            'which samtools', shell=True)
    except:
        logger.error("samtools executable not found in user PATH")
        sys.exit('FATAL ERROR: samtools executable not found')
    logger.info("samtools executable found at: {}".format(pathSAMtools))
    return 0

"""
------------------------------------------------------------
checkPermissions()
------------------------------------------------------------
Function to check for output directory write permissions
Return value is trivial
------------------------------------------------------------
"""
def checkPermissions(ddradseqDir):
    isDir = os.path.isdir(ddradseqDir)
    if isDir:
        logger.info(
            "Confirmed that specified path to output directory {} exists".format(ddradseqDir))
        logger.info(
            "Now testing for user write permissions on {}".format(ddradseqDir))
        try:
            fileName = os.path.join(ddradseqDir, "test")
            f = open(fileName, "w")
            f.close()
            os.remove(fileName)
        except Exception as e:
            logger.error("{}".format(e))
            sys.exit('FATAL ERROR: cannot write to specified output directory')
        logger.info(
            "Confirmed that the user is able to write to {}".format(ddradseqDir))
    return 0
"""
------------------------------------------------------------
checkExistingDirs()
------------------------------------------------------------
Function to check for output directory write permissions
Return value is trivial
------------------------------------------------------------
"""
def checkExistingDirs(ddradseqDir):
    logger.info(
        "Checking for existing output directories in {}".format(ddradseqDir))
    # Log the results of the subdirectory scan
    bamDir = os.path.join(ddradseqDir, 'bam')
    bam_isdir = os.path.isdir(bamDir)
    if bam_isdir:
        logger.info(
            "bwa output directory {} already exists on disk".format(bamDir))
        logger.warning("Files will be overwritten")
        # Remove all existing files in the bam subdirectory
        for root, dirs, files in os.walk(bamDir):
            for name in files:
                os.remove(os.path.join(root, name))
        # Remove existing bam subdirectory
        os.rmdir(bamDir)
    else:
        logger.info(
            "bwa output directory {} does not exist-- it will be created".format(bamDir))
    return 0
"""
------------------------------------------------------------
initializeLog()
------------------------------------------------------------
Function to setup the log file
Takes the login name of the user running the script
Return value is trivial
------------------------------------------------------------
"""
def initializeLog():
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler("ddradseq.log")
    formatter = logging.Formatter(
        '[%(name)s: %(asctime)s] ' +
        '%(levelname)s -- %(message)s', "%c")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    # Get user name
    userName = os.environ['USER']
    logger.info("ddradseq.py program started by user {}".format(userName))
    return 0
"""
------------------------------------------------------------
getCommandLine()
------------------------------------------------------------
Function to read command line arguments using the argparse
module
Returns a Namespace object with command line parameters
------------------------------------------------------------
"""
def getCommandLine(args):
    parser = argparse.ArgumentParser(
        __file__,
            description="Control script for read assembly stage of the ddRadSeq pipeline",
            usage='use "python %(prog)s --help" for more information',
            epilog='For more detailed information on how to run this script, see the overview document.',
            formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-d', '--dir',
                        required=True,
                        dest='ddradseqDir',
                        help="Directory containing input fastQ files",
                        metavar='DIR'
                        )
    parser.add_argument('-r', '--ref',
                        required=False,
                        dest='referenceFile',
                        help="Name of file with reference sequence for read mapping",
                        metavar='FILE'
                        )
    parser.add_argument('-m', '--map',
                        default=1,
                        type=int,
                        required=False,
                        dest='threadsBWA',
                        help="Number of threads available for bwa read mapping",
                        metavar='N'
                        )
    parser.add_argument('-t', '--threads',
                        default=1,
                        type=int,
                        required=False,
                        dest='numThreads',
                        help="Number of threads available for concurrency",
                        metavar='N'
                        )
    parser.add_argument('-v', '--version',
                        action='version',
                        version='ddradseq-bwa v1.2-beta'
                        )
    parameters = parser.parse_args(args)
    # Check if there are enough available CPU threads
    availableCores = multiprocessing.cpu_count()
    requestedCores = parameters.numThreads * parameters.threadsBWA
    if requestedCores > availableCores:
        parser.error(
            'FATAL ERROR: {:d} threads requested, only {:d} threads available on host'.format(
                requestedCores, availableCores))
    # Check if reference genome file is needed
    if not parameters.referenceFile:
        parser.error(
            'FATAL ERROR: reference genome fastA file needed as input')
    else:
        if not os.path.isfile(parameters.referenceFile):
            parser.error(
                'FATAL ERROR: reference fastA file {} not found'.format(
                    parameters.referenceFile))
    return parameters
if __name__ == '__main__':
    main(sys.argv[1:])
