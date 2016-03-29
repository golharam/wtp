#!/usr/bin/python
import os
import sys
import seqUtils

REF_FILE = open( sys.argv[1], 'r', 10000000 )
MAPPING_FILE = open(sys.argv[2], 'r', 10000000 )
ORIGINAL_READ_FILE = open(sys.argv[3], 'r', 10000000)
INDEX_FILE = open(sys.argv[4], 'w' )
readLength = int(sys.argv[5])
readOffset = int(sys.argv[6]) - 1
matchScore = int(sys.argv[7])
mismatchScore = int(sys.argv[8])
matchOnIupac = int(sys.argv[9])
validAdjacentMismatchesCountAsOne = int(sys.argv[10])
maskOfReads = []
for char in sys.argv[11] :  
    maskOfReads.append(int(char)) 
while(len(maskOfReads) < readLength) :
    maskOfReads.append(maskOfReads[-1])

secondHalf = False
if readOffset > 0:
    secondHalf = True


# read reference file
seqNum = 0
refSeqs = [""]
pieces = []
ii = 0
while 1:
    line = REF_FILE.readline()
    if not line: break
    line = line.strip()
    if line.startswith('>'):
        if seqNum > 0:
            bigString = ''.join( pieces )
            refSeqs.append( bigString )
        pieces = []
        seqNum = seqNum + 1
    else:
        pieces.append( line )
    ii = ii + 1
    if ii % 1000000 == 0:
        print >> sys.stderr, "Read ", ii, " lines from reference..."
REF_FILE.close()

bigString = ''.join( pieces )
refSeqs.append( bigString )

print >> sys.stderr, "Read ", ii, " lines from reference..."

indexCounter = 0
strLen = 0
while 1:
    line = MAPPING_FILE.readline()
    if not line: break
    while not line.startswith( ">" ):
        line = MAPPING_FILE.readline()
    line = line.strip()

    originalReadline = ORIGINAL_READ_FILE.readline()
    if not originalReadline: break
    while not originalReadline.startswith( ">" ):
        originalReadline = ORIGINAL_READ_FILE.readline()
    originalReadline = originalReadline.strip()

    if line.startswith( ">" ):
        coords = line.split(",") 

        readName = coords[0]
        readName = readName[1:]
        del coords[0]

        shortRead = MAPPING_FILE.readline()
        shortRead = shortRead.strip()

# read original sequence
        origLineSplit = originalReadline.split(",")
        secondReadName = origLineSplit[0]
        secondReadName = secondReadName[1:]
        #while readName != secondReadName:
	while not readName.startswith(secondReadName):
            read = ORIGINAL_READ_FILE.readline()
            if not read:
                print >> sys.stderr, "Error did not find matching read for ", readName
                exit(-1)
            originalReadline = ORIGINAL_READ_FILE.readline()
            originalReadline = originalReadline.strip()
            origLineSplit = originalReadline.split(",")
            secondReadName = origLineSplit[0]
            secondReadName = secondReadName[1:]

        read = ORIGINAL_READ_FILE.readline()
        read = read.strip()
        
	if indexCounter % 1000000 == 0:
	    INDEX_FILE.write( repr(indexCounter) + '\t' + repr(strLen) + '\n' )
        
        indexCounter += 1

        str = ">" + readName
        sys.stdout.write( str )
        strLen += len( str )

        for coord in coords:
            [foundCoord, numMismatch] = coord.split(".")
            [seqNum, pos ] = foundCoord.split("_")

            if int(pos) < 0:
                strand = '-'
            else:
                strand = '+'

            extend = readLength
            start = int(pos) - readOffset
            end = start + extend
            if strand == '-':
                end = abs(int(pos)) + 1 + readOffset
                start = end - extend 

            reference = refSeqs[ int(seqNum) ]
            
            #Extending a read off the ends of the sequence...skip.
            if start < 0 or end >= len(reference):
               continue
            
            seq = seqUtils.getSeq( reference, start, end, strand )
            seq = seq.upper()
#            print readName, seqName, pos, strand, numMismatch, start, end

            primerBase = read[0:1]

            seq = primerBase + seq


# start at position 23 as we are allowing 3 mismatches
#            ii = 23
            ii = 1
#            print colourSeq
#            print read
            (score, alignmentEnds, numMatches, numMismatches, matchString, numAdjacentValids, adjacentValidsStr) = seqUtils.alignRead( read, seq, secondHalf, matchScore, mismatchScore, matchOnIupac, validAdjacentMismatchesCountAsOne, maskOfReads )
#            alignString = " "
#            if not secondHalf:
#                ii = 1
#                while ii <= alignmentEnds:
#                    alignString += '*'
#                    ii += 1
#            else:
#                ii = 1
#                while ii < alignmentEnds:
#                    alignString += ' '
#                    ii += 1
#                while ii < len(colourSeq):
#                    alignString += '*'
#                    ii += 1

#            print score, alignmentEnds, numMatches, numMismatches
#            print " " + matchString
#            print alignString

            alignmentPosition = 1
            alignmentSize = alignmentEnds
            if secondHalf:
                alignmentPosition = alignmentEnds
# seq is +1 longer than actual seq due to primer base, hence don't need to add one in following
                alignmentSize = len(seq) - alignmentEnds

            refPos = repr(start) 
            if secondHalf:
                  refPos = repr(start + readLength - alignmentSize) 
            if strand == '-': 
                        refPos = '-' + repr(end - alignmentSize) 
                        if secondHalf:
                                     refPos = '-' + repr(start)

            str =  "," + seqNum + "." + refPos + "." + repr(alignmentPosition) + "." + repr(alignmentSize) + "." + repr(score) + "." + repr(numMismatches) + "#" + adjacentValidsStr + "##"
            sys.stdout.write( str )
            strLen += len( str )
####
#            colourSeq = seqUtils.baseToColours( seq )
####

        print
        print read
####
#        print primerBase + colourSeq
#        print seq
####
        strLen += 1 + len(read) + 1

MAPPING_FILE.close()
ORIGINAL_READ_FILE.close()
INDEX_FILE.close()
