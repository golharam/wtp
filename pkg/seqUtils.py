#!/usr/bin/python
from string import maketrans

complementTrans = maketrans( 'ACTGUMRWSYKVHDBXNactgumrwsykvhdbxn.', 'TGACAKYWSRMBDHVXXtgacakywsrmbdhvxx.' )

baseToColourCodes = { 'AA' : '0',
                      'CC' : '0',
                      'GG' : '0',
                      'TT' : '0',
                      'AC' : '1',
                      'CA' : '1',
                      'GT' : '1',
                      'TG' : '1',
                      'AG' : '2',
                      'GA' : '2',
                      'CT' : '2',
                      'TC' : '2',
                      'AT' : '3',
                      'TA' : '3',
                      'CG' : '3',
                      'GC' : '3' }

def baseToColours( seqIn ):

    seqIn = seqIn.upper()
    prevBase = seqIn[0]
    seqIn = seqIn[1:]
    colourSeq = ''
    for base in seqIn:
        transition = prevBase + base
        colourSeq += baseToColourCodes[transition]
        prevBase = base

    return colourSeq

def revComp ( seqIn ):
    seqOut = seqIn[::-1]
    seqOut = seqOut.translate( complementTrans )
    return seqOut


def getSeq( reference, start, end, strand ):
    start = int(start)
    end = int(end)
    savedLine = reference[start:end]

    if strand == '-':
        savedLine = savedLine[::-1]
        savedLine = savedLine.translate( complementTrans )

    return savedLine

def getBasePossibilities( iupacCode ):
    basePossibilities = ""
    if iupacCode == 'A':
        basePossibilities = "A"
    elif iupacCode == 'C':
        basePossibilities = "C"
    elif iupacCode == 'G':
        basePossibilities = "G"
    elif iupacCode == 'T':
        basePossibilities = "T"
    elif iupacCode == 'R':
        basePossibilities = "AG"
    elif iupacCode == 'Y':
        basePossibilities = "CT"
    elif iupacCode == 'S':
        basePossibilities = "GC"
    elif iupacCode == 'W':
        basePossibilities = "AT"
    elif iupacCode == 'K':
        basePossibilities = "GT"
    elif iupacCode == 'M':
        basePossibilities = "AC"
    elif iupacCode == 'B':
        basePossibilities = "CGT"
    elif iupacCode == 'D':
        basePossibilities = "AGT"
    elif iupacCode == 'H':
        basePossibilities = "ACT"
    elif iupacCode == 'V':
        basePossibilities = "ACG"
    elif iupacCode == 'N':
        basePossibilities = "ACGT"

    return basePossibilities

def doTheyMatch(transitionBases, readColour, matchOnIupac):
    if baseToColourCodes.has_key( transitionBases ):
        if baseToColours( transitionBases ) == readColour:
            return True
        else:
            return False
    elif not matchOnIupac:
        return False

    firstBasePossibilities = getBasePossibilities( transitionBases[0] )
    secondBasePossibilities = getBasePossibilities( transitionBases[1] )

    for base1 in firstBasePossibilities:
        for base2 in secondBasePossibilities:
            seq = base1 + base2
            if baseToColours(seq) == readColour:
                return True
    return False

def areAdjacentColoursValid( colour1, colour2, seq1, seq2, seq3 ):
    for base in "AGCT":
        tr1 = seq1 + base
        tr2 = base + seq3
#        print tr1, tr2, baseToColourCodes[tr1], baseToColourCodes[tr2]
        if baseToColourCodes.has_key( tr1 ) and baseToColourCodes.has_key( tr2 ):
            if baseToColourCodes[tr1] == colour1 and baseToColourCodes[tr2] == colour2:
#                print"Adj valid"
                return True

#    print"Adj invalid"
    return False

def alignRead( read, seq, secondHalf, matchScore, mismatchScore, matchOnIupac, validAdjacentMismatchesCountAsOne, maskOfReads ):
    numMismatches = 0
    numMismatchesStore = 0
    numMatches = 0
    numMatchesStore = 0
    alignmentEnds = 0
    score = 0
    scoreStore = 0
    matchString = ""
    currNumMismatches = 0
    numAdjacentValids = 0
    adjacentValidsStr = ""
    ii = 1
    if secondHalf:
        ii = len(seq) - 1
    while (not secondHalf and ii < len(seq)) or (secondHalf and ii > 0):
        if (maskOfReads[ii-1] == 0) :
            currNumMismatches = 0
            matchString += 'M'
        else:
            if doTheyMatch( seq[(ii-1):(ii+1)], read[ii], matchOnIupac ):
                numMatches += 1
                score += matchScore
                matchString += '+'
                currNumMismatches = 0
            else:
                currNumMismatches += 1
                numMismatches += 1
                score += mismatchScore
                matchString += '-'
                if currNumMismatches > 1:
                    if (not secondHalf and areAdjacentColoursValid( read[ii-1], read[ii], seq[ii-2], seq[ii-1], seq[ii] )):
                        numAdjacentValids += 1
                        if numAdjacentValids > 1:
                            adjacentValidsStr += "." + repr(ii-1)
                        else:
                            adjacentValidsStr = repr(ii-1)
                        if validAdjacentMismatchesCountAsOne:
                            numMismatches -= 1
                            score -= mismatchScore
                            score += matchScore
                            currNumMismatches = 0

                    if (secondHalf and areAdjacentColoursValid( read[ii], read[ii+1], seq[ii-1], seq[ii], seq[ii+1] )):
                        numAdjacentValids += 1
                        if numAdjacentValids > 1:
                            adjacentValidsStr += "." + repr(ii)
                        else:
                            adjacentValidsStr = repr(ii)
                        if validAdjacentMismatchesCountAsOne : 
                            numMismatches -= 1
                            score -= mismatchScore
                            score += matchScore
                            currNumMismatches = 0

        if score > scoreStore:
            alignmentEnds = ii
            scoreStore = score
            numMismatchesStore = numMismatches
            numMatchesStore = numMatches

        if not secondHalf:
            ii += 1
        else:
            ii -= 1

    if secondHalf:
        matchString = matchString[::-1]
    return (scoreStore, alignmentEnds, numMatchesStore, numMismatchesStore, matchString, numAdjacentValids, adjacentValidsStr)
