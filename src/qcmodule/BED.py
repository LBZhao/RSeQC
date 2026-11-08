#!/usr/bin/env python
'''
This python3 module was converted from python2.7 code using 2to3
'''
#import built-in modules
import os,sys
#import gc


#import third-party modules
#import psyco_full
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *
from itertools import *
from importlib.metadata import version


#changes to the paths

#changing history to this module
#05/26/2011: suppport multiple spliced mapped reads

__version__ = version("RSeQC")
__author__ = "Liguo Wang"


def _open_bed_records(path):
    """Return the data lines of a BED file as a rewindable text stream.

    Blank lines and ``#``/``track``/``browser`` header lines are dropped, so
    every ParseBED method sees BED records only (UCSC downloads start with a
    ``track`` line, which used to crash several methods).
    """
    import io
    with open(path, 'r') as handle:
        kept = [
            line for line in handle
            if line.strip() and not line.startswith(('#', 'track', 'browser'))
        ]
    return io.StringIO(''.join(kept))


class ParseBED:
    '''manipulate BED (http://genome.ucsc.edu/FAQ/FAQformat.html) format file.'''
    
    def __init__(self,bedFile):
        '''This is constructor of ParseBED'''
        self.transtab = str.maketrans('ACGTNX', 'TGCANX')
        self.f = _open_bed_records(bedFile)
        self.fileName=os.path.basename(bedFile)
        self.ABS_fileName=bedFile
        
    def getUTR(self,utr=35):
        '''Extract UTR regions from input bed file (must be 12-column). output is 6-column bed format.
        When utr=35 [default], extract both 5' and 3' UTR. When utr=3, only extract 3' UTR. When utr=5,
        only extract 5' UTR'''
        
        ret_lst=[]
        for line in self.f:
            if line.startswith('#'):continue
            if line.startswith('track'):continue
            if line.startswith('browser'):continue
            fields=line.rstrip('\r\n').split()
            chrom=fields[0]
            geneName=fields[3]
            strand=fields[5]
            txStart=int(fields[1])
            txEnd=int(fields[2])
            cdsStart=int(fields[6])
            cdsEnd=int(fields[7])       
            if cdsStart == cdsEnd:      # non-coding transcript: no UTRs (see getNoncodingExon)
                continue
            exon_start=list(map(int,fields[11].rstrip(',').split(',')))
            exon_start=list(map((lambda x: x + txStart),exon_start))
                
            exon_end=list(map(int,fields[10].rstrip(',').split(',')))
            exon_end=list(map((lambda x,y:x+y),exon_start,exon_end))
            
            if strand == '+':
                if (utr==35 or utr==5):
                    for st,end in zip(exon_start,exon_end):
                        if st < cdsStart:
                            utr_st = st
                            utr_end = min(end,cdsStart)
                            ret_lst.append([chrom,utr_st,utr_end,geneName,'0',strand])                  
                if (utr==35 or utr==3):
                    for st,end in zip(exon_start,exon_end):
                        if end > cdsEnd:
                            utr_st = max(st, cdsEnd)
                            utr_end = end
                            ret_lst.append([chrom,utr_st,utr_end,geneName,'0',strand])
            if strand == '-':
                if (utr==35 or utr==3):
                    for st,end in zip(exon_start,exon_end):
                        if st < cdsStart:
                            utr_st = st
                            utr_end = min(end,cdsStart)
                            ret_lst.append([chrom,utr_st,utr_end,geneName,'0',strand])                  
                if (utr==35 or utr==5):
                    for st,end in zip(exon_start,exon_end):
                        if end > cdsEnd:
                            utr_st = max(st, cdsEnd)
                            utr_end = end
                            ret_lst.append([chrom,utr_st,utr_end,geneName,'0',strand])
        self.f.seek(0)
        return ret_lst
            
    def getExon(self):
        '''Extract exon regions from input bed file (must be 12-column). output is 6-column Tab 
        separated bed file, each row represents one exon'''
        
        ret_lst=[]
        for f in self.f:
            f = f.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            name = f[4]
            strand = f[5]
            cdsStart = int(f[6])
            cdsEnd = int(f[7])
            blockCount = int(f[9])
            blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
            blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
            for base,offset in zip( blockStarts, blockSizes ):
                ret_lst.append((chrom, base, base+offset))
        self.f.seek(0)
        return ret_lst

    def getNoncodingExon(self):
        """Exons of non-coding transcripts (thickStart == thickEnd) as
        ``[chrom, start, end]`` lists.

        Before RSeQC 5.1.0+lbzhao these exons were reported as 5' or 3' UTR,
        depending on where the BED file put thickStart/thickEnd.
        """
        ret_lst = []
        for line in self.f:
            f = line.split()
            if int(f[6]) != int(f[7]):
                continue
            chrom_start = int(f[1])
            sizes = [int(i) for i in f[10].strip(',').split(',')]
            starts = [chrom_start + int(i) for i in f[11].strip(',').split(',')]
            for base, size in zip(starts, sizes):
                ret_lst.append([f[0], base, base + size])
        self.f.seek(0)
        return ret_lst

    def getCDSExon(self):
        
        '''Extract CDS exon regions from input bed file (must be 12-column).'''     
        ret_lst=[]
        for f in self.f:
            f = f.strip().split()
            chrom = f[0]
            chrom_start = int(f[1])
            name = f[4]
            strand = f[5]
            cdsStart = int(f[6])
            cdsEnd = int(f[7])
            if cdsStart == cdsEnd:      # non-coding transcript: no CDS
                continue
            blockCount = int(f[9])
            blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
            blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
            # grab cdsStart - cdsEnd
            cds_exons = []
            cds_seq = ''
            genome_seq_index = []
            for base,offset in zip( blockStarts, blockSizes ):
                if (base + offset) < cdsStart: continue
                if base > cdsEnd: continue
                exon_start = max( base, cdsStart )
                exon_end = min( base+offset, cdsEnd ) 
                #cds_exons.append( (exon_start, exon_end) )
                ret_lst.append([chrom,exon_start,exon_end]) 
        self.f.seek(0)
        return ret_lst

    def getIntron(self):
        '''Extract Intron regions from input bed file (must be 12-column).  output is 6-column Tab 
        separated bed file, each row represents one intron'''

        ret_lst=[]
        intron_label = 'unknown' # where intron is located (UTR3, UTR5 or CDS)
        for line in self.f:
            try:
                if line.startswith('#'):continue
                if line.startswith('track'):continue
                if line.startswith('browser'):continue   
                # Parse fields from gene tabls
                fields = line.split()
                chrom     = fields[0]
                tx_start  = int( fields[1] )
                tx_end    = int( fields[2] )
                geneName      = fields[3]
                strand    = fields[5].replace(" ","_")
                cds_start = int( fields[6] )
                cds_end   = int( fields[7] )
                if int(fields[9]) == 1:
                    continue
            
                exon_starts = list(map( int, fields[11].rstrip( ',\n' ).split( ',' ) ))
                exon_starts = list(map((lambda x: x + tx_start ), exon_starts))
                exon_ends = list(map( int, fields[10].rstrip( ',\n' ).split( ',' ) ))
                exon_ends = list(map((lambda x, y: x + y ), exon_starts, exon_ends));   
                intron_start = exon_ends[:-1]
                intron_end=exon_starts[1:]
                
                if(strand == '-'):
                    intronNum=len(intron_start)
                    for st,end in zip(intron_start,intron_end):
                        if st > cds_end:
                            intron_label = 'UTR5'
                        elif end < cds_start:
                            intron_label = 'UTR3'
                        else:
                            intron_label = 'CDS'
                        #FO.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t0\t" + strand + '\n')
                        #intronNum -= 1
                        ret_lst.append([chrom, st, end, geneName, intron_label])
                else:
                    intronNum=1
                    for st,end in zip(intron_start,intron_end):
                        if st > cds_end:
                            intron_label = 'UTR3'
                        elif end < cds_start:
                            intron_label = 'UTR5'
                        else:
                            intron_label = 'CDS'

                        #FO.write(chrom + "\t" + str(st) + "\t" + str(end) + "\t" + geneName + "_intron_" + str(intronNum) + "\t0\t" + strand + '\n')
                        #intronNum += 1
                        ret_lst.append([chrom, st, end, geneName, intron_label])
            except (IndexError, ValueError):
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=' ', file=sys.stderr)
                continue
        self.f.seek(0)
        return ret_lst

    def getIntergenic(self,direction='up', size=1000):
        '''get intergenic regions. direction=up or down or both.'''
        
        ret_lst=[]
        for line in self.f:
            if line.startswith(('#','track','browser')):
                continue
            fields = line.split()
            chrom     = fields[0]
            tx_start  = int( fields[1] )
            tx_end    = int( fields[2] )
            strand    = fields[5]
            if(direction=="up" or direction=="both"):
                if strand=='-':
                    region_st=tx_end
                    region_end=tx_end +size
                else:
                    region_st = max(tx_start-size,0)
                    region_end=tx_start
                ret_lst.append([chrom,region_st,region_end])
            if (direction=="down" or direction=="both"):
                if strand == '-':
                    region_st = max(0,tx_start-size)
                    region_end = tx_start
                else:
                    region_st = tx_end
                    region_end = tx_end+size
                ret_lst.append([chrom,region_st,region_end])
        self.f.seek(0)
        return ret_lst

        
def unionBed3(lst):
    '''Take the union of 3 column bed files. return a new list'''
    bitsets = binned_bitsets_from_list(lst)
    ret_lst=[]
    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while 1:
            start = bits.next_set( end )
            if start == bits.size: break
            end = bits.next_clear( start )
            ret_lst.append([chrom, start, end])
    bitsets=dict()
    return ret_lst

def subtractBed3(lst1,lst2):
    '''subtrack lst2 from lst1'''
    bitsets1 = binned_bitsets_from_list(lst1)
    bitsets2 = binned_bitsets_from_list(lst2)
    
    ret_lst=[]
    for chrom in bitsets1:  
        if chrom not in bitsets1:
            continue
        bits1 = bitsets1[chrom]
        if chrom in bitsets2:
            bits2 = bitsets2[chrom]
            bits2.invert()
            bits1.iand( bits2 )
        end=0
        while 1:
            start = bits1.next_set( end )
            if start == bits1.size: break
            end = bits1.next_clear( start )
            ret_lst.append([chrom,start,end])
    bitsets1 = dict()
    bitsets2 = dict()
    return ret_lst

def tillingBed(chrName,chrSize,stepSize=10000):
    '''tilling whome genome into small sizes'''
    #tilling genome
    for start in range(0,chrSize,stepSize):
        end = start + stepSize
        if end < chrSize:
            yield (chrName,start,end)
        else:
            yield (chrName,start,chrSize)
        
