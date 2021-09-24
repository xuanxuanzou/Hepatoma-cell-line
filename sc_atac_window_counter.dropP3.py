import sys
import pysam

if len(sys.argv) != 6:
	sys.exit('Usage: python sc_atac_window_counter.py [Input Bam file] [Input Index table or "NoTags"] [Window BED] [Output file] [Include sites with no reads? (True/False)]')

inbam = sys.argv[1]
indextable = sys.argv[2]
type1bed = sys.argv[3]
outfile = sys.argv[4]
includezeroes = sys.argv[5]

if indextable != 'NoTags':
	descer = open(indextable,'r')
	cells = [x.strip().split()[0] for x in descer.readlines() if '@' not in x]
	descer.close()
else:
	cells = ['Undefined']

cellsdic = {}
for x,cell in enumerate(cells):
	cellsdic[cell] = x

def lister(bedfile):
	currfile = open(bedfile,'r')
	currrecout = [line.strip().split()[0:3] for line in currfile]
	currfile.close()
	return currrecout

print("Building window map...")
rec1list = lister(type1bed)
bamfile = pysam.Samfile(inbam,'rb')

def counter(bedtuple,outsfile,first=False):
	templen = len(cells)
	if first:
		print("peak\t" + "\t".join(cells), file=outsfile)
	if includezeroes:
		for rec in bedtuple:
			recname = rec[0] + "_" + rec[1] + "_" + rec[2]
			currcounts = [0]*templen
			reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
			for read in reads:
				#readname = read.qname.split(':')[0]
				readname = read.get_tag("DB")
				try:
					currcounts[cellsdic[readname]] += 1
				except KeyError:
					pass
			print(recname + "\t" + "\t".join([str(x) for x in currcounts]), file=outsfile)

	else:
		for rec in bedtuple:
			recname = rec[0] + "_" + rec[1] + "_" + rec[2]
			currcounts = [0]*templen
			reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
			for read in reads:
				#readname = read.qname.split(':')[0]
				readname = read.get_tag("DB")
				try:
					currcounts[cellsdic[readname]] += 1
				except KeyError:
					pass
			if sum(currcounts) > 0:
				print(recname + "\t" + "\t".join([str(x) for x in currcounts]), file=outsfile)


outmat = open(outfile,'w')
print("Counting DHS reads...")
#counter(rec1list,outmat,writebinary,first=True)
counter(rec1list,outmat,first=True)
outmat.close()
