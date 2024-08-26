# Code from: https://gencore.bio.nyu.edu/how-to-find-out-what-barcodes-are-in-your-undetermined-reads/

## Usage: python3 count-barcode-freq.py <fastq_file.gz>
## Example: python3 count-barcode-freq.py sample.fastq.gz
from operator import itemgetter
import sys, gzip

barcodes = {}
with gzip.open(sys.argv[1]) as fastq:
        for line in fastq:
                if not line.startswith(b'@'): continue
                bc = line.decode("utf-8").split(':')[-1].strip()
                if bc not in barcodes:
                        barcodes[bc] = 1
                else:
                        barcodes[bc]+=1

total = sum(barcodes.values())
for k, v in sorted(barcodes.items(), key=itemgetter(1)):
        print(k, v, round(v/total*100, 2))
