from Bio import SeqIO
import sys

product = str(sys.argv[3])

with open(sys.argv[2], 'w') as nfh:
        for rec in SeqIO.parse(sys.argv[1], "genbank"):
                if rec.features:
                        for feature in rec.features:
                                if feature.type == "CDS":
					if feature.qualifiers['product'][0] == product:
	                                        nfh.write(">%s|%s from %s\n%s\n" % (
						feature.qualifiers['gene'][0],
        	                                feature.qualifiers['product'][0],
                	                        rec.name,
						feature.location.extract(rec).seq))
