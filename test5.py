from BCBio import GFF
from Bio import SeqIO

input_file = "ccruddi.gbff"
output_file = "ccruddi.gff"

with open(input_file) as in_handle:
    with open(output_file, "w") as out_handle:
        GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
