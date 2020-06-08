from Bio import SeqIO
import Bio.Blast.NCBIWWW as blast
import xml.etree.ElementTree as ET
import sys
import pandas as pd
from Bio import  Entrez
df=pd.DataFrame()
fasta=sys.argv[1:]
for fasta_file in fasta:
    with open(fasta_file) as fasta_open:
        for seq_record in SeqIO.parse(fasta_open, 'fasta'):
            if len(seq_record)>200:
                counter=0
                result_handle = blast.qblast("blastn", "nt", str(seq_record.seq))
                v=result_handle.read()
                tree = ET.fromstring(v)
                for a in tree.iter('Hit_id'):
                    while counter <= 5:
                        cc=a.text.split('|')[3]
                        handle = Entrez.efetch(db="nucleotide", id=cc, rettype="gb")
                        for gb_record in SeqIO.parse(handle, "genbank"):
                            tax=gb_record.annotations['taxonomy']
                            organism=gb_record.annotations['organism']
                            handle = Entrez.esearch(db="Taxonomy", term=organism, )
                            record = Entrez.read(handle)
                            taxid=record["IdList"][0]
                            df=df.append({'sample':fasta_file,'seq_id':str(seq_record.id),'hit_id':cc,'organism':organism,'tax_id':taxid}, ignore_index=True)
df.to_csv('test.tsv', index=False, sep='\t')
