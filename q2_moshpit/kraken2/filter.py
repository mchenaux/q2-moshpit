import gzip
import os

import pandas as pd
from q2_demux._util import read_fastq_seqs
from q2_types.kraken2 import Kraken2OutputDirectoryFormat
from q2_types.per_sample_sequences import \
    CasavaOneEightSingleLanePerSampleDirFmt


def match_sequences_to_hits(
        seqs: CasavaOneEightSingleLanePerSampleDirFmt,
        hits: Kraken2OutputDirectoryFormat):
        
        if seqs.manifest['reverse'].isnull().all():
                sample_dic = {
                        index: {
                        'reads': [row['forward']]
                        }
                        for index, row in seqs.manifest.iterrows()
                }
        else:
                sample_dic = {
                        index: {
                        'reads': [row['forward'], row['reverse']]
                        }
                        for index, row in seqs.manifest.iterrows()
                } 
        hit_fps = []
        for dirpath,_,filenames in os.walk(str(hits)):
                for f in filenames:
                        sample_id = f.split('.')[0]
                        filepath = os.path.abspath(os.path.join(dirpath, f))
                        sample_dic[sample_id]['hits'] = filepath
                        hit_fps.append(filepath)
        return sample_dic


def filter_fastq_by_kraken2(
        sample_dic,
        classified_dir: str,
        unclassified_dir: str):
    
        classifications = {}
        for sample_id, value in sample_dic.items():
                with open(value['hits']) as hits:
                        for line in hits:
                                cols = line.split("\t")
                                classifications[cols[1]] = cols[0]
    
        for read in value['reads']:
                seqs = read_fastq_seqs(read)
                classified_fp = os.path.join(str(classified_dir), os.path.basename(read))
                unclassified_fp = os.path.join(str(unclassified_dir), os.path.basename(read))
                with gzip.open(classified_fp, "wt") as seqs_classified, gzip.open(unclassified_fp, "wt") as seqs_unclassified:
                        for seq in seqs:
                                seq_header, sequence, qual_header, qual = seq
                                classification = classifications[seq_header.strip("@").split(" ")[0]]
                                if classification == "U":
                                        seqs_unclassified.write(seq_header)
                                        seqs_unclassified.write("\n")
                                        seqs_unclassified.write(sequence)
                                        seqs_unclassified.write("\n")
                                        seqs_unclassified.write(qual_header)
                                        seqs_unclassified.write("\n")
                                        seqs_unclassified.write(qual)
                                        seqs_unclassified.write("\n")
                                elif classification == "C":
                                        seqs_classified.write(seq_header)
                                        seqs_classified.write("\n")
                                        seqs_classified.write(sequence)
                                        seqs_classified.write("\n")
                                        seqs_classified.write(qual_header)
                                        seqs_classified.write("\n")
                                        seqs_classified.write(qual)
                                        seqs_classified.write("\n")
                                else:
                                        print("raise an error")



def filter_seqs_by_kraken2_hits(
        seqs: CasavaOneEightSingleLanePerSampleDirFmt,
        hits: Kraken2OutputDirectoryFormat) -> (CasavaOneEightSingleLanePerSampleDirFmt, CasavaOneEightSingleLanePerSampleDirFmt):
    
        #need better error handling
        #refactor write code
        #test with joined file?

        sample_dictionary = match_sequences_to_hits(seqs, hits)
        classified_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
        unclassified_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
        classified_dir_fp = str(classified_seqs)
        unclassified_dir_fp = str(unclassified_seqs)
        filter_fastq_by_kraken2(
                sample_dictionary,
                classified_dir_fp,
                unclassified_dir_fp)

        return classified_seqs, unclassified_seqs
