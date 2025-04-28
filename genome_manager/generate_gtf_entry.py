#!/usr/bin/env python3

"""
Python script to generate all standard GTF line entries (gene, transcript, exon, CDS, etc.)
from data specified in a YAML file.
"""

__author__ = 'Rob Moccia'
__version__ = '0.1'

from argparse import ArgumentParser
import yaml
from typing import Union, Optional
import os
from pydantic import BaseModel, validator

class YamlExon(BaseModel):
    """YAML schema for defining an exon"""
    start: int
    end: int
    strand: Optional[str]

class YamlTranscript(BaseModel):
    """YAML schema for defining a transcript"""
    name: str
    exons: list[YamlExon]
    cds_start: int
    cds_end: int

class YamlGeneModel(BaseModel):
    """YAML schema for defining a gene model"""
    gene_name: str
    gene_id: str
    strand: str
    gene_biotype: str
    transcripts: list[YamlTranscript]

    @validator('strand')
    def validate_strand(cls, value):
        valid_strands = ['+', '-', '.']
        if not value in valid_strands:
            raise ValueError(f'strand must be one of {valid_strands}')
        return value

class YamlGeneCollection(BaseModel):
    """YAML schema for defining a collection of gene models"""
    collection_id: str
    genes: list[YamlGeneModel]

class StrandMismatchError(Exception):
    pass

class Exon:
    """Simple class to describe an exon"""
    def __init__(self, start: int, end: int, strand: str):
        self.start = start
        self.end = end
        self.strand = strand

    @property
    def start(self):
        return self._start
    
    @start.setter
    def start(self, value: int):
        if isinstance(value, int):
            self._start = value
        else:
            raise TypeError('start must be of type int')

    @property
    def end(self):
        return self._end
    
    @end.setter
    def end(self, value: int):
        if isinstance(value, int):
            self._end = value
        else:
            raise TypeError('end must be of type int')

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, value: str):
        valid_strands = ['+', '-', '.']
        if value not in valid_strands:
            raise TypeError(f'strand must be one of {valid_strands}')
        self._strand = value

    def __str__(self):
        return (f"({self.start}, {self.end})")

    def __eq__(self, exon):
        return min([self.start, self.end]) == min([exon.start, exon.end]) and \
            max([self.start, self.end]) == max([exon.start, exon.end])

    def __ne__(self, exon):
        return not self == exon

    def __gt__(self, exon):
        forward_strand = min([self.start, self.end]) > max([exon.start, exon.end])
        if self.strand == '-':
            return not forward_strand
        else:
            return forward_strand

    def __ge__(self, exon):
        return self > exon or self == exon

    def __lt__(self, exon):
        forward_strand =  max([self.start, self.end]) < min([exon.start, exon.end])
        if self.strand == '-':
            return not forward_strand
        else:
            return forward_strand

    def __le__(self, exon):
        return self > exon or self == exon

    def contains(self, position: int):
        """Check if position is contained within the exon"""
        return self.start <= position <= self.end

    def overlap(self, position_range: Union[list, tuple, set]):
        """Report overlapping interval between two exons"""
        overlap_range = (max(self.start, position_range[0]), min(self.end, position_range[1]))
        if overlap_range[0] > overlap_range[1]:
            return None
        else:
            return overlap_range

class ExonicFeature(Exon):
    """Class to describe a feature or portion of a feature contained within an exon like a start codon or a CDS"""
    def __init__(self, start: int, end: int, exon_number: int, strand: Optional[str]=None,
            frame: Optional[Union[str, int]]=None):
        super().__init__(start=start, end=end, strand=strand)
        self.exon_number = exon_number
        self.frame = frame if frame is not None else '.'

    @property
    def exon_number(self):
        return self._exon_number

    @exon_number.setter
    def exon_number(self, value: int):
        if isinstance(value, int):
            self._exon_number = value
        else:
            raise TypeError('exon_number must be of type int')

    @property
    def frame(self):
        return self._frame

    @frame.setter
    def frame(self, value: Union[int, str]):
        valid_frames = [0, 1, 2, '.']
        if value not in valid_frames:
            raise TypeError(f'frame must be one of {valid_frames}')
        else:
            self._frame = value

class Transcript:
    """Class to define a custom transcript"""
    def __init__(self, transcript_name: str, exons: list[Exon], cds_start: int, cds_end: int,
            transcript_biotype: Optional[str]=None, transcript_source: Optional[str]=None):
        self.transcript_name = transcript_name
        self.transcript_id = transcript_name
        self.transcript_source = transcript_source if transcript_source is not None else 'custom'
        self.transcript_biotype = transcript_biotype if transcript_biotype is not None else 'processed_transcript'
        self.exons = sorted(exons)
        self.start = min([exon.start for exon in self.exons])
        self.end = max([exon.end for exon in self.exons])
        exon_strands = set([exon.strand for exon in exons])
        if len(exon_strands) != 1:
            raise StrandMismatchError('All exons in transcript must be on the same strand.')
        self.strand = exon_strands.pop()
        if self.strand == '+':
            self.start_codon_range = (cds_start, cds_start + 2)
            self.stop_codon_range = (cds_end - 2, cds_end)
            self.five_prime_utr_range = (self.start, cds_start - 1) 
            self.three_prime_utr_range = (cds_end + 1, self.end) 
            # Ensembl gtf files do not include stop codon in CDS so subtract 3 from cds_end
            self.cds_range = (cds_start, cds_end - 3)
        else:
            self.start_codon_range = (cds_end - 2, cds_end)
            self.stop_codon_range = (cds_start, cds_start + 2)
            self.five_prime_utr_range = (cds_end + 1, self.end) 
            self.three_prime_utr_range = (self.start, cds_start - 1) 
            # Ensembl gtf files do not include stop codon in CDS so add 3 from cds_start
            self.cds_range = (cds_start + 3, cds_end)
        self.start_codon = []
        self.stop_codon = []
        self.cds = []
        self.five_prime_utr = []
        self.three_prime_utr = []

        # indicator variables to track remaining start/stop codon bases being searched for
        start_codon_remain = 3
        stop_codon_remain = 3
        # variables to track frame of remaining start/stop codon and CDS sequences
        start_codon_frame = 0
        stop_codon_frame = 0
        cds_frame = 0

        for exon_num,exon in enumerate(self.exons, start=1):
            exon.exon_number = exon_num
            if start_codon_remain and (overlap := exon.overlap(self.start_codon_range)):
                self.start_codon.append(
                    ExonicFeature(start=overlap[0], end=overlap[1],
                        exon_number=exon_num, strand=exon.strand, frame=start_codon_frame))
                # determine how much of the 3 bp start codon remains to be found
                bp = overlap[1] - overlap[0] + 1
                start_codon_remain -= bp
                # if all 3 bp of start not found, update frame and next possible exon range to search
                if start_codon_remain:
                    start_codon_frame = bp % 3
                    # use exon_num to index into next exon for transcript because it is being
                    # enumerated with 1-index while the exons for the transcript are 0-index
                    # if this is a + strand transcript, the next position will be the start of the next exon
                    if self.strand == '+':
                        self.start_codon_range = (self.exons[exon_num].start, self.exons[exon_num].start + start_codon_remain - 1)
                    # if it is a - strand transcript, the next position will be the end of the next exon
                    else:
                        self.start_codon_range = (self.exons[exon_num].end - start_codon_remain + 1, self.exons[exon_num].end)

            if stop_codon_remain and (overlap := exon.overlap(self.stop_codon_range)):
                self.stop_codon.append(
                    ExonicFeature(start=overlap[0], end=overlap[1], 
                        exon_number=exon_num, strand=exon.strand, frame=stop_codon_frame))
                # determine how much of the 3 bp stop codon remains to be found
                bp = overlap[1] - overlap[0] + 1
                stop_codon_remain -= bp
                # if all 3 bp of stop not found the remainder has to be in the previous exon just processed
                # assuming no exon is < 2 bp long we can just add the remainder directly using previous exon coordinates
                if stop_codon_remain:
                    # previous exon is explicitly indexed by subtracting one twice to make it
                    # clear that it is the previous exon AND enumeration started at 1 while the indexing
                    # of exons is 0-based
                    prev_exon = self.exons[exon_num - 1 - 1]  
                    if self.strand == '+':
                        overlap = prev_exon.overlap([prev_exon.end - stop_codon_remain + 1, prev_exon.end])
                        # this also means that the last cds end will be too big by len(overlap)
                        bp = overlap[1] - overlap[0] + 1
                        self.cds[-1].end -= bp 
                    else:
                        overlap = prev_exon.overlap([prev_exon.start, prev_exon.start + stop_codon_remain - 1])
                        # this also means that the last cds start will be too small by len(overlap)
                        bp = overlap[1] - overlap[0] + 1
                        self.cds[-1].start += bp
                    stop_codon_remain -= bp
                    stop_codon_frame =  bp % 3
                    self.stop_codon.append(
                        ExonicFeature(start=overlap[0], end=overlap[1],
                            exon_number=exon_num - 1, strand=prev_exon.strand, frame=stop_codon_frame))

            if (overlap := exon.overlap(self.cds_range)):
                self.cds.append(
                    ExonicFeature(start=overlap[0], end=overlap[1],
                        exon_number=exon_num, strand=exon.strand, frame=cds_frame))
                # calculate frame of next piece of coding sequence
                cds_frame = (overlap[1] - overlap[0]) % 3
            if (overlap := exon.overlap(self.five_prime_utr_range)):
                self.five_prime_utr.append(
                    ExonicFeature(start=overlap[0], end=overlap[1],
                        exon_number=exon_num, strand=exon.strand))
            if (overlap := exon.overlap(self.three_prime_utr_range)):
                self.three_prime_utr.append(
                    ExonicFeature(start=overlap[0], end=overlap[1],
                        exon_number=exon_num, strand=exon.strand))

class Gene:
    """Class to define a custom gene added to the reference genome"""
    def __init__(self, gene_name: str, gene_id: str, strand: str, transcripts: list[Transcript],
            gene_biotype: Optional[str]=None, gene_source: Optional[str]=None):
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.strand = strand
        self.gene_biotype = gene_biotype if gene_biotype is not None else 'exogenous_protein_coding'
        self.gene_source = gene_source if gene_source is not None else 'custom'
        self.transcripts = transcripts
        # use gene_name as the artifically added chromosome name by convention
        self.chrom = self.gene_id 
        self.start = min([transcript.start for transcript in self.transcripts])
        self.end = max([transcript.end for transcript in self.transcripts])

    @property
    def gene_name(self):
        return self._gene_name

    @gene_name.setter
    def gene_name(self, value: str):
        if isinstance(value, str):
            self._gene_name = value
        else:
            raise TypeError('gene_name must of type str')

    @property
    def gene_id(self):
        return self._gene_id

    @gene_id.setter
    def gene_id(self, value: str):
        if isinstance(value, str):
            self._gene_id = value
        else:
            raise TypeError('gene_id must of type str')

    @property
    def strand(self):
        return self._strand

    @strand.setter
    def strand(self, value: str):
        valid_strands = ['+', '-', '.']
        if value not in valid_strands:
            raise TypeError(f'strand must be one of {valid_strands}')
        self._strand = value

    @property
    def transcripts(self):
        return self._transcripts
    
    @transcripts.setter
    def transcripts(self, value: list[Transcript]):
        if not isinstance(value, list):
            raise TypeError('transcripts must be a list of type Transcript')
        for val in value:
            if not isinstance(val, Transcript):
                raise TypeError('transcripts must be a list of type Transcript')
            if val.strand != self.strand:
                raise StrandMismatchError('all Transcript strands must match Gene strand')
        self._transcripts = value

    def to_gtf(self, string: bool=False) -> str:
        """Generate a full string representation of gene including all features in GTF format"""
        gene_attributes = (f'gene_id "{self.gene_id}"; gene_name "{self.gene_name}"; '
            f'gene_source "{self.gene_source}"; gene_biotype "{self.gene_biotype}";')
        gene = f'{self.chrom}\t{self.gene_source}\tgene\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{gene_attributes}'
        result = [gene]

        for transcript in self.transcripts:
            # transcript line
            transcript_attributes = (f'gene_id "{self.gene_id}"; transcript_id "{transcript.transcript_id}"; '
                f'gene_name "{self.gene_name}"; transcript_name "{transcript.transcript_name}"; '
                f'gene_source "{self.gene_source}"; gene_biotype "{self.gene_biotype}"; '
                f'transcript_source "{transcript.transcript_source}"; transcript_biotype "{transcript.transcript_biotype}";')
            result.append(f'{self.chrom}\t{transcript.transcript_source}\ttranscript\t{transcript.start}\t{transcript.end}\t.\t{transcript.strand}\t.\t{transcript_attributes}')
            
            # exon lines
            for exon in transcript.exons:
                exon_attributes = f'{transcript_attributes} exon_number "{exon.exon_number}";'
                result.append(f'{self.chrom}\t{transcript.transcript_source}\texon\t{exon.start}\t{exon.end}\t.\t{exon.strand}\t.\t{exon_attributes}')

            # start codon line
            for codon in transcript.start_codon:
                start_codon_attributes = f'{transcript_attributes} exon_number "{codon.exon_number}";'
                result.append(f'{self.chrom}\t{transcript.transcript_source}\tstart_codon\t{codon.start}\t{codon.end}\t.\t{codon.strand}\t{codon.frame}\t{start_codon_attributes}')

            # CDS lines
            for cds in transcript.cds:
                cds_attributes = f'{transcript_attributes} exon_number "{cds.exon_number}";'
                result.append(f'{self.chrom}\t{transcript.transcript_source}\tCDS\t{cds.start}\t{cds.end}\t.\t{cds.strand}\t{cds.frame}\t{cds_attributes}')

            # UTR lines
            for utr in transcript.five_prime_utr:
                result.append(f'{self.chrom}\t{transcript.transcript_source}\tfive_prime_utr\t{utr.start}\t{utr.end}\t.\t{utr.strand}\t.\t{transcript_attributes}')
            for utr in transcript.three_prime_utr:
                result.append(f'{self.chrom}\t{transcript.transcript_source}\tthree_prime_utr\t{utr.start}\t{utr.end}\t.\t{utr.strand}\t.\t{transcript_attributes}')

            # stop codon line
            for codon in transcript.stop_codon:
                start_codon_attributes = f'{transcript_attributes} exon_number "{codon.exon_number}";'
                result.append(f'{self.chrom}\t{transcript.transcript_source}\tstop_codon\t{codon.start}\t{codon.end}\t.\t{codon.strand}\t{codon.frame}\t{start_codon_attributes}')

        result = sorted([record.split('\t') for record in result], key=lambda x:(x[0], int(x[3])))
        if string:
            return '\n'.join(['\t'.join(record) for record in result])
        else:
            return result

def read_genes_from_yaml(yaml_path: Union[str, bytes, os.PathLike]) -> list[Gene]:
    """Create Gene objects containing Transcripts and Exons as defined in a yaml file"""
    genes = [] # gene accumulator

    with open(yaml_path, 'r') as f:
        yaml_args = yaml.load(f, Loader=yaml.CLoader)
    # allow alternate YAML format with top-level key 'genes' to enable the addition of other types of
    # information in the future; TODO: deprecate if not used soon
    if isinstance(yaml_args, dict):
        collection = YamlGeneCollection.parse_obj(yaml_args)
        gene_list = collection.genes

    elif isinstance(yaml_args, list):
        gene_list = [YamlGeneModel.parse_obj(arg) for arg in yaml_args]

    for gene_model in gene_list:
        transcripts = []
        for transcript_model in gene_model.transcripts:
            # Build a list of Exon objects while ignoring exon 'strand'
            # This is a complete hack to prevent users from providing a strand as part of the exon
            # entries in the yaml file. Instead, strand will be pulled from the gene model.
            exon_list = []
            for exon in transcript_model.exons:
                exon_dict = exon.dict()
                del exon_dict['strand']
                exon_list.append(Exon(strand = gene_model.strand, **exon_dict))
            transcripts.append(
                Transcript(
                    transcript_name = transcript_model.name,
                    cds_start = transcript_model.cds_start,
                    cds_end = transcript_model.cds_end,
                    exons = exon_list))
        gene = Gene(
            gene_name=gene_model.gene_name,
            gene_id=gene_model.gene_id,
            strand=gene_model.strand,
            transcripts = transcripts,
            gene_biotype=gene_model.gene_biotype)
        genes.append(gene)
    return genes

def generate_gtf(genes: list[Gene]) -> str:
    """Concatenate a list of Gene objects in gtf format"""
    result = [gene.to_gtf() for gene in genes]
    result = [gene for sub_gtf in result for gene in sub_gtf]
    result = sorted(result, key=lambda x:(x[0], int(x[3])))
    return '\n'.join(['\t'.join(record) for record in result])

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('yaml',
        help='path to YAML file containing gene specifications')
    parser.add_argument('--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()
    genes = read_genes_from_yaml(args.yaml)
    print(generate_gtf(genes))


# scrap test code:
# test_transcript = Transcript('test_transcript',
#     [Exon(10, 101), Exon(173, 502), Exon(717, 900)], 50, 481)
# test_gene = Gene('test', 'test_id', [test_transcript])

# test_transcript = Transcript('test_transcript',
#     [Exon(10, 101, strand='-'), Exon(173, 502, strand='-'), Exon(717, 900, strand='-')], 25, 493)

# test_transcript = Transcript('test_transcript',
#     [Exon(10, 101), Exon(173, 502), Exon(717, 900)], 100, 718)

# test_transcript = Transcript('test_transcript',
#     [Exon(10, 101, strand='-'), Exon(173, 502, strand='-'), Exon(717, 900, strand='-')], 100, 718)

# with open('/home/moccir/Share/code/rnaseq_pipeline/smg/bag3_vector.json', 'r') as f:
#     json_args = json.load(f)
# test_exon = Exon(10, 101)
# test_exon2 = Exon(57, 101)
# test_exon3 = Exon(63, 127)
# test_exon4 = Exon(102, 119)
# test_exon5 = Exon(12, 100)

# test_exon.overlap(test_exon)
# test_exon.overlap(test_exon2)
# test_exon.overlap(test_exon3)
# test_exon.overlap(test_exon4)
# test_exon.overlap(test_exon5)

# test_exon = Exon(10, 101, strand='what')
# test_exon.strand = 'hello'

# genes_json = json.load(args.json[0])

# genes = read_genes_from_json(genes_json)


# print(genes_yaml)
# print(type(genes_yaml))
