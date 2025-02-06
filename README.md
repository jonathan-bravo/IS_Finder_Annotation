# ISFinder Annotator

This script will take ACLAME ids and convert them to the associated genbank id
and then submit the fasta sequence of the genes to IS Finder and store the
results in an excel workbook.

## ACLAME to GeneBank

This is an example of the ACLAME database (v0.4) that the data was pulled from.

```
>gene:vir:3245 Length: 1335 # Range: Scmbb::Annotations::FuzzyRange=HASH(0x7fc859f4c4d8) # LocusTag: 933Wp01 # GeneName: int # MgeID: mge:66 # MgeName: 933W # Cross-refs: genbank:GeneID:1261964
GTGCTCTTGGACGCAGGAGGAACAATGGCGAATTCAGCCTATCCAGCCGGCGTTGAAAATCACGGAGGAAAACTCCGAATAACGTTTAAGTACAGGGGTAAACGAGTGCGCGAAAATCTTCGCGTGCCCGATACACCGAAAAACAGAAAGATCGCTGGTGAGTTAAGGGCTTCGGTCTGCTTTGCAATCAGAACAGGAACGTTTGATTATGCCGATCGATTCCCTGACTCACCTAACCTGAAGCTATTTGGCCTGGTAAAAAAAGATATCACCGTCGGTGAACTGGCACAGAAATGGCTTACTCTGAAAGCAATGGAAATCGGTAGTAACGCCTTAAATCGTTATCAATCAGTGATGAAAAATATGCTACCGAGGCTTGGTCCTGGCAGGCTGGCGTCATCGATTACAAAAGAAGATCTGCTGTTTATCAGGAAAGATTTACTGACCGGGGAAAAGGGAAGCAGGAAAACCAGCACGTCCCGAAAAGGAAGAACCGTACCCACAGTGAACTATTACATGACAACAACAGCCGGAATGTTCAGCTTTGCCGCCGAAAACGGGTATCTGGAGAAAAACCCGTTTAATTCAATAACACCGCTGAGGAAATCAAAACCAGTGCCGGATCCACTGACCAGAGATGAGTTTAGCCGTCTCATTGATGCCTGCCATCATCAACAGACCAAAAACCTCTGGACAGTGGCTGTTTTTACAGGGATGCGACACGGTGAAATTGCCGCACTTGCATGGGAGGATATCGACCTGAAAGCTGGCACGATAACAGTGCGACGAAATTTTACAAAAATAGGTGATTTTACGCTACCAAAGACCGACGCAGGCACTAACCGGGTTATACATCTTCTGGCACCAGCAATTGAAGCACTTAAAAACCAGGCGATGCTTACTCGTCTTAGCAGGCAGCATCAGATCACTGTTCAATTACGCGAGTACGGAAGAACAATTTTGCACGAGTGCACTTTTGTTTTCTGTCCGCAAATCGTTCGCAAGAATCACAAGGCGGGTATTAACTACGCGGTAAGCTCCATCGGAGCGACATGGGATTCAGCAATAAAAAGAGCGGGTATCCGATCCCGTAAAGCGTATCAGTCACGCCATACCTATGCGTGCTGGGCTTTATCTTCCGGAGCAAACCCGACATTTATTGCATCACAGATGGGGCACTCCAGCGCCAGCATGGTCTACAATGTTTATGGTGCATGGATGCCTGAGTGCAGCGTGACTCAAGTTGCCATGTTGAATAATGTCCTTAATGCCCGTGCCCCAGACGTGCCCCAAAGTGACCAGGAGGATGAAATAAAATTATATTTTTCAAAATGA
>gene:vir:3246 Length: 300 # Range: Scmbb::Annotations::FuzzyRange=HASH(0x7fc859fabce8) # LocusTag: 933Wp02 # GeneName: xis # MgeID: mge:66 # MgeName: 933W # Cross-refs: genbank:GeneID:1261967
ATGCAGAGGCAACTTATGCGCGAGTTAGTAAACCAACATAACCATGGCATTCAGCCAGTCATCACACCTGTTGTACAGATAAATGCGAATGAATGGGTAACACTGGAGCTTTTAATGGCTGTAACAGGCCTGAGAAAAGGAACAATATTACGCGCCAGGGACAGTGCGTGGATGAACGGCAGAGAATATAAACAAATCGCCCCCGACGGAACGCCAAAGAAAAACAGCGAATGTCTCTATCATCTTCCTACCATCAACACTTGGATCAAAAACCAACCCTTACCATCTCAGGATGTTTAA
```

The full fasta file is ~120 MB and the resulting parquet file is only 1.4 MB.
Below is the script I used to grab the ACLAME id and genebank id:

```python
from Bio import SeqIO
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

aclame_fasta = 'aclame_genes_all_0.4.fasta'
parquet_file = 'aclame_genes_plasmids_0.4.parquet'

data_rows = [(record.id, record.description.split('Cross-refs:')[-1]) for record in SeqIO.parse(aclame_fasta, 'fasta')]


df = pd.DataFrame(data_rows, columns=["ACLAME ID","Cross-references"])

parquet_schema = pa.Table.from_pandas(df).schema

parquet_writer = pq.ParquetWriter(
    parquet_file,
    parquet_schema,
    compression='snappy'
)

table = pa.Table.from_pandas(df, schema=parquet_schema)

parquet_writer.write_table(table)

parquet_writer.close()
```

## Input Flags

The script has only two required inputs `-i`, your input `XLSX` file and `-e`,
the email to use for `Entrez`.

| Flag | MetaVar | Purpose | Default | Required |
| - | - | - | - | - |
| -i | --IN_FILE | The excel file containing the ACLAME IDs | NA | True |
| -c | --COL | The column containing the ACLAME IDs | 0 | False |
| -k | --CHUNK_SIZE | The number of genes to process at once | 300 | False |
| -e | --EMAIL | The email for use with Entrez | NA | True |
| -o | --OUTFILE | The excel file containing the results of is_finder.py | IS_finder_results.xlsx | False |

## Usage

To test execution of the code create the conda environment, activate it,
then run the python script as follows:

```
conda env create -f ifa.yaml

conda activate ifa

python is_finder.py -i test.xlsx -e email@gmail.com
```

## Data Sources & Citations

An archived version of the [ACLAME v0.4](https://ngdc.cncb.ac.cn/databasecommons/database/id/1179) database was used to generate the `parquet` file.

The [test.xlsx](https://slizovskiy.com/) file came from our collaberators Slizovskiy et. al. at Purdue University.

```
@article{leplae2010aclame,
  title={ACLAME: a CLAssification of Mobile genetic Elements, update 2010},
  author={Leplae, Raphael and Lima-Mendez, Gipsi and Toussaint, Ariane},
  journal={Nucleic acids research},
  volume={38},
  number={suppl\_1},
  pages={D57--D61},
  year={2010},
  publisher={Oxford University Press}
}
```