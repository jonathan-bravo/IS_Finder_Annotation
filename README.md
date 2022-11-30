# Annotator

This script will take ACLAME ids and convert them to the associated genbank id
and then submit the fasta sequence of the genes to IS Finder and store the
results in an excel workbook.

## ACLAME to GeneBank

This is an example of the csv file that was originally used to convert the
ACLAME ids to GeneBank ids.


| ACLAME ID | Length | Genomic range | Gene name/locus tag | MGE ID | MGE Name | Cross-references |
| - | - | - | - | - | - | - |
| gene:plasmid:10466 | 2895 | complement(1458..4352) | tnpA - pAC5_p1 | mge:185 | pAC5 | genbank:GeneID:874852 |
| gene:plasmid:10467 | 549 | 4547..5095 | tnpR - pAC5_p2 | mge:185 | pAC5 | genbank:GeneID:874851 |
| gene:plasmid:10468 | 1011 | 158..1168 | pAP12875p01 | mge:186 | pAP12875 | genbank:GeneID:3276701 |
| gene:plasmid:10469 | 1026 | complement(173..1198) | pAP12875p02 | mge:186 | pAP12875 | genbank:GeneID:3276702 |
| gene:plasmid:10470 | 663 | complement(197..859) | mobE - pTC-F14_p01 | mge:187 | pTC-F14 | genbank:GeneID:1076408 |
| gene:plasmid:10471 | 681 | complement(852..1532) | mobD - pTC-F14_p02 | mge:187 | pTC-F14 | genbank:GeneID:1076409 |
| gene:plasmid:10472 | 396 | complement(1529..1924) | mobC - pTC-F14_p03 | mge:187 | pTC-F14 | genbank:GeneID:1076410 |
| gene:plasmid:10473 | 312 | 2209..2520 | mobB - pTC-F14_p04 | mge:187 | pTC-F14 | genbank:GeneID:1076411 |
| gene:plasmid:10474 | 2670 | 2510..5179 | mobA - pTC-F14_p05 | mge:187 | pTC-F14 | genbank:GeneID:1076412 |
| gene:plasmid:10475 | 1059 | 4121..5179 | repB - pTC-F14_p06 | mge:187 | pTC-F14 | genbank:GeneID:1076404 |

This table contains 96641 entries at came to 9 MB so it was converted to
parquet format using the following short script reducing its size considerably.

```python
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

csv_file = 'aclame_genes_plasmids_0.4.csv'
parquet_file = 'aclame_genes_plasmids_0.4.parquet'
chunksize = 100_000

csv_stream = pd.read_csv(
    csv_file,
    sep=',',
    chunksize=chunksize,
    low_memory=False
)

for i, chunk in enumerate(csv_stream):
    print("Chunk", i)
    if i == 0:
        # Guess the schema of the CSV file from the first chunk
        parquet_schema = pa.Table.from_pandas(df=chunk).schema
        # Open a Parquet file for writing
        parquet_writer = pq.ParquetWriter(
            parquet_file,
            parquet_schema,
            compression='snappy'
        )
    # Write CSV chunk to the parquet file
    table = pa.Table.from_pandas(chunk, schema=parquet_schema)
    parquet_writer.write_table(table)

parquet_writer.close()
```
