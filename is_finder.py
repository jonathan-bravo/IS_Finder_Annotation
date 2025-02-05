from Bio import Entrez
import pandas as pd
import math
import argparse
import gc
from tqdm import tqdm
from sys import platform
from io import BytesIO
import subprocess
from dataclasses import dataclass
from selenium import webdriver
from selenium.webdriver.firefox.options import Options as FirefoxOptions
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import xml.etree.ElementTree as ET


@dataclass
class GeneInfo:
    accession_id:str = 'NR_182530.1'
    start:int = 0
    stop:int = 0
    gene_description:str = ''
    locus:str = 'Manual'


@dataclass
class ISBrowserEntry:
    sig_seq:str = ''
    is_family:str = ''
    group:str = ''
    origin:str = ''
    bits_score:str = ''
    e_value:str = ''
    final_call:str = 'UNCLASSIFIED'


def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        metavar = '--IN_FILE',
        help = 'The excel file containing the ACLAME IDs',
        required = True
    )
    parser.add_argument(
        '-c',
        metavar = '--COL',
        type = int,
        help = 'The column containing the ACLAME IDs',
        default = 0
    )
    parser.add_argument(
        '-k',
        metavar = '--CHUNK_SIZE',
        type = int,
        help = 'The number of genes to process at once',
        default = 300
    )
    parser.add_argument(
        '-e',
        metavar = '--EMAIL',
        help = 'The email for use with Entrez',
        required = True
    )
    parser.add_argument(
        '-o',
        metavar = '--OUTFILE',
        help = 'The excel file containing the results of is_finder.py',
        default = 'IS_finder_results.xlsx'
    )
    return parser.parse_args()


def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def gen_xlsx(outfile):
    headers_df = pd.DataFrame(columns=[
        'aclame_id',
        'genebank_id',
        'accession_id',
        'start',
        'stop',
        'gene_description',
        'socus',
        'sig_seq',
        'is_family',
        'group',
        'origin',
        'bits_score',
        'e_value',
        'final_call'
    ])
    with pd.ExcelWriter(outfile, mode='w') as writer:
        headers_df.to_excel(writer, index=False)
    

def gen_chunk_list(infile, data_col, chunk_size):
    telseq_df = pd.read_excel(infile)
    aclame_id_strings = [id_string for id_string in telseq_df.iloc[:, data_col]]
    aclame_id_chunk_list = list(divide_chunks(aclame_id_strings, chunk_size))
    print(f'Chunked {len(aclame_id_strings)} IDs into {len(aclame_id_chunk_list)} chunks')
    return aclame_id_chunk_list


def fetch_gene_info(gene_ids_list):
    gene_ids = ','.join(gene_ids_list)
    handle = Entrez.efetch(db="gene", retmode="xml", id=gene_ids)
    records = handle.read()
    handle.close()
    entrez_data = ET.iterparse(BytesIO(records), events=("end",))

    del gene_ids
    del handle
    del records
    gc.collect()

    return entrez_data


def fetch_sequences(accession_ids_list):
    accession_ids = ','.join(accession_ids_list)
    handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="xml", id=accession_ids)
    records = handle.read()
    handle.close()
    entrez_data = ET.iterparse(BytesIO(records), events=("end",))

    del accession_ids
    del handle
    del records
    gc.collect()

    return entrez_data


def parse_gene_info(entrez_data):
    gene_info_list = []
    for event, gene in entrez_data:
        if event == "end" and gene.tag == 'Entrezgene':
            gene_info = GeneInfo()
            loci = gene.find('Entrezgene_locus')
            for locus in loci.iter('Gene-commentary'):
                gene_locus = gene.find('Entrezgene_gene/Gene-ref/Gene-ref_locus-tag') # ******
                accession = locus.find('Gene-commentary_accession')
                version = locus.find('Gene-commentary_version')
                descript = locus.find('Gene-commentary_products/Gene-commentary/Gene-commentary_label')
                seqs_commentary = locus.find('Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval')
                if seqs_commentary != None and accession != None and version != None:
                    gene_info.accession_id = f'{accession.text}.{version.text}'
                    gene_info.start = int(seqs_commentary.find('Seq-interval_from').text)
                    gene_info.stop = int(seqs_commentary.find('Seq-interval_to').text)
                    if descript != None:
                        gene_info.gene_description = descript.text
                    if gene_locus != None:
                        gene_info.locus = gene_locus.text
                    break
                else:
                    gene_info = GeneInfo()
            gene_info_list.append(gene_info)
        elif event == "end" and gene.tag == 'Error':
            gene_info = GeneInfo()
            gene_info.locus = 'ERROR'
            gene_info_list.append(gene_info)
    return gene_info_list


def parse_sequences(entrez_data, info_df):
    seq_list = []
    index = 0
    for event, seq in entrez_data:
        if event == "end" and seq.tag == "TSeq_sequence":
            info_entry = info_df.loc[[index]]
            seq_start = info_entry['start'].iloc[0]
            seq_stop = info_entry['stop'].iloc[0] + 1
            sequence = seq.text
            if seq_start == 0 and seq_stop == 1:
                seq_list.append('BLANK')
            else:
                seq_list.append(sequence[seq_start:seq_stop])
            index += 1

    del index
    del info_entry
    del seq_start
    del seq_stop
    del sequence
    gc.collect()

    return seq_list


def is_parser(driver, wait_time, poll_freq, seq):
    url = "https://www-is.biotoul.fr/blast.php"
    driver.get(url)
    seq_input_box = WebDriverWait(driver, wait_time, poll_frequency=poll_freq).until(EC.presence_of_element_located((By.CLASS_NAME, "seq")))
    blast_button = driver.find_element(By.CLASS_NAME, "boutonblast")
    seq_input_box.send_keys(seq)
    blast_button.click()
    results_table = WebDriverWait(driver, wait_time, poll_frequency=poll_freq).until(EC.presence_of_element_located((By.TAG_NAME, "table")))
    rows = results_table.find_elements(By.TAG_NAME, "tr")
    row1_html = rows[1].get_attribute('innerHTML').split('</td>')
    row1 = [x[4:] for x in row1_html]
    if float(row1[5]) > math.pow(10, -10):
        f_call = 'plasmid'
    else:
        f_call = 'likely IS/TE'
    is_result = ISBrowserEntry(
        row1[0].split('>')[1].split('<')[0],
        row1[1],
        row1[2],
        row1[3].split('>')[1].split('<')[0],
        row1[4].split('>')[1].split('<')[0],
        row1[5],
        f_call
    )

    del seq_input_box
    del blast_button
    del results_table
    del rows
    del row1_html
    del row1
    gc.collect()

    return is_result


def is_browser(seq_list):
    is_finder_results = []

    # need to modify the windows path, but don't have windows machine to test currently
    if platform == "win32":
        cmd = ""
        conda_path = subprocess.run(["powershell", "-Command", cmd], capture_output=True)
        binary_path = r'\bin\FirefoxApp\Contents\MacOS\firefox-bin'
    else:
        cmd = "conda info --envs | grep 'ifa' | awk '{print $(NF)}'"
        conda_path = subprocess.run([cmd], shell=True, capture_output=True).stdout.decode('ascii').strip()
        binary_path = '/bin/FirefoxApp/Contents/MacOS/firefox'

    firefox_options = FirefoxOptions()
    firefox_options.add_argument("--headless")
    firefox_options.binary_location = f'{conda_path}{binary_path}'
    driver = webdriver.Firefox(options=firefox_options)
    wait_time = 60
    poll_freq = 0.25

    for seq in tqdm(seq_list):
        if seq == 'BLANK':
            is_result = ISBrowserEntry()
        else:
            try_count = 0
            while try_count < 5:
                try:
                    is_result = is_parser(driver, wait_time, poll_freq, seq)
                    try_count = 5
                except:
                    try_count += 1
                    driver.quit()
                    driver = webdriver.Firefox(options=firefox_options)
        is_finder_results.append(is_result)
    
    driver.quit()

    del cmd
    del conda_path
    del binary_path
    del firefox_options
    del driver
    del wait_time
    del poll_freq
    gc.collect()

    return is_finder_results


def write_out(out_df, outfile):
    with pd.ExcelWriter(outfile, if_sheet_exists='overlay', mode='a') as writer:
        out_df.to_excel(
            writer,
            index = False,
            header = False,
            startrow = writer.sheets['Sheet1'].max_row
        )


def process_chunks(aclame_id_chunk_list, aclame_to_geneid_df, outfile):
    for i, chunk in enumerate(aclame_id_chunk_list):
        print(f'Working on chunk: {i+1}')

        aclame_to_geneid_df_subset = aclame_to_geneid_df[aclame_to_geneid_df['ACLAME ID'].isin(chunk)]
        gene_ids_strings = aclame_to_geneid_df_subset['Cross-references']

        info_df = pd.DataFrame()
        info_df['aclame_id'] = aclame_to_geneid_df_subset['ACLAME ID']
        info_df['genebank_id'] = aclame_to_geneid_df_subset['Cross-references']
        info_df['genebank_id'] = info_df['genebank_id'].apply(lambda x: x.split(':')[-1])
        info_df = info_df.reset_index(drop=True)
        del aclame_to_geneid_df_subset

        gene_ids_list = [gene_ids_string.split(':')[-1] for gene_ids_string in gene_ids_strings]
        entrez_data = fetch_gene_info(gene_ids_list)
        gene_info = parse_gene_info(entrez_data)
        accession_ids_list = [gene.accession_id for gene in gene_info]
        info_df = pd.concat([info_df, pd.DataFrame(gene_info)], axis=1)
        del gene_ids_strings
        del gene_ids_list
        del entrez_data
        del gene_info

        entrez_data = fetch_sequences(accession_ids_list)
        seq_list = parse_sequences(entrez_data, info_df)
        is_finder_results = is_browser(seq_list)
        info_df = pd.concat([info_df, pd.DataFrame(is_finder_results)], axis=1)
        del accession_ids_list
        del entrez_data
        del seq_list

        info_df.loc[(info_df['locus'] == 'ERROR') | (info_df['locus'] == 'Manual'), 'accession_id'] = 'NC_000000.0'
        info_df.loc[(info_df['aclame_id'].str.contains('gene:proph')) & (info_df['final_call'] == 'plasmid'), 'final_call'] = 'prophage'
        write_out(info_df, outfile)
        del info_df
        
        gc.collect()


def main():
    args = parse_inputs()
    Entrez.email = args.e
    
    gen_xlsx(args.o)
    
    aclame_id_chunk_list = gen_chunk_list(
        infile = args.i,
        data_col = args.c,
        chunk_size = args.k
    )

    aclame_to_geneid_df = pd.read_parquet(
        'aclame_genes_plasmids_0.4.parquet',
        engine='pyarrow',
        columns=['ACLAME ID', 'Cross-references']
    )

    process_chunks(
        aclame_id_chunk_list = aclame_id_chunk_list,
        aclame_to_geneid_df = aclame_to_geneid_df,
        outfile = args.o
    )


if __name__ == '__main__':
    main()