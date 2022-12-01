# conda install -c conda-forge geckodriver
# conda install -c conda-forge firefox
# binary path is then the conda env directory + /bin/FirefoxApp/Contents/MacOS/firefox-bin
# /Users/johnny/opt/anaconda3/envs/genomics
# Will probably need to set up conda anv for this step

# imports
from Bio import Entrez
import pandas as pd
import math
import argparse
from tqdm import tqdm
from sys import platform
import subprocess
from selenium import webdriver
from selenium.webdriver.firefox.options import Options as FirefoxOptions
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import xml.etree.ElementTree as ET


def parse_inputs():
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        '-i',
        metavar = '--IN_FILE',
        type = str,
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
        type = str,
        help = 'The email for use with Entrez',
        required = True
    )
    parser.add_argument(
        '-o',
        metavar = '--OUTFILE',
        type = str,
        help = 'The column containing the ACLAME IDs',
        default = 'IS_finder_results.xlsx'
    )
    return parser.parse_args()


def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def gen_xlsx(outfile):
    headers_df = pd.DataFrame(columns=[
        'ACLAME ID',
        'Cross-references',
        'Accession Id',
        'Start',
        'Stop',
        'Gene Description',
        'Locus',
        'Sequences producing significant alignments',
        'IS Family',
        'Group',
        'Origin',
        'Score (bits)',
        'E. value',
        'Final Call'
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

    entrez_data = ET.fromstring(records.decode('ascii'))

    return entrez_data


def fetch_sequences(accession_ids_list):
    accession_ids = ','.join(accession_ids_list)

    handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="xml", id=accession_ids)
    records = handle.read()
    handle.close()

    entrez_data = ET.fromstring(records.decode('ascii'))

    return entrez_data


def parse_gene_info(entrez_data):
    accession_ids_list = []
    start_list = []
    stop_list = []
    descript_list = []
    locus_list = []

    for gene in entrez_data:
        if gene.tag != 'Error':
            accession_val = 'NR_182530'
            version_val = 1
            start_val = 0
            stop_val = 0
            descript_val = ''
            gene_locus_val = 'Manual'
            loci = gene.find('Entrezgene_locus')
            for locus in loci.iter('Gene-commentary'):
                gene_locus = gene.find('Entrezgene_gene/Gene-ref/Gene-ref_locus-tag') # ******
                accession = locus.find('Gene-commentary_accession')
                version = locus.find('Gene-commentary_version')
                descript = locus.find('Gene-commentary_products/Gene-commentary/Gene-commentary_label')
                seqs_commentary = locus.find('Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval')
                if seqs_commentary != None and accession != None and version != None:
                    accession_val = accession.text
                    version_val = version.text
                    start_val = seqs_commentary.find('Seq-interval_from').text
                    stop_val = seqs_commentary.find('Seq-interval_to').text
                    if descript != None:
                        descript_val = descript.text
                    if gene_locus != None:
                        gene_locus_val = gene_locus.text
                    break
                else:
                    accession_val = 'NR_182530'
                    version_val = 1
                    start_val = 0
                    stop_val = 0
                    descript_val = ''
                    gene_locus_val = 'Manual'
        else:
            accession_val = 'NR_182530' # NC_006500.2
            version_val = 1
            start_val = 0
            stop_val = 0
            descript_val = ''
            gene_locus_val = 'ERROR'
        accession_ids_list.append(f'{accession_val}.{version_val}')
        start_list.append(start_val)
        stop_list.append(stop_val)
        descript_list.append(descript_val)
        locus_list.append(gene_locus_val)

    return (
        accession_ids_list,
        start_list,
        stop_list,
        descript_list,
        locus_list
    )


def parse_sequences(entrez_data, info_df):
    seq_list = []

    for i, seq in enumerate(entrez_data):
        info_entry = info_df.loc[[i]]
        seq_start = int(info_entry['Start'])
        seq_stop = int(info_entry['Stop']) + 1
        sequence = seq.find('TSeq_sequence').text
        if seq_start == 0 and seq_stop == 1:
            seq_list.append('BLANK')
        else:
            seq_list.append(sequence[seq_start:seq_stop])

    return seq_list
    

def is_browser(seq_list):
    is_finder_results = []

    # need to modify the windows path, but don't have windows machine to test currently
    if platform == "win32":
        cmd = ""
        conda_path = subprocess.run(["powershell", "-Command", cmd], capture_output=True)
        binary_path = r'\bin\FirefoxApp\Contents\MacOS\firefox-bin'
    else:
        cmd = "conda info --envs | grep 'genomics' | awk '{print $(NF)}'"
        conda_path = subprocess.run([cmd], shell=True, capture_output=True).stdout.decode('ascii').strip()
        binary_path = '/bin/FirefoxApp/Contents/MacOS/firefox-bin'

    firefox_options = FirefoxOptions()
    firefox_options.add_argument("--headless")
    firefox_options.binary_location = f'{conda_path}{binary_path}'
    driver = webdriver.Firefox(options=firefox_options)
    wait_time = 30

    for seq in tqdm(seq_list):
        if seq == 'BLANK':
            is_result = {
                'Sequences producing significant alignments': '',
                'IS Family': '',
                'Group': '',
                'Origin': '',
                'Score (bits)': '',
                'E. value': '',
                'Final Call': 'UNCLASSIFIED'
            }
        else:    
            url = "https://www-is.biotoul.fr/blast.php"
            driver.get(url)
            seq_input_box = WebDriverWait(driver, wait_time, poll_frequency=1).until(EC.presence_of_element_located((By.CLASS_NAME, "seq")))
            blast_button = driver.find_element(By.CLASS_NAME, "boutonblast")
            seq_input_box.send_keys(seq)
            blast_button.click()
            results_table = WebDriverWait(driver, wait_time, poll_frequency=1).until(EC.presence_of_element_located((By.TAG_NAME, "table")))
            rows = results_table.find_elements(By.TAG_NAME, "tr")
            row1_html = rows[1].get_attribute('innerHTML').split('</td>')
            row1 = [x[4:] for x in row1_html]
            if float(row1[5]) > math.pow(10, -10):
                f_call = 'plasmid'
            else:
                f_call = 'likely IS/TE'
            is_result = {
                'Sequences producing significant alignments': row1[0].split('>')[1].split('<')[0],
                'IS Family': row1[1],
                'Group': row1[2],
                'Origin': row1[3].split('>')[1].split('<')[0],
                'Score (bits)': row1[4].split('>')[1].split('<')[0],
                'E. value': row1[5],
                'Final Call': f_call
            }
        is_finder_results.append(is_result)

    driver.quit()
    return is_finder_results


def process_chunks(aclame_id_chunk_list, aclame_to_geneid_df, outfile):
    for i, chunk in enumerate(aclame_id_chunk_list):
        print(f'Working on chunk: {i}')

        aclame_to_geneid_df_subset = aclame_to_geneid_df[aclame_to_geneid_df['ACLAME ID'].isin(chunk)]
        gene_ids_strings = aclame_to_geneid_df_subset['Cross-references']

        info_df = pd.DataFrame()
        info_df['ACLAME ID'] = aclame_to_geneid_df_subset['ACLAME ID']
        info_df['Cross-references'] = aclame_to_geneid_df_subset['Cross-references']
        info_df['Cross-references'] = info_df['Cross-references'].apply(lambda x: x.split(':')[-1])
        info_df = info_df.reset_index(drop=True)

        gene_ids_list = [gene_ids_string.split(':')[-1] for gene_ids_string in gene_ids_strings]

        entrez_data = fetch_gene_info(gene_ids_list)
        gene_info = parse_gene_info(entrez_data)
        accession_ids_list = gene_info[0]
        start_list = gene_info[1]
        stop_list = gene_info[2]
        descript_list = gene_info[3]
        locus_list = gene_info[4]


        info_df['Accession Id'] = accession_ids_list
        info_df['Start'] = start_list
        info_df['Stop'] = stop_list
        info_df['Gene Description'] = descript_list
        info_df['Locus'] = locus_list

        entrez_data = fetch_sequences(accession_ids_list)
        seq_list = parse_sequences(entrez_data, info_df)

        is_finder_results = is_browser(seq_list)

        out_df = pd.concat([info_df, pd.DataFrame(is_finder_results)], axis=1)

        out_df.loc[(out_df['Locus'] == 'ERROR') | (out_df['Locus'] == 'Manual'), 'Accession Id'] = 'NC_000000.0'

        write_out(out_df, outfile)


def write_out(out_df, outfile):
    with pd.ExcelWriter(outfile, if_sheet_exists='overlay', mode='a') as writer:
        out_df.to_excel(
            writer,
            index = False,
            header = False,
            startrow = writer.sheets['Sheet1'].max_row
        )

def main():
    args = parse_inputs()
    infile = args.i
    outfile = args.o
    data_col = args.c
    chunk_size = args.k
    Entrez.email = args.e #"hakmonkey@gmail.com"
    gen_xlsx(outfile)
    aclame_id_chunk_list = gen_chunk_list(infile, data_col, chunk_size)
    aclame_to_geneid_df = pd.read_parquet('aclame_genes_plasmids_0.4.parquet', engine='pyarrow')
    process_chunks(aclame_id_chunk_list, aclame_to_geneid_df, outfile)


if __name__ == '__main__':
    main()