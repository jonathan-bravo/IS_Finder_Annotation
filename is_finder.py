# conda install -c conda-forge geckodriver
# conda install -c conda-forge firefox
# binary path is then the conda env directory + /bin/FirefoxApp/Contents/MacOS/firefox-bin
# /Users/johnny/opt/anaconda3/envs/genomics
# Will probably need to set up conda anv for this step

# imports
from Bio import Entrez
import pandas as pd
import math
from tqdm import tqdm
from selenium import webdriver
#from selenium.webdriver.chrome.service import Service
#from webdriver_manager.chrome import ChromeDriverManager
#from selenium.webdriver.firefox.service import Service
from selenium.webdriver.firefox.options import Options as FirefoxOptions
from selenium.webdriver.common.by import By
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import xml.etree.ElementTree as ET


# function to chunk the data into smaller requests
def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


# read in the ids that need to be translated - turn into arg later
telseq_df = pd.read_excel("WORKER_MOBILOME_REMAINING_PLASMID_PROPHAGE_IS_TE_PARSING.xlsx")
aclame_id_strings = [id_string for id_string in telseq_df['Genes']]

# read in the translation "database" - Turn into internal data later?
#aclame_to_geneid_df = pd.read_csv("aclame_genes_plasmids_0.4.csv")
aclame_to_geneid_df = pd.read_parquet('example_pa.parquet', engine='pyarrow')

# chunk the data into smaller requests
aclame_id_chunk_list = list(divide_chunks(aclame_id_strings, 500))

print(f'Chunked {len(aclame_id_strings)} IDs into {len(aclame_id_chunk_list)} chunks')

# set email for entrez - turn into arg later
Entrez.email = "hakmonkey@gmail.com"

for i, chunk in enumerate(aclame_id_chunk_list):
    aclame_to_geneid_df_subset = aclame_to_geneid_df[aclame_to_geneid_df['ACLAME ID'].isin(chunk)]
    gene_ids_strings = aclame_to_geneid_df_subset['Cross-references']

    # start the info database
    info_df = pd.DataFrame()
    info_df['ACLAME ID'] = aclame_to_geneid_df_subset['ACLAME ID']
    info_df['Cross-references'] = aclame_to_geneid_df_subset['Cross-references']
    info_df['Cross-references'] = info_df['Cross-references'].apply(lambda x: x.split(':')[-1])
    info_df = info_df.reset_index(drop=True)

    # get list of gene ids
    gene_ids_list = [gene_ids_string.split(':')[-1] for gene_ids_string in gene_ids_strings]
    gene_ids = ','.join(gene_ids_list)
    print(f'Got {len(gene_ids_list)} genes for chunk: {i}')

    print(f'Fetching gene info for chunk: {i}')

    # fetch gene info
    handle = Entrez.efetch(db="gene", retmode="xml", id=gene_ids)
    records = handle.read()
    handle.close()

    print(f'Parsing gene info for chunk: {i}')

    accession_ids_list = []
    start_list = []
    stop_list = []
    descript_list = []
    locus_list = []

    # This is the raw parsing of the xml data
    entrez_data = ET.fromstring(records.decode('ascii'))
    for gene in tqdm(entrez_data):
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



    info_df['Accession Id'] = accession_ids_list
    info_df['Start'] = start_list
    info_df['Stop'] = stop_list
    info_df['Gene Description'] = descript_list
    info_df['Locus'] = locus_list


    accession_ids = ','.join(accession_ids_list)

    print(f'Fetching fasta sequences info for chunk: {i}')

    handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="xml", id=accession_ids)
    records = handle.read()
    handle.close()

    entrez_data = ET.fromstring(records.decode('ascii'))



    seq_list = []

    for j,seq in enumerate(entrez_data):
        info_entry = info_df.loc[[j]]
        seq_start = int(info_entry['Start'])
        seq_stop = int(info_entry['Stop']) + 1
        sequence = seq.find('TSeq_sequence').text
        if seq_start == 0 and seq_stop == 1:
            seq_list.append('BLANK')
        else:
            seq_list.append(sequence[seq_start:seq_stop])


    # IS BROWSER APP

    is_finder_results = []


    # TESTING USING FIREFOX INSTEAD BECAUSE CHROME IS UNSTABLE?
    # chrome_options = webdriver.ChromeOptions()
    # chrome_options.headless = True
    # driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=chrome_options)

    firefox_options = FirefoxOptions()
    firefox_options.add_argument("--headless")
    firefox_options.binary_location = r'/Users/johnny/opt/anaconda3/envs/genomics/bin/FirefoxApp/Contents/MacOS/firefox-bin'


    driver = webdriver.Firefox(options=firefox_options)


    print(f'Annotating with ISFinder for chunk: {i}')


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
            seq_input_box = WebDriverWait(driver, 10, poll_frequency=1).until(EC.presence_of_element_located((By.CLASS_NAME, "seq")))
            blast_button = driver.find_element(By.CLASS_NAME, "boutonblast")
            seq_input_box.send_keys(seq)
            blast_button.click()
            results_table = WebDriverWait(driver, 10, poll_frequency=1).until(EC.presence_of_element_located((By.TAG_NAME, "table")))
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

    out_df = pd.concat([info_df, pd.DataFrame(is_finder_results)], axis=1)

    out_df.loc[(out_df['Locus'] == 'ERROR') | (out_df['Locus'] == 'Manual'), 'Accession Id'] = 'NC_000000.0'

    out_df.to_excel('IS_finder_results.xlsx', index=False, sheet_name=i)