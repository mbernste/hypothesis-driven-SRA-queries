from optparse import OptionParser
import os
from os.path import isdir, join
import json
from collections import defaultdict
import sqlite3

import load_ontology

QUERY_SRA_DB = """SELECT experiment_accession, run_accession, sample_accession, study_accession
    FROM experiment JOIN sample USING (sample_accession) JOIN run USING (experiment_accession) WHERE 
    library_strategy = 'RNA-Seq' AND scientific_name = 'Homo sapiens' AND platform = 'ILLUMINA'"""

def main():
    parser = OptionParser()
    parser.add_option("-o", "--out_dir", help="Directory in which to write output")
    (options, args) = parser.parse_args()

    metasra_f = args[0]
    sra_db_f = args[1]
    out_dir = options.out_dir

    og = load_ontology.the_ontology()

    sample_to_runs = defaultdict(lambda: [])
    sample_to_study = {}
    print('Querying the SRAdb file {}...'.format(sra_db_f))
    with sqlite3.connect(sra_db_f) as conn:
        curs = conn.cursor()
        ret = curs.execute(QUERY_SRA_DB)
        for r in ret:
            row = [x for x in r]
            exp = row[0]
            run = row[1]
            sample = row[2]
            study = row[3]
            sample_to_runs[sample].append(run)
            sample_to_study[sample] = study
    print('done.')        

    term_name_to_id = {}
    sample_to_real_val_props = defaultdict(lambda: [])
    sample_to_terms = {}
    sample_to_type = {}
    with open(metasra_f, 'r') as f:
        metasra = json.load(f)
        for sample in metasra:
            term_ids = metasra[sample]['mapped ontology terms']
            sample_type = metasra[sample]['sample type']
            real_val_props = metasra[sample]['real-value properties']
            if len(real_val_props) > 0:
                sample_to_real_val_props[sample].append(real_val_props)
            term_names = []
            for term_id in term_ids:
                term_name = og.id_to_term[term_id].name
                term_name_to_id[term_name] = term_id
                term_names.append(term_name)
            sample_to_terms[sample] = term_names
            sample_to_type[sample] = sample_type

    with open(join(out_dir, 'term_name_to_id.json'), 'w') as f:
        json.dump(term_name_to_id, f, indent=True)
    with open(join(out_dir, 'sample_to_runs.json'), 'w') as f:
        json.dump(sample_to_runs, f, indent=True)
    with open(join(out_dir, 'sample_to_study.json'), 'w') as f:
        json.dump(sample_to_study, f, indent=True)
    with open(join(out_dir, 'sample_to_real_value.json'), 'w') as f:
        json.dump(sample_to_real_val_props, f, indent=True)
    with open(join(out_dir, 'sample_to_type.json'), 'w') as f:
        json.dump(sample_to_type, f, indent=True)
    with open(join(out_dir, 'sample_to_terms.json'), 'w') as f:
        json.dump(sample_to_terms, f, indent=True)  

if __name__ == "__main__":
    main()


