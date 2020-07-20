from optparse import OptionParser
import os
from os.path import isdir, join
import json
from collections import defaultdict
import sqlite3

import load_ontology
from onto_lib_py3 import ontology_graph

QUERY_SRA_DB = """SELECT experiment_accession, run_accession, sample_accession, study_accession
    FROM experiment JOIN sample USING (sample_accession) LEFT JOIN run USING (experiment_accession) WHERE 
    library_strategy = 'RNA-Seq' AND scientific_name = 'Homo sapiens' AND platform = 'ILLUMINA'"""

def main():
    parser = OptionParser()
    parser.add_option("-o", "--out_dir", help="Directory in which to write output")
    (options, args) = parser.parse_args()

    metasra_f = args[0]
    sra_db_f = args[1]
    out_dir = options.out_dir

    og = load_ontology.the_ontology()
    u_og = load_ontology.unit_ontology()

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
            if run is not None:
                sample_to_runs[sample].append(run)
            sample_to_study[sample] = study
    print('done.')        

    term_name_to_id = {}
    term_id_to_name = {}
    sample_to_real_val_props = {}
    sample_to_ms_terms = {}
    sample_to_all_terms = {}
    sample_to_type = {}
    with open(metasra_f, 'r') as f:
        metasra = json.load(f)
        for sample in metasra:
            # Map sample to term ID's
            term_ids = metasra[sample]['mapped ontology terms']
            all_term_ids = set(term_ids)
            for term_id in term_ids:
                all_term_ids.update(og.recursive_superterms(term_id))
            ms_term_ids = ontology_graph.most_specific_terms(all_term_ids, og)
            #sample_to_all_terms[sample] = all_term_ids
            #sample_to_ms_terms[sample] = ms_term_ids

            sample_type = metasra[sample]['sample type']
            raw_real_val_props = metasra[sample]['real-value properties']
            real_val_props = []
            for real_val_prop in raw_real_val_props:
                prop_id = real_val_prop['property_id']
                unit_id = real_val_prop['unit_id']
                prop_name = og.id_to_term[prop_id].name
                if unit_id != 'missing' and unit_id is not None:
                    unit_name = u_og.id_to_term[unit_id].name
                else:
                    unit_name = unit_id
                new_real_val_prop = {}
                new_real_val_prop['property'] = prop_name
                new_real_val_prop['unit'] = unit_name
                new_real_val_prop['value'] = real_val_prop['value']
                real_val_props.append(new_real_val_prop)
            if len(real_val_props) > 0:
                sample_to_real_val_props[sample] = real_val_props
            ms_term_names = []
            all_term_names = []
            for term_id in ms_term_ids:
                term_name = og.id_to_term[term_id].name
                term_id_to_name[term_id] = term_name
                if term_name in term_name_to_id:
                    if not ('EFO' in term_id and 'DOID' in term_name_to_id[term_name]):
                        term_name_to_id[term_name] = term_id
                else:
                    term_name_to_id[term_name] = term_id
                ms_term_names.append(term_name)
            for term_id in all_term_ids:
                term_name = og.id_to_term[term_id].name
                term_id_to_name[term_id] = term_name
                if term_name in term_name_to_id:
                    if not ('EFO' in term_id and 'DOID' in term_name_to_id[term_name]):
                        term_name_to_id[term_name] = term_id
                else:
                    term_name_to_id[term_name] = term_id
                all_term_names.append(term_name)
            sample_to_ms_terms[sample] = ms_term_names
            sample_to_all_terms[sample] = all_term_names
            sample_to_type[sample] = sample_type

    assert set(sample_to_all_terms.keys()) <= set(sample_to_study.keys())

    with open(join(out_dir, 'term_name_to_id.json'), 'w') as f:
        json.dump(term_name_to_id, f, indent=True)
    with open(join(out_dir, 'term_id_to_name.json'), 'w') as f:
        json.dump(term_id_to_name, f, indent=True)
    with open(join(out_dir, 'sample_to_runs.json'), 'w') as f:
        json.dump(sample_to_runs, f, indent=True)
    with open(join(out_dir, 'sample_to_study.json'), 'w') as f:
        json.dump(sample_to_study, f, indent=True)
    with open(join(out_dir, 'sample_to_real_value.json'), 'w') as f:
        json.dump(sample_to_real_val_props, f, indent=True)
    with open(join(out_dir, 'sample_to_type.json'), 'w') as f:
        json.dump(sample_to_type, f, indent=True)
    with open(join(out_dir, 'sample_to_all_terms.json'), 'w') as f:
        json.dump(sample_to_all_terms, f, indent=True) 
    with open(join(out_dir, 'sample_to_ms_terms.json'), 'w') as f:
        json.dump(sample_to_ms_terms, f, indent=True)

if __name__ == "__main__":
    main()


