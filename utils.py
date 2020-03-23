#   Functions called from the Jupyter notebooks
#   for implementing the Case-Control Finder and
#   the Series Finder.

from collections import defaultdict
import json
import pandas as pd
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns

TREATMENT_TERMS = set([
    'treatment',
    'transfection'
])

def term_to_samples(sample_to_terms, term):
    samples_with_term = []
    samples_without_term = []
    for samples, terms in sample_to_terms.items():
        if term in terms:
            samples_with_term.append(samples)
        else:
            samples_without_term.append(samples)
    assert len(frozenset(samples_with_term) & frozenset(samples_without_term)) == 0
    return samples_with_term, samples_without_term


def _is_poor_quality(terms, term_name_to_id):
    found_tissue = False
    #found_cell_line = False
    found_cell_type = False
    for term in terms:
        term_id = term_name_to_id[term]
        if 'UBERON' in term_id and term != 'male organism' \
            and term != 'female organism' and term != 'adult organism' \
            and term != 'organ':
            found_tissue = True
        #elif 'CVCL' in term_id:
        #    found_cell_line = True
        elif 'CL' in term_id and 'CVCL' not in term_id \
            and term != 'cultured cell' \
            and term != 'cell' and term != 'eukaryotic cell' \
            and term != 'animal cell' and term != 'native cell':
            found_cell_type = True
    #return not (found_tissue or found_cell_line or found_cell_type)
    return not (found_tissue or found_cell_type)  

def _is_diseased(terms, term_name_to_id):
    """
    Determine whether sample may be a diseased
    sample
    """
    for term in terms:
        if term == 'disease' or 'DOID' in term_name_to_id[term]:
            return True
    return False


def _is_treated(terms, term_name_to_id):
    """
    Determine whether sample may be a treated
    sample
    """
    for term in terms:
        if term in TREATMENT_TERMS:
            return True
    return False


def _is_cell_line(terms, term_name_to_id):
    """
    Determine whether a given ontology term is describing
    a cell line.
    """
    for term in terms:
        if 'CVCL' in term_name_to_id[term]:
            return True
    return False


def series(term, target_property, sample_to_real_val, sample_to_terms, sample_to_type, 
        sample_to_study, term_name_to_id, filter_disease=True, 
        filter_poor=True, filter_cell_line=True, filter_differentiated=True,
        target_unit=None, value_limit=None, skip_missing_unit=False
    ):
    age_to_samples = defaultdict(lambda: set())
    poor_samples = set()
    cell_line_samples = set()
    differentiated_samples = set()
    disease_samples = set()
    for sample, real_val_infos in sample_to_real_val.items():
        if sample not in sample_to_terms:
            continue
        for real_val_info in real_val_infos:
            property_ = real_val_info['property']
            unit = real_val_info['unit']
            value = int(real_val_info['value'])
            if target_unit and unit != target_unit:
                continue
            if property_ == target_property:
                if value_limit and value > value_limit:
                    continue
                terms = sample_to_terms[sample]
                #if len(blacklist_terms & set(sample_to_terms[sample])) > 0:
                #    continue
                if term in terms:
                    age_to_samples[value].add(sample)
                if _is_poor_quality(terms, term_name_to_id):
                    poor_samples.add(sample)
                if _is_diseased(terms, term_name_to_id):
                    disease_samples.add(sample)
                if sample_to_type[sample] == 'cell line':
                    cell_line_samples.add(sample)
                if filter_differentiated \
                    and (sample_to_type[sample] == 'in vitro differentiated cells'
                    or sample_to_type[sample] == 'induced pluripotent stem cell line'):
                    differentiated_samples.add(sample)

    for age in age_to_samples:
        if filter_poor:
            age_to_samples[age] -= poor_samples
        if filter_cell_line:
            age_to_samples[age] -= cell_line_samples
        if filter_disease:
            age_to_samples[age] -= disease_samples
        if filter_differentiated:
            age_to_samples[age] -= differentiated_samples
    

    da = []
    for age in sorted(age_to_samples.keys()):
        for sample in age_to_samples[age]:
            da.append((
                sample,
                sample_to_study[sample],
                age,
                sample in poor_samples,
                sample in cell_line_samples,
                sample in differentiated_samples,
                sample in disease_samples
            ))
    df = pd.DataFrame(data=da, columns=[
        'sample', 'study',
        'age', 'missing_metadata',
        'cell_line', 'differentiated',
        'diseased'
    ])
    return age_to_samples, df


def _create_key_terms(terms, term_name_to_id):
    #print('original terms: ', terms)
    term_set = set([
        term for term in terms
        if ('UBERON' in term_name_to_id[term]
        or 'CL' in term_name_to_id[term])
        and 'CVCL' not in term_name_to_id[term]
    ])
    #print('Now its: ', term_set)
    term_set -= set([
        'male organism',
        'female organism',
        'adult organism',
        'organ'    
    ])
    term_set -= set([
        'cultured cell',
        'cell',
        'eukaryotic cell',
        'animal cell',
        'native cell'
    ])
    assert len(term_set) > 0
    return '\n'.join(sorted(term_set))
    
    
def match_case_to_controls(term, control_samples, case_samples, sample_to_terms, 
    sample_to_study, term_name_to_id, sample_to_type, 
    filter_poor=True, filter_treated=True, filter_disease=True, 
    filter_cell_line=True, filter_differentiated=True, 
    by_run=False, sample_to_runs=None):
    filtered = set()
    control_samples = set(control_samples)
    case_samples = set(case_samples)

    #for sample in control_samples:
    #    if len(blacklist_terms & set(sample_to_terms[sample])) == 0:
    #        filtered.add(sample)
    #control_samples = filtered

    control_term_to_samples = defaultdict(lambda: set())
    case_term_to_samples = defaultdict(lambda: set())

    # Identify poor quality, in vitro differentiated,
    # and cell line samples
    poor_samples = set()
    cell_line_samples = set()
    differentiated_samples = set()
    disease_samples = set()
    treated_samples = set()
    for sample in set(control_samples) | set(case_samples):
        terms = sample_to_terms[sample]
        if _is_poor_quality(terms, term_name_to_id):
            poor_samples.add(sample)
        if _is_diseased(terms, term_name_to_id):
            disease_samples.add(sample)
        if _is_treated(terms, term_name_to_id):
            treated_samples.add(sample)
        if sample_to_type[sample] == 'cell line':
            cell_line_samples.add(sample)
        if filter_differentiated \
            and (sample_to_type[sample] == 'in vitro differentiated cells'
            or sample_to_type[sample] == 'induced pluripotent stem cell line'):
            differentiated_samples.add(sample)

    # Filter samples using filtering parameters
    if filter_poor:
        control_samples -= poor_samples
        case_samples -= poor_samples
    if filter_disease:
        control_samples -= disease_samples
    if filter_treated:
        control_samples -= treated_samples
        case_samples -= treated_samples
    if filter_cell_line:
        control_samples -= cell_line_samples
        case_samples -= cell_line_samples
    if filter_differentiated:
        control_samples -= differentiated_samples
        case_samples -= differentiated_samples

    # Partition each term into case and control samples
    control_term_set_to_samples = defaultdict(lambda: set())
    case_term_set_to_samples = defaultdict(lambda: set())
    for sample in case_samples:
        terms = sample_to_terms[sample]
        for term in terms:
            case_term_to_samples[term].add(sample)
        key_term_set = _create_key_terms(terms, term_name_to_id)
        case_term_set_to_samples[key_term_set].add(sample)
    for sample in control_samples:
        terms = sample_to_terms[sample]
        for term in terms:
            control_term_to_samples[term].add(sample)
        key_term_set = _create_key_terms(terms, term_name_to_id)
        control_term_set_to_samples[key_term_set].add(sample)

    # Search for confounding variables
    control_confound = set()
    case_confound = set()
    for term, samples in control_term_to_samples.items():
        if control_samples == control_term_to_samples[term]:
            control_confound.add(term)
    for term, samples in case_term_to_samples.items():
        if case_samples == case_term_to_samples[term]:
            case_confound.add(term)

    # Find common variables between case and control
    # identify tissue common variables
    tissue_intersections = set(control_term_set_to_samples.keys()) \
        & set(case_term_set_to_samples.keys())
    term_to_partition = {}
    for term_set in tissue_intersections:
        #term_id = term_name_to_id[term]
        #if ('UBERON' in term_id \
        #    and term != 'male organism' \
        #    and term != 'female organism' \
        #    and term != 'adult organism' \
        #    and term != 'organ') or \
        #    ('CL' in term_id \
        #    and term != 'cultured cell' \
        #    and term != 'cell' \
        #    and term != 'eukaryotic cell' \
        #    and term != 'animal cell' \
        #    and term != 'native cell'):
        #    tissue_intersections.add(term)
        term_to_partition[term_set] = {
            'case': list(case_term_set_to_samples[term_set]),
            'control': list(control_term_set_to_samples[term_set])
        }

    da = []
    for tissue_term in tissue_intersections:
        partition = term_to_partition[tissue_term]
        if by_run:
            for sample in partition['case']:
                if sample not in sample_to_runs:
                    continue
                for run in sample_to_runs[sample]:
                    da.append((
                        run,
                        sample_to_study[sample],
                        'case',
                        tissue_term,
                        sample in poor_samples,
                        sample in cell_line_samples,
                        sample in differentiated_samples,
                        sample in disease_samples
                    ))
            for sample in partition['control']:
                if sample not in sample_to_runs:
                    continue
                for run in sample_to_runs[sample]:
                    da.append((
                        run,
                        sample_to_study[sample],
                        'control',
                        tissue_term,
                        sample in poor_samples,
                        sample in cell_line_samples,
                        sample in differentiated_samples,
                        sample in disease_samples
                    ))
        else:
            for sample in partition['case']:
                da.append((
                    sample,
                    sample_to_study[sample],
                    'case',
                    tissue_term,
                    sample in poor_samples,
                    sample in cell_line_samples,
                    sample in differentiated_samples,
                    sample in disease_samples
                ))
            for sample in partition['control']:
                da.append((
                    sample,
                    sample_to_study[sample],
                    'control',
                    tissue_term,
                    sample in poor_samples,
                    sample in cell_line_samples,
                    sample in differentiated_samples,
                    sample in disease_samples
                ))
    if by_run:
        df = pd.DataFrame(data=da, columns=[
            'sample', 'project',
            'condition',
            'type', 'missing_metadata',
            'cell_line', 'differentiated',
            'diseased'
        ])
    else:
        df = pd.DataFrame(data=da, columns=[
            'sample', 'project', 
            'condition',
            'type', 'missing_metadata',
            'cell_line', 'differentiated',
            'diseased'
        ])
    return (
        df, 
        control_confound, 
        case_confound, 
        tissue_intersections
    )


def select_case_control_subset(df, case_control, term):
    if term is None:
        return list(df.loc[
            (df['condition'] == case_control), \
            'sample'
        ])
    else:
        return list(df.loc[
            (df['condition'] == case_control) \
            & (df['type'] == term), \
            'sample'
        ])


def create_barplot_most_common_coterms_match(
        df, view_cases, targ_term, sample_to_terms
    ):
    if targ_term is not None:
        targ_term = targ_term.replace(',', '\n')

    if view_cases:
        case_control = 'case'
    else:
        case_control = 'control'
    view_samples = select_case_control_subset(
        df, case_control, targ_term
    )    
    _create_barplot_most_common_coterms(
        view_samples, 
        sample_to_terms,
        skip_terms=set([targ_term])
    )

def create_barplot_most_common_coterms_series(
        val_to_samples, val, sample_to_terms
    ):
    if val in val_to_samples:
        view_samples = list(val_to_samples[val])
        print("Displaying data for %d sample with propert=%d" % (len(view_samples), val))
    else:
        print("Value {} was not found in the longitudinal query. Please try another query.".format(val))
    _create_barplot_most_common_coterms(
        view_samples, 
        sample_to_terms
    )   

def _create_barplot_most_common_coterms(
        view_samples, sample_to_terms, skip_terms=None 
    ):
    sample_to_terms = {
        sample: sample_to_terms[sample]
        for sample in view_samples
    }
    term_to_samples = defaultdict(lambda: set())
    for sample, terms in sample_to_terms.items():
        for term in terms:
            if skip_terms is None or term not in skip_terms:
                term_to_samples[term].add(sample)

    term_counts_df = pd.DataFrame(
        data = [
            (term, len(term_to_samples[term])/len(view_samples))
            for term in term_to_samples
        ],
        columns=['Term', 'Fraction of Samples']
    )

    term_counts_df = term_counts_df.sort_values(
        by='Fraction of Samples', 
        ascending=False
    )
    if term_counts_df.shape[0] > 20:
        term_counts_df = term_counts_df.iloc[:20]
    fig, ax = plt.subplots(
        1,
        1,
        sharey=False,
        figsize=(0.4*len(term_counts_df['Term'].unique()), 4)
        #figsize=(10,4)
    )
    sns.barplot(
        data=term_counts_df, 
        x='Term', 
        y='Fraction of Samples',
        color='#1057e5',
        ax=ax
    )
    ax.set_ylim((0.0, 1.0))
    plt.ylabel('Fraction of Samples', fontsize=12)
    plt.xlabel('Term', fontsize=12)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=90)
    plt.tight_layout()
    return term_counts_df


def create_pie_charts_matched(df, view_cases, targ_term, sample_to_terms):
    if view_cases:
        case_control = 'case'
    else:
        case_control = 'control'
    view_samples = select_case_control_subset(
        df, case_control, targ_term
    ) 
    _create_pie_charts(df, view_samples, sample_to_terms, skip_terms=set([targ_term]))


def create_pie_charts_series(
        df, val_to_samples, val, sample_to_terms
    ):
    if val in val_to_samples:
        view_samples = list(val_to_samples[val])
        print("Displaying most frequent co-occuring terms for %d sample with property = %d" % (len(view_samples), val))
    else:
        print("Value {} was not found in the longitudinal query. Please try another query.".format(val))
    _create_pie_charts(
        df, 
        view_samples,
        sample_to_terms
    )

def _create_pie_charts(df, view_samples, sample_to_terms, skip_terms=None):
    sample_to_terms = {
        sample: sample_to_terms[sample]
        for sample in view_samples
    }
    term_to_samples = defaultdict(lambda: set())
    for sample, terms in sample_to_terms.items():
        for term in terms:
            if skip_terms is None or term not in skip_terms:
                term_to_samples[term].add(sample)

    fig, axarr = plt.subplots(
        2,
        2,
        sharey=False,
        figsize=(8,8)
    )
    
    # Cell line pie chart
    n_cell_line = len(
        df.set_index('sample').loc[view_samples].loc[
            df.set_index('sample').loc[view_samples]['cell_line'] == True
        ]
    )
    n_no_cell_line = len(view_samples) - n_cell_line
    sizes = [n_cell_line, n_no_cell_line]
    labels = ['cell line', 'no cell line']
    axarr[0][0].set_title('Cell Line')
    patches, x, y = axarr[0][0].pie(
        sizes, 
        autopct=lambda p: '{:.1f}%'.format(round(p)) if p > 0 else '',
        shadow=False, startangle=90
    )
    axarr[0][0].legend(
        patches, 
        labels, 
        loc='upper right', 
        #bbox_to_anchor=(0.1, 1.),
        fontsize=12
    ) 

    # Sex pie chart
    n_male = len(term_to_samples['male organism'])
    n_female = len(term_to_samples['female organism'])
    n_unknown_sex = len(view_samples) - n_male - n_female
    sizes = [n_female, n_male, n_unknown_sex]
    labels = ['female', 'male', 'unknown']
    axarr[0][1].set_title('Sex')
    patches, x, y  = axarr[0][1].pie(
        sizes,
        autopct=lambda p: '{:.1f}%'.format(round(p)) if p > 0 else '',
        shadow=False, startangle=90
    )
    axarr[0][1].legend(
        patches, 
        labels, 
        loc='upper right', 
        #bbox_to_anchor=(0.1, 1.),
        fontsize=12
    )

    # Developmental stage pie chart
    n_adult = len(term_to_samples['adult organism'])
    n_embryo = len(term_to_samples['embryo']) + len(term_to_samples['embryonic cell'])
    n_unknown_dev = len(view_samples) - n_adult - n_embryo
    sizes = [n_adult, n_embryo, n_unknown_dev]
    labels = ['adult', 'embryonic', 'unknown']
    axarr[1][0].set_title('Developmental Stage')
    patches, x, y  = axarr[1][0].pie(
        sizes,
        autopct=lambda p: '{:.1f}%'.format(round(p)) if p > 0 else '',
        shadow=False, startangle=90
    )
    axarr[1][0].legend(
        patches,
        labels,
        loc='upper right',
        #bbox_to_anchor=(0.1, 1.),
        fontsize=12
    )

    # Developmental stage pie chart
    n_treat = len(term_to_samples['treatment'])
    n_no_treat = len(view_samples) - n_treat
    sizes = [n_treat, n_no_treat]
    labels = ['treatment', 'no treatment']
    axarr[1][1].set_title('Treatment')
    patches, x, y = axarr[1][1].pie(
        sizes,
        autopct=lambda p: '{:.1f}%'.format(round(p)) if p > 0 else '',
        shadow=False, startangle=90
    )
    axarr[1][1].legend(
        patches,
        labels,
        loc='upper right',
        #bbox_to_anchor=(0.1, 1.),
        fontsize=12
    )
    plt.tight_layout()

def create_series_plots(val_to_samples, target_property):
    df = pd.DataFrame(
        data=[
            (k, len(v)) 
            for k,v in val_to_samples.items()
        ], 
        columns=[
            target_property, 
            'Number of samples'
        ]
    )
    df.sort_values(target_property)
    plt.figure(figsize=(0.2*len(df),5.0))
    sns.barplot(
        x=target_property, 
        y="Number of samples", 
        data=df, 
        color='#1057e5'
    )
    plt.ylabel('Number of samples', fontsize=20)
    plt.xlabel(target_property, fontsize=20)
    plt.tight_layout()

    
def create_summary_plots(df):
    """
    Create bar-plots for the Case-Control Finder.
    """
    # The labels can be very long. We need to get
    # the maximum length to figure out a good height
    # for the plots.
    types = df['type'].unique()
    max_len = max([len(x) for x in types])
    grouped = df.groupby(by='type')
    da_n_studies = []
    for name, group in grouped:
        da_n_studies.append((
            name, 
            len(group.loc[(df['condition'] == 'case')]['project'].unique()), 
            'case'
        ))
        da_n_studies.append((
            name, 
            len(group.loc[(df['condition'] == 'control')]['project'].unique()), 
            'control'
        ))
    df_n_studies = pd.DataFrame(
        data=da_n_studies,
        columns=[
            'Tissue/Cell type',
            'Number of studies',
            'Condition'
        ]
    )
    fig, axarr = plt.subplots(
        1,
        2,
        sharey=False,
        figsize=(2*0.9*len(df_n_studies['Tissue/Cell type'].unique())+2.5, max_len/15+2.5)
    ) 
    sns.barplot(
        data=df_n_studies, 
        x='Tissue/Cell type', 
        y='Number of studies', 
        hue='Condition', 
        ax=axarr[0]
    )
    axarr[0].set_title('Number of studies\nper tissue/cell type')
    axarr[0].legend(
        loc='center left',
        bbox_to_anchor=(1, 0.5)
    )
    for p in axarr[0].patches:
        height = p.get_height()
        y_lim = axarr[0].get_ylim()[1]
        if height > 1000:
            x_offset = -0.1* p.get_width()
        else:
            x_offset = 0.1 * p.get_width()
        axarr[0].text(
            p.get_x() + x_offset,
            height + 0.015 * y_lim,
            '%d' % height,
            fontsize=9
        )
    axarr[0].set_ylim(0, axarr[0].get_ylim()[1] + 0.05*axarr[0].get_ylim()[1])
    plt.setp(axarr[0].xaxis.get_majorticklabels(), rotation=90)

    da_n_samples = []
    for name, group in grouped:
        da_n_samples.append((
            name, 
            len(group.loc[(df['condition'] == 'case')]), 
            'case'
        ))
        da_n_samples.append((
            name, 
            len(group.loc[(df['condition'] == 'control')]), 
            'control'
        ))
    df_n_samples = pd.DataFrame(
        data=da_n_samples,
        columns=[
            'Tissue/Cell type',
            'Number of samples',
            'Condition'
        ]
    )
    sns.barplot(
        data=df_n_samples, 
        x='Tissue/Cell type', y='Number of samples', 
        hue='Condition', 
        ax=axarr[1]
    )
    axarr[1].set_title('Number of samples\nper tissue/cell type')
    axarr[1].legend(
        loc='center left',
        bbox_to_anchor=(1, 0.5)
    )
    for p in axarr[1].patches:
        height = p.get_height()
        y_lim = axarr[1].get_ylim()[1]
        if height > 1000:
            x_offset = -0.1* p.get_width()
        else:
            x_offset = 0.1 * p.get_width()
        axarr[1].text(
            p.get_x() + x_offset,
            height + 0.015 * y_lim,
            '%d' % height,
            fontsize=9
        )
    axarr[1].set_ylim(0, axarr[1].get_ylim()[1] + 0.015*axarr[1].get_ylim()[1])
    plt.setp(axarr[1].xaxis.get_majorticklabels(), rotation=90)
    plt.tight_layout()
    
    
def load_metadata(available_data_f=None):
    """
    Load the SRA metadata.
    """
    sample_to_terms_f_json = './data/sample_to_terms.json'
    term_name_to_id_f = './data/term_name_to_id.json'
    sample_to_study_f = './data/sample_to_study.json'
    sample_to_real_value_terms_f = './data/sample_to_real_value.json'
    sample_to_runs_f = './data/sample_to_runs.json'
    sample_to_type_f = './data/sample_to_type.json'
    with open(sample_to_terms_f_json, 'r') as f:
        sample_to_terms = json.load(f)    
    with open(term_name_to_id_f, 'r') as f:
        term_name_to_id = json.load(f)
    with open(sample_to_type_f, 'r') as f:
        sample_to_type = json.load(f)
    with open(sample_to_study_f, 'r') as f:
        sample_to_study = json.load(f)
    with open(sample_to_real_value_terms_f, 'r') as f:
        sample_to_real_val = json.load(f)
    with open(sample_to_runs_f, 'r') as f:
        sample_to_runs = json.load(f)    
    if available_data_f:
        with open(available_data_f, 'r') as f:
            available = set(json.load(f))
        sample_to_terms = {
            k:v for k,v in sample_to_terms.items()  
            if k in available
        }    
    return (
        sample_to_terms,
        term_name_to_id,
        sample_to_type,
        sample_to_study,
        sample_to_runs,
        sample_to_real_val
    )

def main():
    with open('./data/experiment_to_terms.json', 'r') as f:
        sample_to_terms = json.load(f)

    with open('./data/term_name_to_id.json', 'r') as f:
        term_name_to_id = json.load(f)

    #with open('./data/experiments_in_hackathon_data.json', 'r') as f:
    with open('./data/my_available_exps.json', 'r') as f:
        available = set(json.load(f))

    with open('./data/experiment_to_type.json', 'r') as f:
        sample_to_type = json.load(f)

    with open('./data/experiment_to_study.json', 'r') as f:
        sample_to_study = json.load(f)

    with open('./data/experiment_to_real_value_terms.json', 'r') as f:
        sample_to_real_val = json.load(f)

    with open('./data/experiment_to_runs.json', 'r') as f:
        sample_to_runs = json.load(f)

    filter_available = False
    if filter_available:
        sample_to_terms = {
            k:v for k,v in sample_to_terms.items()
            if k in available
        }

    #term = 'glioblastoma multiforme' # A good one
    #term = 'cystic fibrosis' # okay

    """    
    term = 'blood'
    #term = 'brain'
    case, control = term_to_run(sample_to_terms, term)
    blacklist_terms = set(['disease', 'disease of cellular proliferation'])
    age_to_samples, df = series(term, 'age', sample_to_real_val, sample_to_terms,             
        sample_to_type, sample_to_study, term_name_to_id, blacklist_terms, 
        filter_poor=False, filter_cell_line=True, filter_differentiated=True,
        value_limit=100, target_unit=None
    )
    print(df)
    for age in sorted(age_to_samples.keys()):
        print("%d\t%d" % (age, len(age_to_samples[age])))
    """

    """
    r = match_case_to_controls(term, control, case, sample_to_terms, 
        sample_to_study, blacklist_terms, term_name_to_id, sample_to_type, 
        filter_poor=True, filter_cell_line=True, filter_differentiated=True,
        sample_to_runs=sample_to_runs)
    df = r[0]
    control_confound = r[1]
    case_confound = r[2]
    tissue_intersections = r[3]
    #df.to_csv('diabetes_case_control.csv')
    print(df)
    print('Tissue intersections: %s' % tissue_intersections)
    """

    #term = 'glioblastoma multiforme' # A good one
    #term = 'systemic lupus erythematosus'
    term = 'breast cancer'
    case, control = term_to_run(sample_to_terms, term)
    blacklist_terms = set(['disease', 'disease of cellular proliferation'])
    r = match_case_to_controls(term, control, case, sample_to_terms,
        sample_to_study, blacklist_terms, term_name_to_id, sample_to_type,
        filter_poor=True, filter_cell_line=True, filter_differentiated=True,
        sample_to_runs=sample_to_runs, by_run=False)
    df = r[0]
    control_confound = r[1]
    case_confound = r[2]
    tissue_intersections = r[3]
    #df.to_csv('glioblastoma.tsv', sep='\t')
    df.to_csv('breast_cancer.tsv', sep='\t')
    print(df)
    #print(df.loc[(df['type'] == 'brain')])
    #print(df.loc[(df['type'] == 'brain') & (df['condition'] == 'control')])
    #print('Tissue intersections: %s' % tissue_intersections)
    #print(select_case_control_experiment_set(df, 'case', 'blood'))

if __name__ == "__main__":
    main() 

    

