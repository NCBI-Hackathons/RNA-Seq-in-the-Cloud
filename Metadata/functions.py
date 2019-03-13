from collections import defaultdict
import json
import pandas as pd

def term_to_run(sample_to_terms, term):

    runs_with_term = []
    runs_without_term = []
    for samples, terms in sample_to_terms.items():
        if term in terms:
            runs_with_term.append(samples)
        else:
            runs_without_term.append(samples)

    assert set(runs_with_term) != set(runs_without_term), 'Oops! The set of samples with the term and without the term overlap!'

    return runs_with_term, runs_without_term

def _is_poor_quality(terms, term_name_to_id):
    found_tissue = False
    found_cell_line = False
    found_cell_type = False
    for term in terms:
        term_id = term_name_to_id[term]
        if 'UBERON' in term_id and term != 'male organism' \
            and term != 'female organism' and term != 'adult organism':
            found_tissue = True
        elif 'CVCL' in term_id:
            found_cell_line = True
        elif 'CL' in term_id and term != 'cultured cell':
            found_cell_type = True
    return not (found_tissue or found_cell_line or found_cell_type)
       
def _is_cell_line(terms, term_name_to_id):
    for term in terms:
        if 'CVCL' in term_name_to_id[term]:
            return True
    return False


def series(term, target_property, sample_to_real_val, sample_to_terms, sample_to_type, 
        sample_to_study, term_name_to_id, blacklist_terms, 
        filter_poor=True, filter_cell_line=True, filter_differentiated=True,
        target_unit=None, value_limit=None, skip_missing_unit=False
    ):
    age_to_samples = defaultdict(lambda: set())
    for sample, real_val_infos in sample_to_real_val.items():
        if sample not in sample_to_terms:
            continue
        for real_val_info in real_val_infos:
            property_ = real_val_info['property']
            unit = real_val_info['unit']
            value = int(real_val_info['value'])
            if target_unit and unit != target_unit:
                continue
            poor_samples = set()
            cell_line_samples = set()
            differentiated_samples = set()
            if property_ == target_property:
                if value_limit and value > value_limit:
                    continue
                terms = sample_to_terms[sample]
                if len(blacklist_terms & set(sample_to_terms[sample])) > 0:
                    continue
                if term in terms:
                    age_to_samples[value].add(sample)
                if _is_poor_quality(terms, term_name_to_id):
                    poor_samples.add(sample)
                    if filter_poor:
                        continue
                if sample_to_type[sample] == 'cell line':
                    cell_line_samples.add(sample)
                    if filter_cell_line:
                        continue
                if filter_differentiated \
                    and (sample_to_type[sample] == 'in vitro differentiated cells'
                    or sample_to_type[sample] == 'induced pluripotent stem cell line'):
                    differentiated_samples.add(sample)
                    if filter_differentiated:
                        continue

    da = []
    for age in sorted(age_to_samples.keys()):
        for sample in age_to_samples[age]:
            da.append((
                sample,
                sample_to_study[sample],
                age,
                sample in poor_samples,
                sample in cell_line_samples,
                sample in differentiated_samples
            ))
    df = pd.DataFrame(data=da, columns=[
        'experiment_accession', 'study_accession',
        'age', 'missing_metadata',
        'cell_line', 'differentiated'    
    ])

    return age_to_samples, df

def match_case_to_controls(term, control_samples, case_samples, sample_to_terms, 
    sample_to_study, blacklist_terms, term_name_to_id, sample_to_type, 
    filter_poor=True, filter_cell_line=True, filter_differentiated=True, 
    by_run=False, sample_to_runs=None):
    filtered = set()
    control_samples = set(control_samples)
    case_samples = set(case_samples)

    for sample in control_samples:
        if len(blacklist_terms & set(sample_to_terms[sample])) == 0:
            filtered.add(sample)
    control_samples = filtered

    control_term_to_samples = defaultdict(lambda: set())
    case_term_to_samples = defaultdict(lambda: set())

    # Identify poor quality, in vitro differentiated,
    # and cell line samples
    poor_samples = set()
    cell_line_samples = set()
    differentiated_samples = set()
    for sample in set(control_samples) | set(case_samples):
        terms = sample_to_terms[sample]
        if _is_poor_quality(terms, term_name_to_id):
            poor_samples.add(sample)
        if sample_to_type[sample] == 'cell line':
        #if _is_cell_line(terms, term_name_to_id):
            cell_line_samples.add(sample)
        if filter_differentiated \
            and (sample_to_type[sample] == 'in vitro differentiated cells'
            or sample_to_type[sample] == 'induced pluripotent stem cell line'):
            differentiated_samples.add(sample)

    # Filter samples using filtering parameters
    if filter_poor:
        control_samples -= poor_samples
        case_samples -= poor_samples
    if filter_cell_line:
        control_samples -= cell_line_samples
        case_samples -= cell_line_samples
    if filter_differentiated:
        control_samples -= differentiated_samples
        case_samples -= differentiated_samples

    # Partition each term into case and control samples
    for sample in control_samples:
        terms = sample_to_terms[sample]
        for term in terms:
            control_term_to_samples[term].add(sample)
    for sample in case_samples:
        terms = sample_to_terms[sample]
        print(terms)
        for term in terms:
            case_term_to_samples[term].add(sample)    

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
    intersect_terms = set(control_term_to_samples.keys()) \
        & set(case_term_to_samples.keys())
    print('Intersect terms: %s' % intersect_terms)
    term_to_partition = {}
    tissue_intersections = set()
    for term in intersect_terms:
        term_id = term_name_to_id[term]
        if 'UBERON' in term_id \
            and term != 'male organism' \
            and term != 'female organism' \
            and term != 'adult organism':
            tissue_intersections.add(term)
        term_to_partition[term] = {
            'case': list(case_term_to_samples[term]),
            'control': list(control_term_to_samples[term])
        }

    da = []
    for tissue_term in tissue_intersections:
        partition = term_to_partition[tissue_term]
        if by_run:
            for sample in partition['case']:
                for run in sample_to_runs[sample]:
                    da.append((
                        run,
                        sample_to_study[sample],
                        'case',
                        tissue_term,
                        sample in poor_samples,
                        sample in cell_line_samples,
                        sample in differentiated_samples
                    ))
            for sample in partition['control']:
                for run in sample_to_runs[sample]:
                    da.append((
                        run,
                        sample_to_study[sample],
                        'control',
                        tissue_term,
                        sample in poor_samples,
                        sample in cell_line_samples,
                        sample in differentiated_samples
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
                    sample in differentiated_samples
                ))
            for sample in partition['control']:
                da.append((
                    sample,
                    sample_to_study[sample],
                    'control',
                    tissue_term,
                    sample in poor_samples,
                    sample in cell_line_samples,
                    sample in differentiated_samples
                ))
    if by_run:
        df = pd.DataFrame(data=da, columns=[
            'Run', 'project',
            'condition',
            'type', 'missing_metadata',
            'cell_line', 'differentiated'
        ])
    else:
        df = pd.DataFrame(data=da, columns=[
            'experiment', 'project', 
            'condition',
            'type', 'missing_metadata',
            'cell_line', 'differentiated'
        ])
    return df, control_confound, case_confound, tissue_intersections

def select_case_control_experiment_set(df, case_control, term):
    return list(df.loc[(df['condition'] == case_control) & (df['type'] == term), 'experiment_accession'])

def main():
    with open('./data/experiment_to_terms.json', 'r') as f:
        sample_to_terms = json.load(f)

    with open('./data/term_name_to_id.json', 'r') as f:
        term_name_to_id = json.load(f)

    with open('./data/experiments_in_hackathon_data.json', 'r') as f:
        available = set(json.load(f))

    with open('./data/experiment_to_type.json', 'r') as f:
        sample_to_type = json.load(f)

    with open('./data/experiment_to_study.json', 'r') as f:
        sample_to_study = json.load(f)

    with open('./data/experiment_to_real_value_terms.json', 'r') as f:
        sample_to_real_val = json.load(f)

    with open('./data/experiment_to_runs.json', 'r') as f:
        sample_to_runs = json.load(f)

    filter_available = True
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

    term = 'glioblastoma multiforme' # A good one
    case, control = term_to_run(sample_to_terms, term)
    blacklist_terms = set(['disease', 'disease of cellular proliferation'])
    r = match_case_to_controls(term, control, case, sample_to_terms,
        sample_to_study, blacklist_terms, term_name_to_id, sample_to_type,
        filter_poor=True, filter_cell_line=True, filter_differentiated=True,
        sample_to_runs=sample_to_runs, by_run=True)
    df = r[0]
    control_confound = r[1]
    case_confound = r[2]
    tissue_intersections = r[3]
    df.to_csv('glioblastoma_multiforme_case_control.csv')
    print(df)
    #print(df.loc[(df['type'] == 'brain')])
    #print('Tissue intersections: %s' % tissue_intersections)
    #print(select_case_control_experiment_set(df, 'case', 'blood'))

if __name__ == "__main__":
    main() 

    

