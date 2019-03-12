from collections import defaultdict
import json

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
        if 'UBERON' in term_id and term_id != 'male organism' and term_id != 'female organism':
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
 

def match_case_to_controls(term, control_samples, case_samples, sample_to_terms, blacklist_terms, term_name_to_id, filter_poor=True, filter_cell_line=True):
    filtered = set()
    for sample in control_samples:
        if len(blacklist_terms & set(sample_to_terms[sample])) == 0:
            filtered.add(sample)
    control_samples = filtered

    control_term_to_samples = defaultdict(lambda: set())
    case_term_to_samples = defaultdict(lambda: set())
    poor_case = set()
    poor_control = set()
    cell_line_case = set()
    cell_line_control = set()
    for sample in control_samples:
        terms = sample_to_terms[sample]
        if filter_poor and _is_poor_quality(terms, term_name_to_id):
            print("Sample %s is poor" % sample)
            continue
        if filter_cell_line and _is_cell_line(terms, term_name_to_id):
            continue
        for term in terms:
            control_term_to_samples[term].add(sample)
    
    for sample in case_samples:
        terms = sample_to_terms[sample]
        terms = sample_to_terms[sample]
        if filter_poor and _is_poor_quality(terms, term_name_to_id):
            continue
        if filter_cell_line and _is_cell_line(terms, term_name_to_id):
            continue
        for term in terms:
            case_term_to_samples[term].add(sample)

    control_confound = set()
    case_confound = set()
    for term, samples in control_term_to_samples.items():
        if control_samples == control_term_to_samples[term]:
            control_confound.add(term)
    for term, samples in case_term_to_samples.items():
        if case_samples == case_term_to_samples[term]:
            case_confound.add(term)

    intersect_terms = set(control_term_to_samples.keys()) \
        & set(case_term_to_samples.keys())

    print("Term intersection: %s" % intersect_terms)
   
    term_to_partition = {}
    tissue_intersections = set()
    for term in intersect_terms:
        term_id = term_name_to_id[term]
        if 'UBERON' in term_id and term != 'male organism' and term != 'female organism':
            tissue_intersections.add(term)
        term_to_partition[term] = {
            'case': list(case_term_to_samples[term]),
            'control': list(control_term_to_samples[term])
        }

    return (
        term_to_partition,
        control_confound,
        case_confound,
        tissue_intersections
    )


def generate_output(term_to_partition, tissue_intersections):
    table = []
    for tissue_term in tissue_intersections:
        partition = term_to_partition[tissue_term]
        for sample in partition['case']:
            table.append((
                sample,
                'case',
                tissue_term
            ))
        for sample in partition['control']:
            table.append([
                sample,
                'control',
                tissue_term
            ])
    with open('output.txt', 'w') as f:
        for entry in table:
            f.write("%s\n" % "\t".join(entry))

def main():
    with open('./data/experiment_to_terms.json', 'r') as f:
        sample_to_terms = json.load(f)

    with open('./data/term_name_to_id.json', 'r') as f:
        term_name_to_id = json.load(f)

    with open('./data/experiments_in_hackathon_data.json', 'r') as f:
        available = set(json.load(f))


    filter_available = False
    if filter_available:
        sample_to_terms = {
            k:v for k,v in sample_to_terms.items()
            if k in available
        }

    term = 'glioblastoma multiforme'
    case, control = term_to_run(sample_to_terms, term)

    blacklist_terms = set(['disease', 'disease of cellular proliferation'])
    r = match_case_to_controls(term, control, case, sample_to_terms, 
        blacklist_terms, term_name_to_id, filter_poor=True, 
        filter_cell_line=True)
    term_to_partition = r[0]
    control_confound = r[1] 
    case_confound = r[2] 
    tissue_intersections = r[3] 
    
    print(json.dumps({
        k: {k1: len(v1) for k1,v1 in v.items()} for k,v in term_to_partition.items()
    }, indent=True))

    print("Terms that confound the control: %s" % control_confound)
    print("Terms that confound the case: %s" % case_confound)
    print("Tissue intersections: %s" % tissue_intersections)

    generate_output(term_to_partition, tissue_intersections)

if __name__ == "__main__":
    main() 

    

