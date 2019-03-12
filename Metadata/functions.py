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

def match_case_to_controls(term, control_samples, case_samples, sample_to_terms):
    extra_terms_with = set()
    control_term_to_samples = defaultdict(lambda: set())
    case_term_to_samples = defaultdict(lambda: set())
    for sample in control_samples:
        for term in sample_to_terms[sample]:
            control_term_to_samples[term].add(sample)
    for sample in case_samples:
        for term in sample_to_terms[sample]:
            case_term_to_samples[term].add(sample)

    intersect_terms = set(control_term_to_samples.keys()) \
        & set(case_term_to_samples.keys())
    
    term_to_partition = {}
    for term in intersect_terms:
        term_to_partition[term] = {
            'case': list(case_term_to_samples[term]),
            'control': list(control_term_to_samples[term])
        }
    return term_to_partition



def main():
    with open('./data/experiment_to_terms.json', 'r') as f:
        sample_to_terms = json.load(f)
    term = 'breast cancer'
    case, control = term_to_run(sample_to_terms, term)

    #print case
    #print control
     
    print(json.dumps({
        k: {k1: len(v1) for k1,v1 in v.items()} for k,v in match_case_to_controls(term, control, case, sample_to_terms).items()
    }, indent=True))

    #match_case_to_controls(term, control, case, sample_to_terms)

if __name__ == "__main__":
    main() 

    

