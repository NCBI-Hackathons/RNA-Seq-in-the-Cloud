def term_to_run(data_dict, term):

    runs_with_term = []
    runs_without_term = []
    for samples, terms in data_dict.items():
        if term in terms:
            runs_with_term.append(samples)
        else:
            runs_without_term.append(samples)

    assert set(runs_with_term) != set(runs_without_term), 'Oops! The set of samples with the term and without the term overlap!'

    return runs_with_term, runs_without_term