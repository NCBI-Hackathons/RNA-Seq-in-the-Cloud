import sys
import json


def term_to_run(data_dict, term):

    runs_with_term = [k for k, v in data_dict.items() if term in v]
    return runs_with_term


def main():

    json_file = sys.argv[1]
    term = sys.argv[2]

    out_file = '{}_term-{}.json'.format(json_file.split('.json')[0], term)

    print('finding {} in {}'.format(term, json_file))
    print('printing runs with {} to {}'.format(term, out_file))

    with open(json_file) as f:
        data_dict = json.load(f)

    run_to_term_dict = data_dict['run_to_term']
    runs_with_term = term_to_run(run_to_term_dict, term)
    print(runs_with_term)

    with open(out_file, 'w') as of:
        json.dump(runs_with_term, of, indent=True)


if __name__ == "__main__":
    main()
