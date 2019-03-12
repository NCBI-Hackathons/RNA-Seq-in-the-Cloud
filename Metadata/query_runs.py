import sys
import json

from functions import *


def main():

    json_file = sys.argv[1]
    term = sys.argv[2]

    out_file_in = '{}_term-{}.json'.format(json_file.split('.json')[0], term.replace(' ', '-'))
    out_file_notin = '{}_term-no-{}.json'.format(json_file.split('.json')[0], term.replace(' ', '-'))

    print('finding {} in {}'.format(term, json_file))
    print('printing runs with {} to {}'.format(term, out_file_in))
    print('printing runs without {} to {}'.format(term, out_file_notin))

    with open(json_file) as f:
        data_dict = json.load(f)

    runs_with_term, runs_without_term = term_to_run(data_dict, term)
    print(runs_with_term[0:10])
    print(runs_without_term[0:10])
    print("There are {} samples with the term '{}'".format(len(runs_with_term), term))
    print("There are {} samples without the term '{}'".format(len(runs_without_term), term))

    with open(out_file_in, 'w') as of:
        json.dump(runs_with_term, of, indent=True)
    with open(out_file_notin, 'w') as of:
        json.dump(runs_without_term, of, indent=True)


if __name__ == "__main__":
    main()
