import random
import json


def main():

    data_dict = {}
    run_to_term = {}
    run_to_study = {}

    runs = range(100)
    terms = 15
    studies = 9

    for run in runs:
        num_terms = random.randint(1, 5)
        run_to_term[run] = str(random.sample(range(1, terms), num_terms))
        run_to_study[run] = random.randint(1, studies)

    data_dict['run_to_term'] = run_to_term
    data_dict['run_to_study'] = run_to_study

    print(data_dict)

    with open('result.json', 'w') as f:
        json.dump(data_dict, f, indent=True)


if __name__ == "__main__":
    main()