from optparse import OptionParser
import json

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    #parser.add_option("-a", "--a_descrip", action="store_true", help="This is a flat")
    #parser.add_option("-b", "--b_descrip", help="This is an argument")
    (options, args) = parser.parse_args()

    with open('transfected_sample_raw_metadata.json', 'r') as f:
        sample_to_key_to_val = json.load(f)

    gene_id_to_symbol = {}
    gene_id_to_name = {}
    with open('genes.tsv', 'r') as f:
        for i,l in enumerate(f):
            if i == 0:
                continue
            toks = l.split()       
            g_id = toks[0]
            symbol = toks[1]
            gene_id_to_symbol[g_id] = symbol
            if len(toks) == 3:
                name = toks[2]
                gene_id_to_name[g_id] = name 
            
    


if __name__ == "__main__":
    main()
