import argparse
import requests
import sys

def parse_args(argv):
    """ Command line interface for the the module """
    parser = argparse.ArgumentParser()
    parser.add_argument("--species_list",
                        help="<Path to an existing FILE>",
                        dest="species_list",
                        action="store",
                        required=True)

    parser.add_argument("--sep",
                        help="<Field separator>",
                        dest="sep",
                        action="store",
                        required=True)

    parser.add_argument("--h_column",
                        help="<Column for proteome ID of species in healthy condition>",
                        dest="h_column",
                        action="store",
                        type=int,
                        required=True)

    parser.add_argument("--d_column",
                        help="<Column for proteome ID of species in diseased condition>",
                        dest="d_column",
                        action="store",
                        type=int,
                        required=True)

    parser.add_argument("--output_h",
                        help="<Path to the output files (healthy condition)>",
                        dest="output_h",
                        action="store",
                        required=True)

    parser.add_argument("--output_d",
                        help="<Path to the output files (diseased condition)>",
                        dest="output_d",
                        action="store",
                        required=True)

    results = parser.parse_args(argv)
    return results


# Creating a list of proteome ids to download from Uniprot
def read_proteome_ids(file, separator, id_column):
    with open(file) as proteome_list:
        proteome_list.readline()
        proteomeids = []
        for line in proteome_list:
            line = line.strip().split(separator)

            #Python starts to count by 0, therefore the user-provided column number should be decreased by one
            id_column = int(id_column) - 1
            proteomeids.append(line[id_column])
    return proteomeids


def download_proteome(id):
    url = 'https://www.uniprot.org/uniprot/?query=proteome:' + id + '&format=tab&columns=id,database(Pfam)'
    #proteome_url = "http://www.uniprot.org/uniprot/?sort=&desc=&compress=no&query=proteome:" + id + "&fil=&force=no&preview=true&format=fasta"
    r = requests.get(url)
    return r.text

def write_output(proteomeids, output_file):
    with open(output_file, 'w') as output:
        for id in proteomeids:
            print(id)
            output.write(download_proteome(id))


def main(argv):
    """ Main method and logic """

    # Read args
    args = parse_args(argv)

    # Download proteomes of healthy species
    proteome_ids = read_proteome_ids(args.species_list, args.sep, args.h_column)
    write_output(proteome_ids, args.output_h)

    # Download proteomes of species in periodontitis
    proteome_ids = read_proteome_ids(args.species_list, args.sep, args.d_column)
    write_output(proteome_ids, args.output_d)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


