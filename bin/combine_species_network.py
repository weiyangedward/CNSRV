import sys
import argparse

def combine_networks(network_files):
    for netw in network_files:
        with open(netw, 'r') as F:
            spe = F.readline().strip()
            num_genes = F.readline()
            for line in F:
                g1, g2, edge = line.strip().split()
                print "%s\t%s\t%s\t%s" % (spe, g1, g2, edge)

            F.close()


def main():
    parser = argparse.ArgumentParser(description='combine different species coexpression network into one file.')
    parser.add_argument('network_files', type=str, nargs='+', help='network pair files')
    args = parser.parse_args()

    combine_networks(args.network_files)


if __name__ == '__main__':
    main()