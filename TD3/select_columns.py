#!/usr/bin/env python3

import sys
import csv

if __name__ == '__main__':
    if len(sys.argv) <= 2:
        print(f'Usage: {sys.argv[0]} csv_file index1 index2 ...')
        sys.exit(1)

    indices = [ int(sys.argv[k]) for k in range(2, len(sys.argv)) ]

    with open(sys.argv[1], newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')

        data = [ [field for field in row] for row in csv_reader ]

        for row in data:
            for k in range(len(row)):
                if k in indices or k == len(row)-1:
                    print(row[k], end='')
                    if k < len(row)-1:
                        print('\t', end='')
            print('\n',end='')
