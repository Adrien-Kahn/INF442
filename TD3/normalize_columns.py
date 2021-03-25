#!/usr/bin/env python3

import sys
import csv

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f'Usage: {sys.argv[0]} csv_file')
        sys.exit(1)

    with open(sys.argv[1], newline='') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')

        data = [ [field for field in row] for row in csv_reader ]

        first_row = data[0]

        column_min = [float('+inf') for j in range(len(first_row)-1)]
        column_max = [float('-inf') for j in range(len(first_row)-1)]

        for row in data[1:]:
            for j in range(len(row)-1):
                column_min[j] = min(column_min[j], float(row[j]))
                column_max[j] = max(column_max[j], float(row[j]))

        for row in data[1:]:
            for j in range(len(row)-1):
                row[j] = (float(row[j])-column_min[j]) / (column_max[j]-column_min[j])

        for row in data:
            for j in range(len(row)):
                print(row[j], end='')
                if j < len(row)-1:
                    print('\t', end='')
            print('\n',end='')
