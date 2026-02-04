import sys

def read_cells(filename):
    i = 0
    cell_ids = []
    for line in open(filename,'r'):
        i += 1
        if i == 1:
             continue
        cell_ids.append(int(line.split()[0]))
    return cell_ids

duplicate_cells = set([0,1]) & set([1,2])
print(duplicate_cells)

duplicate_cells = set(read_cells(sys.argv[1])) & set(read_cells(sys.argv[2]))
print(duplicate_cells)
print(f'Number of duplicates: {len(duplicate_cells)}')

print('Done')
