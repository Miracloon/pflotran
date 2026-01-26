# convert_ex_to_cell_ids.py

import sys
import re

if len(sys.argv) != 2:
    sys.exit("ERROR: A filename must be supplied in "+sys.argv[0])


pattern = r'\.[^.]+$' # finds the filename suffix
in_filename = sys.argv[1]
out_filename = re.sub(pattern,'.cell_ids',in_filename)

print('Converting .ex formatted file "{}" to a list of cell IDs in "{}".'.
      format(in_filename,out_filename))

fout = open(out_filename,'w')
i = 0
for line in open(in_filename,'r'):
    i += 1
    if i > 1: # skip the first line
        # grab the cell id which is the first integer of each line
        fout.write(line.split()[0].strip()+'\n')
fout.close

print('Done')
