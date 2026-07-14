import os
import sys

input_file = sys.argv[1]
out_dir = sys.argv[2]

os.makedirs(out_dir, exist_ok=True)

with open(input_file, 'r') as f:
    seq_id = ''
    seq_data = []
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if seq_id:
                with open(os.path.join(out_dir, f"{seq_id}.fasta"), 'w') as out:
                    out.write(f'>{seq_id}\n' + '\n'.join(seq_data) + '\n')
            seq_id = line[1:].split()[0]
            seq_data = []
        else:
            seq_data.append(line)
    if seq_id:
        with open(os.path.join(out_dir, f"{seq_id}.fasta"), 'w') as out:
            out.write(f'>{seq_id}\n' + '\n'.join(seq_data) + '\n')
