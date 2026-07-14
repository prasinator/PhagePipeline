import sys
import os

def format_fasta(input_path, output_path, phage_id="phage"):
    with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
        count = 0
        seq_lines = []
        header = None
        
        def write_record(h, s_lines):
            nonlocal count
            if not h:
                return
            count += 1
            seq_str = "".join(s_lines).replace("\n", "").replace(" ", "")
            
            # Check if it already has the ref| format
            if h.startswith(">ref|"):
                f_out.write(h + "\n" + "\n".join(s_lines) + "\n")
                return
                
            # Parse header
            parts = h[1:].split()
            id_part = parts[0]
            
            # Parse gene call header: >ID # START # END # STRAND # ...
            start = count * 1000
            end = start + len(seq_str) * 3
            strand_sign = "+"
            product = "unknown"
            
            if len(parts) >= 7 and parts[1] == "#" and parts[3] == "#" and parts[5] == "#":
                try:
                    start = int(parts[2])
                    end = int(parts[4])
                    strand = int(parts[6])
                    strand_sign = "+" if strand >= 0 else "-"
                except ValueError:
                    pass
            
            # Construct standard header
            new_header = f">ref|{id_part}|{count}|[{start}:{end}]({strand_sign})|{product}|{phage_id}"
            f_out.write(new_header + "\n" + "\n".join(s_lines) + "\n")
            
        for line in f_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                write_record(header, seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        write_record(header, seq_lines)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: format_fasta.py <input.faa> <output.faa> [phage_id]")
        sys.exit(1)
    
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    phage_id = sys.argv[3] if len(sys.argv) > 3 else "phage"
    format_fasta(in_path, out_path, phage_id)
