## Script to split the output files of the code into separate files based on metallicity
Zvals = [-2,-1,0]

def split_by_metallicity(input_file, output_prefix):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    header = lines[0]  # <-- This is the added line to capture the header
    blocks = []
    current_block = []

    for line in lines:
        if line.strip() == "":
            if current_block:
                blocks.append(current_block)
                current_block = []
        else:
            current_block.append(line)
    
    # Add the final block if file doesn't end with a blank line
    if current_block:
        blocks.append(current_block)

    for i, block in enumerate(blocks):
        output_file = "{}_Z{}.dat".format(output_prefix, Zvals[i])
        with open(output_file, 'w') as f_out:
            f_out.writelines([header] + block)  # <-- This line ensures header is included
        print(f"Written {output_file} with {len(block)} lines")

# Example usage
split_by_metallicity("fort.22", "CollapseAB")
split_by_metallicity("fort.31", "CollapseCools")
split_by_metallicity("fort.911", "CollapseHeats")