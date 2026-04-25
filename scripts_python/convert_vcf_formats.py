import sys
import os

def parse_bsdp_acgt(bsdp_str, ref, alt):
    """
    Parses BSDP. Supports:
    1. A,C,G,T format (4 values)
    2. Legacy Alt,Ref format (2 values)
    Returns (ref_count, alt_count).
    """
    if bsdp_str == "." or bsdp_str == ".,.":
        return ".", "."

    counts = bsdp_str.split(',')
    
    if len(counts) == 4:
        bases = ['A', 'C', 'G', 'T']
        try:
            ref_idx = bases.index(ref)
            alt_idx = bases.index(alt.split(',')[0])
            return counts[ref_idx], counts[alt_idx]
        except ValueError:
            return 0, 0
    elif len(counts) >= 2:
        # Legacy NGSEP: Alt,Ref
        return counts[1], counts[0]
    
    return 0, 0

def convert_vcf(input_vcf, target_format):
    out_vcf = input_vcf.replace('.vcf', f'_{target_format}.vcf')
    
    with open(input_vcf, 'r') as fin, open(out_vcf, 'w') as fout:
        for line in fin:
            if line.startswith('##source='):
                fout.write(f'##source={target_format}\n')
            elif line.startswith('##FORMAT=<ID=BSDP'):
                if target_format == 'GATK':
                    fout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
                elif target_format == 'freebayes':
                    fout.write('##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count">\n')
                    fout.write('##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">\n')
            elif line.startswith('#'):
                fout.write(line)
            else:
                parts = line.strip().split('\t')
                ref = parts[3]
                alt = parts[4]
                format_field = parts[8]
                
                format_tokens = format_field.split(':')
                try:
                    bsdp_idx = format_tokens.index('BSDP')
                except ValueError:
                    fout.write(line) # No BSDP found
                    continue
                
                if target_format == 'GATK':
                    format_tokens[bsdp_idx] = 'AD'
                elif target_format == 'freebayes':
                    format_tokens[bsdp_idx] = 'RO:AO'
                    
                parts[8] = ':'.join(format_tokens)
                
                for i in range(9, len(parts)):
                    sample_data = parts[i].split(':')
                    if len(sample_data) > bsdp_idx:
                        bsdp_str = sample_data[bsdp_idx]
                        ref_c, alt_c = parse_bsdp_acgt(bsdp_str, ref, alt)
                        
                        if target_format == 'GATK':
                            sample_data[bsdp_idx] = f'{ref_c},{alt_c}'
                        elif target_format == 'freebayes':
                            sample_data[bsdp_idx] = f'{ref_c}:{alt_c}'
                                
                    parts[i] = ':'.join(sample_data)
                    
                fout.write('\t'.join(parts) + '\n')
                
    return out_vcf

if __name__ == '__main__':
    input_file = sys.argv[1]
    convert_vcf(input_file, 'GATK')
    convert_vcf(input_file, 'freebayes')
