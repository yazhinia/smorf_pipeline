import sys
from collections import defaultdict
import re

"""
Merge GFF files by seqid, start, end, strand
In-frame overlapping or within smOFRs with long ORFs are retained 
"""

def normalize_seqname(seqname):
    """
    Strips prefixes/suffixes to extract the core contig identifier.

    contig_1       -> 1 bakta
    1              -> 1 Genelm
    filename-1     -> 1 Smorfinder
    filename-1_2   -> 1 gmscmapper (drop the ORF suffix after last underscore)
    1_2            -> 1 gmscmapper without a filename prefix (drop the ORF suffix)
    """

    # Current HRGM2 naming: "<genome>_c<NNN>" (bakta, smorfinder) and
    # "<genome>_c<NNN>_<orf_index>" (gmscmapper). The contig id is the c-token.
    m = re.match(r'^.+_(c\d+)(?:_\d+)?$', seqname)
    if m:
        return m.group(1)

    # gmscmapper without a filename prefix: bare "<contig>_<orf_index>".
    # Both sides are numeric here, unlike bakta's "contig_1" (word prefix),
    # so the trailing part is an ORF index to drop, not the contig id.
    if re.match(r'^\d+_\d+$', seqname):
        return seqname.split('_')[0]

    # Tool 1: contig_1 — take part after underscore
    if '-' in seqname and re.search(r'_\d+$', seqname):
        seqname = re.sub(r'_\d+$', '', seqname)   # filename-1_2 -> filename-1

    if 'GENOME' in seqname:
        if '_' in seqname:
            seqname = seqname.split('_')[1]
    # Handle NCBI accessions with ENA| prefix
    if 'ENA|' in seqname:
        seqname = seqname.split('.')[-1]
    # Tool 3 & 4: has hyphen — take after last hyphen
    if '-' in seqname:
        seqname = seqname.split('-')[-1]           # filename-1 -> 1

    # Tool 1: contig_1 — take after underscore
    elif '_' in seqname:
        seqname = seqname.split('_')[-1]           # contig_1 -> 1
    
    # Tool 2: already bare — "1" stays "1" 
    return seqname

def natural_sort_key(seq_str):
    
    """Sort numeric strings numerically, non-numeric alphabetically."""
    return int(seq_str) if seq_str.isdigit() else seq_str

def update_id_with_contig(attrs, contig_id):
    """Append contig_id to the ID field in the attributes column."""
    attr_parts = attrs.rstrip('\n').split(';')
    new_parts = []
    for kv in attr_parts:
        if kv.startswith('ID='):
            old_id = kv.split('=', 1)[1]
            new_parts.append(f"ID={contig_id}:{old_id}")
        else:
            new_parts.append(kv)
    return ';'.join(new_parts)

def get_source(cols):
    """Return the source tool from column 1 (index 1)."""
    return cols[1]  # e.g. 'Pyrodigal' or 'GeneLM'

# GFF column-2 source string -> tool name used in the tools= attribute.
# Keys are matched lowercased: bakta writes 'Pyrodigal' for ab-initio CDS but 'Bakta' for the
# CDS its own sORF/expert-system pass adds, and both have to fold to the same tool name.
SOURCE_TO_TOOL = {
    'pyrodigal': 'bakta',
    'bakta': 'bakta',
    'prodigal_v2.6.3': 'smorfinder',
    'gmscmapper': 'gmscmapper',
    'genelm': 'genelm',
}

def source_to_tool(source):
    """Map a GFF column-2 source string to a canonical tool name.

    Unknown sources are lowercased rather than returned verbatim. Passing the raw string
    through is what let bakta's 'Bakta'-sourced CDS become a phantom tool alongside 'bakta',
    splitting ~0.2% of records into tool combinations that do not exist.
    """
    return SOURCE_TO_TOOL.get(source.lower(), source.lower())

def nt_to_aa(nt_length):
    return nt_length / 3

def get_length(cols):
    return int(cols[4]) - int(cols[3]) + 1

def overlaps(s1, e1, s2, e2):
    """Check if two CDS intervals overlap at all."""
    return s1 <= e2 and s2 <= e1

def is_contained(s1, e1, s2, e2):
    """Check if one interval is fully contained within the other."""
    return (s1 >= s2 and e1 <= e2) or (s2 >= s1 and e2 <= e1)

def shared_boundary(s1, e1, s2, e2):
    """Check if intervals share a start or end coordinate."""
    return s1 == s2 or e1 == e2

def same_stop(s1, e1, s2, e2, strand):
    """Whether two same-strand CDS share a stop codon.

    The stop sits at the end coordinate on '+' and at the start coordinate on '-'.
    Sharing it means one ORF called with two different start codons; differing on it
    means two different ORFs that merely overlap.
    """
    return e1 == e2 if strand == '+' else s1 == s2

def _absorb(entries, survivor, dropped, same_gene):
    """Drop `dropped`, folding its tools into `survivor` only if it is the same gene.

    Every drop site goes through here. When the two share a stop codon they are one ORF
    with two start calls, and discarding the loser's tag is what made support_n undercount
    agreement -- for smORFs, start-codon disagreement is the common case. When they do not
    share a stop they are different ORFs, and carrying the tag across would invent support
    for a call that tool never made.
    """
    if same_gene:
        entries[survivor]['tools'] |= entries[dropped]['tools']
    entries[dropped]['keep'] = False

def resolve_overlaps(cds_records):
    AA_CUTOFF = 100

    entries = [{'cols': c, 'len_nt': l, 'keep': True,
                'tools': {source_to_tool(get_source(c))}} for c, l in cds_records]

    n = len(entries)
    for i in range(n):
        if not entries[i]['keep']:
            continue
        s1   = int(entries[i]['cols'][3])
        e1   = int(entries[i]['cols'][4])
        aa1  = nt_to_aa(entries[i]['len_nt'])
        tool1 = source_to_tool(get_source(entries[i]['cols']))
        strand1 = entries[i]['cols'][6]
        
        for j in range(i + 1, n):
            if not entries[i]['keep']:
                break
            if not entries[j]['keep']:
                continue
            s2   = int(entries[j]['cols'][3])
            e2   = int(entries[j]['cols'][4])
            aa2  = nt_to_aa(entries[j]['len_nt'])
            tool2 = source_to_tool(get_source(entries[j]['cols']))
            strand2 = entries[j]['cols'][6]
            
            if strand1 != strand2:
                continue
                
            if not overlaps(s1, e1, s2, e2):
                continue
                
            # Sharing a stop codon is what makes two records the same gene, and so what
            # decides whether a dropped record's tool tag transfers to the survivor.
            same_gene = same_stop(s1, e1, s2, e2, strand1)

            # Rule 1: exact duplicate — keep first, drop second.
            # Same coordinates + strand = the same ORF called by another tool.
            if s1 == s2 and e1 == e2:
                _absorb(entries, i, j, same_gene)
                continue
            
            # Rule 2: partial overlap — neither boundary shared, not contained
            # Keep both — legitimate overlapping genes in compact genomes
            if not shared_boundary(s1, e1, s2, e2) and not is_contained(s1, e1, s2, e2):
                continue  # keep both, skip further checks
                
            aa_diff  = abs(aa1 - aa2)
            larger   = i if entries[i]['len_nt'] >= entries[j]['len_nt'] else j
            smaller  = j if larger == i else i
            aa_larger  = nt_to_aa(entries[larger]['len_nt'])
            aa_smaller = nt_to_aa(entries[smaller]['len_nt'])

            if aa1 > AA_CUTOFF and aa2 > AA_CUTOFF:
                # Rule 3: both large — prefer bakta's call. Compare mapped tool names,
                # not raw sources: bakta emits 'Pyrodigal' for ab-initio CDS but 'Bakta'
                # for its own sORF pass, and a raw src == 'Pyrodigal' test loses the latter.
                # These usually differ in their stop codon, in which case they are separate
                # ORFs and the loser's tool tag is not carried across.
                if tool1 == 'bakta':
                    _absorb(entries, i, j, same_gene)
                elif tool2 == 'bakta':
                    _absorb(entries, j, i, same_gene)
                else:
                    # Neither is bakta — keep the larger call.
                    _absorb(entries, larger, smaller, same_gene)
                
            elif aa_larger >= AA_CUTOFF and aa_smaller < AA_CUTOFF and aa_diff > 0:
                # Rule 4: one large, one small, sufficiently different — distinct ORFs.
                # Keep both, and keep their tool sets separate: unioning here would
                # overcount support by crediting one gene's callers to the other.
                pass

            else:
                # Rule 5: keep the longer. Where the two share a stop codon this is one
                # gene with two start calls and the shorter one's tool carries across —
                # the dominant case for smORFs. Where they do not, only the coordinates
                # are dropped, not the attribution.
                if tool1 != 'orfipy' and tool2 != 'orfipy':
                    _absorb(entries, larger, smaller, same_gene)

    return [(e['cols'], e['tools']) for e in entries if e['keep']]

def merge_gff_files(input_files, output_file):
    
    fileid=input_files[0].split('/')[-1].replace('.gff3','').replace('.gff','').replace('.assembly','')

    all_cds    = {}
    non_cds = []
    
    region_lines = {}
    seq_region_lines = {}
    active_contigs = set()
    
    if "-" in fileid:
        fileid = fileid.split('-')[1]

    if "_collective_cds" in fileid:
        fileid = fileid.replace("_collective_cds","")
    if "_intergenic_cds" in fileid:
        fileid = fileid.replace("_intergenic_cds", "")

    for filepath in input_files:
        with open(filepath) as f:
            for line in f:
                line = line.rstrip('\n')
                
                # Collect header/comment lines from first file only
                if line.startswith('#'):
                    if line.startswith('##sequence-region'):
                        # ##sequence-region contig_1 1 1000
                        parts = line.split()
                        if len(parts) >= 2:
                            seqname = normalize_seqname(parts[1])
                            parts[1] = fileid+"-"+seqname
                            if seqname not in seq_region_lines:
                                seq_region_lines[seqname] = ' '.join(parts)  # store with normalized contig name
                    continue
                
                if not line.strip():
                    continue

                cols = line.split('\t')
                if len(cols) < 9:
                    continue

                feature_type = cols[2]

                # Dedup key: seqname, start, end, strand
                seqname = normalize_seqname(cols[0])
                
                if feature_type == 'region':
                    cols[0] = seqname
                    cols[8] = update_id_with_contig(cols[8], fileid+"-"+seqname)
                    sort_key = (natural_sort_key(seqname), int(cols[3]), int(cols[4]))
                    region_lines[seqname] = (sort_key, '\t'.join(cols))
                    continue
                
                if feature_type != 'CDS':
                   # Keep non-CDS features as-is (no dedup)
                    cols[0] = seqname
                    cols[8] = update_id_with_contig(cols[8], fileid+"-"+seqname)
                    sort_key = (natural_sort_key(seqname), int(cols[3]), int(cols[4]))
                    non_cds.append((sort_key, '\t'.join(cols)))
                    active_contigs.add(seqname)
                    continue
                    
                cols[0] = seqname
                strand  = cols[6]
                len_nt  = get_length(cols)
                key     = (seqname, strand) # group by contig + strand
                
                if key not in all_cds:
                    all_cds[key] = []
                all_cds[key].append((cols, len_nt))
    

    # Resolve overlaps per contig+strand group
    records = []
    resolved = []

    for (norm_seq, strand), cds_list in all_cds.items():
        # Sort by start position before resolving
        cds_list.sort(key=lambda x: int(x[0][3]))
        for cols, tools in resolve_overlaps(cds_list):
            cols[8] = update_id_with_contig(cols[8], fileid+"-"+norm_seq)
            active_contigs.add(norm_seq)
            sort_key = (natural_sort_key(norm_seq), int(cols[3]), int(cols[4]))
            resolved.append((sort_key, norm_seq, cols, sorted(tools)))

    # Number ORFs in output order so the GFF and the caller's table agree
    resolved.sort(key=lambda x: x[0])
    kept_cds = []
    for n, (sort_key, norm_seq, cols, tool_list) in enumerate(resolved, 1):
        orf_id = f"{fileid}_ORF{n:05d}"
        cols[8] = cols[8].rstrip(';') + \
            f";orf_id={orf_id};n_tools={len(tool_list)};tools={','.join(tool_list)}"
        records.append((sort_key, '\t'.join(cols)))
        kept_cds.append((orf_id, norm_seq, int(cols[3]), int(cols[4]), cols[6], tool_list))

    records.extend(non_cds)
    
    for norm_seq, (sort_key, line) in region_lines.items():
        if norm_seq in active_contigs:
            if norm_seq in seq_region_lines:
                records.append(((natural_sort_key(norm_seq), -1, 0), seq_region_lines[norm_seq]))
            # Use start=0 so region sorts before all features on that contig
            region_sort_key = (natural_sort_key(norm_seq), 0, 0)
            records.append((region_sort_key, line))

    digit_contigs = {c for c in active_contigs if c.isdigit()}
    nondigit_contigs = active_contigs - digit_contigs
    if digit_contigs and nondigit_contigs:
        raise ValueError(
            "Mixed numeric and non-numeric contig names after normalization — "
            "cannot sort consistently. Numeric examples: "
            f"{sorted(digit_contigs)[:5]}; non-numeric examples: {sorted(nondigit_contigs)[:5]}"
        )

    # Sort all together by contig, then start, then end
    records.sort(key=lambda x: x[0])
    with open(output_file, 'w') as out:
        out.write('##gff-version 3\n')
        for _, line in records:
            if line.startswith("##"):
                out.write(line + '\n')
            else:
                out.write(fileid+"-"+line + '\n')

    return kept_cds


if __name__ == '__main__':
    # Usage: python merge_gff.py file1.gff file2.gff file3.gff file4.gff merged.gff
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_file1> [<input_file2> ...] <output_file>")
        sys.exit(1)

    input_files = sys.argv[1:-1]
    output_file = sys.argv[-1]

    kept = merge_gff_files(input_files, output_file)
    print(f"Done. Written to: {output_file}")
    print(f"Total CDS records after resolving overlaps: {len(kept)}")
