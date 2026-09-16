"""Donor-level sequence-footprint analysis of Nanaware et al. 2025 Table S5.

Input: authors' background-filtered MHC-I peptide list, not raw hormone fragments.
Counterfactual: replace native INS1-24 with a donor signal peptide. Sequence
retention is a deterministic annotation, not predicted HLA display after editing.
"""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import warnings
import openpyxl

warnings.filterwarnings('ignore', message='Unknown extension is not supported')
BASE = Path(__file__).resolve().parent
S4 = BASE / 'NIHMS2107237-supplement-4.xlsx'
S5 = BASE / 'NIHMS2107237-supplement-5.xlsx'
INS = 'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN'
CHGA = 'MRSAAVLALLLCAGQVTALPVNSPMNKGDTEVMKCIVEVISDTLSKPSPMPVSQECFETLRGDERILSILRHQNLLKELQDLALQGAKERAHQQKKHSGFEDELSEVLENQSSQAELKEAVEEPSSKDVMEKREDSKEAEKSGEATDGARPQALPEPMQESKAEGNNQAPGEEEEEEEEATNTHPPASLPSQKYPGPQAEGDSEGLSQGLVDREKGLSAEPGWQAKREEEEEEEEEAEAGEEAVPEEEGPTVVLNPHPSLGYKEIRKGESRSEALAVDGAGKPGAEEAQDPEGKGEQEHSQQKEEEEEMAVVPQGLFRGGKSGELEQEEERLSKEWEDSKRWSKMDQLAKELTAEKRLEGQEEEEDNRDSSMKLSFRARAYGFRGPGPQLRRGWRPSSREDSLEAGLPLQVRGYPEEKKEEEGSANRRPEDQELESLSAIEAELEKVAHQLQALRRG'
REFERENCE = {'INS': INS, 'CHGA': CHGA}
# Amino acid strings verified against the protein-reference rows in Table S4.
COUNTERFACTUAL = CHGA[:18] + INS[24:]


def dump_csv(name, rows):
    with (BASE / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def shared_six(a, b):
    return sorted(set(a[i:i + 6] for i in range(len(a) - 5)) &
                  set(b[i:i + 6] for i in range(len(b) - 5)))


hormones = []
w = openpyxl.load_workbook(S4, read_only=True, data_only=True)
for sheet, gene in [('Insulin', 'INS'), ('Chromogranin', 'CHGA')]:
    rows = list(w[sheet].values)
    assert rows[1][2] == REFERENCE[gene]
    for excel_row, r in enumerate(rows[2:], start=3):
        if not r[2]:
            continue
        assert r[1] == gene
        record = dict(zip(['protein', 'gene', 'peptide', 'length', 'blank', 'class_i', 'dr', 'dp', 'dq', 'donor', 'start', 'end'], r[:12]))
        record['excel_row'] = excel_row
        assert REFERENCE[gene][record['start'] - 1:record['end']] == record['peptide']
        hormones.append(record)

w = openpyxl.load_workbook(S5, read_only=True, data_only=True)
all_rows = []
for excel_row, r in enumerate(list(w['Class I'].values)[1:], start=2):
    assert 8 <= r[3] <= 13 and r[3] == len(r[2]) and r[4] > 0
    all_rows.append({'excel_row': excel_row, 'protein': r[0], 'gene': r[1], 'peptide': r[2], 'length': r[3], 'ion_intensity': r[4], 'donor': r[5], 'predicted_hla_as_published': r[6]})
assert len(all_rows) == 3766
assert len({(r['donor'], r['peptide']) for r in all_rows}) == len(all_rows)

selected = []
for r in all_rows:
    if r['gene'] not in REFERENCE:
        continue
    r = dict(r)
    reference = REFERENCE[r['gene']]
    assert reference.count(r['peptide']) == 1
    r['start'] = reference.index(r['peptide']) + 1
    r['end'] = r['start'] + r['length'] - 1
    if r['gene'] == 'INS':
        if r['end'] <= 24:
            r['footprint'] = 'inside_native_INS_signal'
        elif r['start'] <= 24:
            r['footprint'] = 'crosses_native_INS_signal_boundary'
        else:
            r['footprint'] = 'downstream_INS_retained_sequence'
        r['sequence_in_canonical_CHGA1_18_INS25_110'] = r['peptide'] in COUNTERFACTUAL
    else:
        r['footprint'] = 'inside_native_CHGA_signal' if r['end'] <= 18 else 'downstream_CHGA'
        r['sequence_in_canonical_CHGA1_18_INS25_110'] = r['peptide'] in COUNTERFACTUAL
    matched = [h for h in hormones if (h['gene'], h['peptide'], h['donor']) == (r['gene'], r['peptide'], r['donor'])]
    assert len(matched) == 1
    h = matched[0]
    assert h['class_i'] == r['ion_intensity']
    r['S4_exact_blank_intensity'] = h['blank']
    same_donor_blanks = [h for h in hormones if h['gene'] == r['gene'] and h['donor'] == r['donor'] and h['blank'] is not None and h['blank'] > 0]
    overlaps = [h['peptide'] for h in same_donor_blanks if shared_six(h['peptide'], r['peptide'])]
    r['same_donor_blank_overlap_at_least_6aa'] = ';'.join(overlaps)
    selected.append(r)
assert len(selected) == 9
ins = [r for r in selected if r['gene'] == 'INS']
assert len(ins) == 6 and len({r['peptide'] for r in ins}) == 4
assert all((r['footprint'] == 'downstream_INS_retained_sequence') == r['sequence_in_canonical_CHGA1_18_INS25_110'] for r in ins)
dump_csv('annotated_INS_CHGA_mhci.csv', selected)

by_donor = []
for donor in sorted({r['donor'] for r in all_rows}):
    ir = [r for r in ins if r['donor'] == donor]
    cr = [r for r in selected if r['donor'] == donor and r['gene'] == 'CHGA']
    counts = Counter(r['footprint'] for r in ir)
    altered = counts['inside_native_INS_signal'] + counts['crosses_native_INS_signal_boundary']
    retained = counts['downstream_INS_retained_sequence']
    by_donor.append({'donor': donor, 'all_curated_mhci_species': sum(r['donor'] == donor for r in all_rows), 'INS_species': len(ir), 'INS_inside_signal': counts['inside_native_INS_signal'], 'INS_boundary_crossing': counts['crosses_native_INS_signal_boundary'], 'INS_retained_downstream': retained, 'INS_changed_fraction_of_detected_species': altered / len(ir) if ir else '', 'INS_retained_fraction_of_detected_species': retained / len(ir) if ir else '', 'retained_INS_sequences': ';'.join(sorted(r['peptide'] for r in ir if r['footprint'] == 'downstream_INS_retained_sequence')), 'CHGA_signal_species': sum(r['footprint'] == 'inside_native_CHGA_signal' for r in cr)})
dump_csv('donor_footprint.csv', by_donor)


def summarize(label, rows):
    peptides = {r['peptide'] for r in rows}
    altered = {r['peptide'] for r in rows if r['footprint'] != 'downstream_INS_retained_sequence'}
    kept = {r['peptide'] for r in rows if r['footprint'] == 'downstream_INS_retained_sequence'}
    d = {r['donor'] for r in rows}
    retained_d = {r['donor'] for r in rows if r['footprint'] == 'downstream_INS_retained_sequence'}
    return {'analysis': label, 'INS_donor_peptide_observations': len(rows), 'INS_distinct_species': len(peptides), 'changed_distinct_species': len(altered), 'retained_distinct_species': len(kept), 'INS_positive_donors': len(d), 'donors_with_retained_INS_species': len(retained_d)}

sensitivity = [summarize('authors_curated_S5', ins),
    summarize('drop_exact_blank_positive_keep_unknown', [r for r in ins if not r['S4_exact_blank_intensity']]),
    summarize('require_explicit_zero_exact_blank', [r for r in ins if r['S4_exact_blank_intensity'] == 0]),
    summarize('drop_same_donor_6aa_blank_overlap_keep_unknown', [r for r in ins if not r['same_donor_blank_overlap_at_least_6aa']]),
    summarize('exclude_lowest_depth_donor_HP18101_16_total_peptides', [r for r in ins if r['donor'] != 'HP18101'])]
dump_csv('background_sensitivity.csv', sensitivity)
loo = [summarize('exclude_' + d, [r for r in ins if r['donor'] != d]) for d in sorted({r['donor'] for r in all_rows})]
dump_csv('leave_one_donor_out.csv', loo)

# Quantify why the unfiltered hormone-elution table cannot be used as antigen list.
raw = [h for h in hormones if h['gene'] == 'INS' and h['class_i'] and h['class_i'] > 0]
raw_stats = {'S4_INS_classI_positive_donor_peptide_rows': len(raw), 'S4_INS_classI_positive_distinct_sequences': len({h['peptide'] for h in raw}), 'S4_INS_length8_13_rows': sum(8 <= h['length'] <= 13 for h in raw), 'S5_curated_INS_rows': len(ins), 'S5_curated_INS_distinct_sequences': len({r['peptide'] for r in ins})}

summary = {'total_classI_rows': len(all_rows), 'total_distinct_peptide_strings_across_donors': len({r['peptide'] for r in all_rows}), 'total_donors': len(by_donor), 'primary': sensitivity[0], 'raw_background_audit': raw_stats, 'native_INS_signal_range': [1, 24], 'canonical_CHGA_signal_range_for_sequence_only_counterfactual': [1, 18], 'counterfactual_sequence_is_not_verified_paper_construct': True, 'intensity_policy': 'Absolute recorded ion intensities archived for traceability. No pooling across donors/instruments, no comparison of different peptide intensities as molecular abundance, and no intensity-based protection estimate.', 'provenance': [{'file': f.name, 'sha256': hashlib.sha256(f.read_bytes()).hexdigest()} for f in [S4, S5]]}
(BASE / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
print('\nDONORS', json.dumps(by_donor, indent=2))
print('\nSENSITIVITY', json.dumps(sensitivity, indent=2))
