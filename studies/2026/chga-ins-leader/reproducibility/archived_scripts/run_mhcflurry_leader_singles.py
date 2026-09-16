"""Independent MHCflurry2.2.1 cross-model check of the canonical boundary study.

Official 20200611 presentation weights, distribution release2.2.0.
Every HLA is a singleton sample to preserve allele-specific predictions.
CPU execution with two threads; full sequence scans8-11, no output filtering.
"""
import os
from pathlib import Path
HERE=Path(__file__).resolve().parent
os.environ['MHCFLURRY_DATA_DIR']=str(HERE/'mhcflurry_data')
import datetime,hashlib,json,platform,subprocess
import mhcflurry,torch,numpy,pandas
from mhcflurry import Class1PresentationPredictor
from mhcflurry.common import configure_pytorch
configure_pytorch(backend='cpu',num_threads=2)
OUT=HERE/'mhcflurry_leader_singles';OUT.mkdir(exist_ok=True)
base=json.loads((HERE/'scan_manifest.json').read_text())
leader_manifest=json.loads((HERE.parent/'leader_design/scan_manifest.json').read_text())
constructs=leader_manifest['constructs']
full_sequences={c['id']:c['full_sequence'] for c in constructs}
# Retained peptides start1-18, end<=28; five C-flank residues requireprefix33.
sequences={name:seq[:33] for name,seq in full_sequences.items()}
alleles={a:[a] for a in base['alleles']}
model_dir=HERE/'mhcflurry_data/2.2.0/models_class1_presentation/models'
predictor=Class1PresentationPredictor.load(str(model_dir))
model_hashes={str(f.relative_to(model_dir)):hashlib.sha256(f.read_bytes()).hexdigest() for f in model_dir.rglob('*') if f.is_file()}
metadata=dict(code_version=mhcflurry.__version__,torch_version=torch.__version__,python=platform.python_version(),
    numpy_version=numpy.__version__,pandas_version=pandas.__version__,platform=platform.platform(),
    model_distribution='2.2.0',model_download='https://github.com/openvax/mhcflurry/releases/download/pre-2.0/models_class1_presentation.20200611.tar.bz2',
    backend='cpu',threads=2,peptide_lengths=[8,9,10,11],alleles=alleles,
    supported_peptide_lengths=list(predictor.supported_peptide_lengths),
    all_requested_alleles_supported=all(a in predictor.supported_alleles for a in alleles),
    requested_allele_support={a:a in predictor.supported_alleles for a in alleles},
    model_file_sha256=model_hashes,started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    scope='Canonical CHGAleader134single-variant sensitivity only; exact experimental construct notclaimed; endogenous CHGA remains a separate peptide source')
(OUT/'manifest.json').write_text(json.dumps(metadata,indent=2))
(OUT/'input_sequences.json').write_text(json.dumps(sequences,indent=2))
(OUT/'full_sequences.json').write_text(json.dumps(full_sequences,indent=2))
(OUT/'supported_alleles.txt').write_text('\n'.join(predictor.supported_alleles)+'\n')
print(json.dumps({k:v for k,v in metadata.items() if k not in ['model_file_sha256']}),flush=True)
assert metadata['all_requested_alleles_supported']
expected=sum(sum(len(s)-n+1 for n in [8,9,10,11]) for s in sequences.values())*len(alleles)
for use_flanks in [True,False]:
    name='with_flanks' if use_flanks else 'without_flanks'
    path=OUT/f'{name}.csv'
    if path.exists():print(f'Cached {path}',flush=True);continue
    result=predictor.predict_sequences(sequences=sequences,alleles=alleles,result='all',peptide_lengths=[8,9,10,11],
        use_flanks=use_flanks,include_affinity_percentile=True,verbose=0,throw=True)
    assert len(result)==expected
    result=result[result.pos < 18].copy()
    assert len(result)==len(constructs)*len(alleles)*18*4
    for row in result.itertuples():
        assert sequences[row.sequence_name][row.pos:row.pos+len(row.peptide)]==row.peptide
        assert row.best_allele==row.sample_name
        if use_flanks:
            full=full_sequences[row.sequence_name]
            assert row.n_flank==full[max(0,row.pos-5):row.pos]
            assert row.c_flank==full[row.pos+len(row.peptide):row.pos+len(row.peptide)+5]
    assert result[['affinity','affinity_percentile','processing_score','presentation_score','presentation_percentile']].notna().all().all()
    result.to_csv(path,index=False)
    print(f'Wrote {len(result)} validated rows to {path}',flush=True)
metadata['evaluated_prefix_grid_rows_per_context']=expected
metadata['retained_grid_rows_per_context']=len(constructs)*len(alleles)*18*4
metadata['retained_region']='peptide start1-18 inclusive,8-11mers'
metadata['flank_validation']='five flanks exactly matchfull protein for every retained row'
metadata['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(OUT/'manifest.json').write_text(json.dumps(metadata,indent=2))
