#!/usr/bin/env python3
"""Optional exact/I-L-equivalent rescan of separately retrieved human proteomes."""
import argparse,gzip,hashlib,json
from pathlib import Path
R=Path(__file__).resolve().parent

def records(stream):
    head=None;parts=[]
    for line in stream:
        if line.startswith('>'):
            if head is not None:yield head,''.join(parts)
            head=line[1:].strip();parts=[]
        else:parts.append(line.strip())
    if head is not None:yield head,''.join(parts)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--directory',type=Path,required=True,help='Directory containing the two FASTA.gz files listed in retrieval_manifest.json');ap.add_argument('--output',type=Path,default=R/'results/reference_match_rescan.json');a=ap.parse_args()
    expected=json.loads((R/'provenance/novelty/proteome_exact_matches.json').read_text())
    result={p:{'exact':[],'IL_equivalent':[]} for p in expected['matches']};count=residues=0;provenance=[]
    for source in expected['source_files']:
        f=a.directory/source['file'];h=hashlib.sha256();md5=hashlib.md5()
        with f.open('rb') as stream:
            for b in iter(lambda:stream.read(1048576),b''):h.update(b);md5.update(b)
        assert h.hexdigest()==source['sha256'],f'Wrong source snapshot: {f.name}'
        assert md5.hexdigest()==source['published_md5']
        n=0
        with gzip.open(f,'rt') as stream:
            for header,seq in records(stream):
                n+=1;count+=1;residues+=len(seq);collapsed=seq.replace('I','L')
                for p in result:
                    for label,needle,hay in [('exact',p,seq),('IL_equivalent',p.replace('I','L'),collapsed)]:
                        i=hay.find(needle)
                        while i>=0:
                            result[p][label].append({'protein':header,'position_1based':i+1,'actual_sequence':seq[i:i+len(p)]});i=hay.find(needle,i+1)
        assert n==source['sequence_count'];provenance.append(source)
    assert count==expected['protein_sequence_records'] and residues==expected['total_residues']
    assert result==expected['matches'],'Match result differs from the manuscript snapshot'
    report={'status':'PASS','reference_release':expected['UniProt_release'],'protein_sequence_records':count,'total_residues':residues,'matches':result,'source_files':provenance}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','protein_sequence_records':count,'counts':{p:{k:len(v) for k,v in d.items()} for p,d in result.items()}},indent=2))

if __name__=='__main__':main()
