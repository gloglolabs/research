"""Render manuscript figures from the adjacent source-data CSV files.

Requires matplotlib and numpy. Reads only adjacent source-data tables.
"""
from pathlib import Path
import csv, hashlib, json
HERE = Path(__file__).resolve().parent
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.ticker import NullFormatter

plt.rcParams.update({
    'font.family':'DejaVu Sans','font.size':7.5,'axes.titlesize':8.5,
    'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,
    'axes.linewidth':.65,'lines.linewidth':1.,'pdf.fonttype':42,'ps.fonttype':42,
    'svg.fonttype':'none','svg.hashsalt':'chga-ins-manuscript-figures-v1','savefig.facecolor':'white','axes.spines.top':False,
    'axes.spines.right':False,'figure.facecolor':'white',
})
BLUE='#24689B'; ORANGE='#C86D2B'; PURPLE='#78659E'; RETAIN='#25768B'
SIGNAL='#C8793C'; MUT='#A63637'; GRAY='#D8DADD'; DARK='#252A30'
outputs=[]
def read(name):return list(csv.DictReader((HERE/name).open()))
def save(fig,name):
    for ext in ['pdf','svg','png']:
        path=HERE/f'{name}.{ext}'
        metadata={'CreationDate':None} if ext=='pdf' else {'Date':None} if ext=='svg' else {}
        fig.savefig(path,dpi=300,metadata=metadata)
        outputs.append(dict(file=path.name,width_inches=float(fig.get_figwidth()),height_inches=float(fig.get_figheight()),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    plt.close(fig)
def panel(ax,letter,title,pad=8):
    ax.set_title(title,loc='left',pad=pad,fontweight='normal')
    ax.text(-.13,1.03+max(0,pad-8)/100,letter,transform=ax.transAxes,fontsize=10,fontweight='bold',ha='left',va='bottom')
def bare(ax):
    for s in ax.spines.values():s.set_visible(False)
    ax.tick_params(length=0)

# FIGURE 1: architecture and sequence alignment, not proportional assay yields.
seqs=read('figure1_sequences.csv');design=read('figure1_scan_design.csv')
fig=plt.figure(figsize=(7,6.5));gs=fig.add_gridspec(3,1,height_ratios=[1.8,1.7,1.6],hspace=.65,left=.21,right=.97,top=.94,bottom=.08)
ax=fig.add_subplot(gs[0]);panel(ax,'A','Precursor architecture and replacement boundary')
for i,r in enumerate(seqs):
    y=2-i;end=int(r['signal_end']);length=len(r['sequence']);leadercolor=SIGNAL if i==0 else PURPLE
    ax.add_patch(Rectangle((0,y-.19),end,.38,facecolor=leadercolor,edgecolor='none'))
    ax.add_patch(Rectangle((end,y-.19),length-end,.38,facecolor=RETAIN,edgecolor='none'))
    ax.text((end+length)/2,y,'INS25–110',color='white',ha='center',va='center',fontsize=7.2)
    ax.text(end,y+.27,f'{end}/{end+1}',ha='center',va='bottom',fontsize=7)
    ax.text(length+1.5,y,str(length),ha='left',va='center',fontsize=7)
    if i==2:
        for p in [9,17]:ax.plot(p-.5,y,'|',color='white',markersize=10,markeredgewidth=1.6)
ax.set_xlim(0,115);ax.set_ylim(-.5,2.65);ax.set_yticks([2,1,0],[r['construct'] for r in seqs]);ax.tick_params(axis='y',length=0,pad=7)
ax.set_xticks([1,18,24,50,75,100]);ax.set_xlabel('Precursor residue (construct coordinates)')
ax.spines['left'].set_visible(False)
ax.legend(handles=[Patch(color=SIGNAL,label='Native INS signal'),Patch(color=PURPLE,label='CHGA signal'),Patch(color=RETAIN,label='Retained proinsulin')],loc='upper center',bbox_to_anchor=(.5,-.22),ncol=3,frameon=False,fontsize=6.8,handlelength=1.2,columnspacing=1.4)
ax=fig.add_subplot(gs[1]);panel(ax,'B','Alignment anchored at the cleavage boundary')
for i,r in enumerate(seqs):
    y=2-i;end=int(r['signal_end']);leader=r['sequence'][:end];offset=24-end
    for j,char in enumerate(leader):
        x=offset+j;color=DARK
        if i==2 and j+1 in [9,17]:
            ax.add_patch(Rectangle((x-.46,y-.28),.92,.56,facecolor='#F5DEDA',edgecolor=MUT,lw=.7));color=MUT
        ax.text(x,y,char,ha='center',va='center',family='DejaVu Sans Mono',fontsize=7.6,color=color)
    for j,char in enumerate(r['sequence'][end:end+6]):ax.text(24+j,y,char,ha='center',va='center',family='DejaVu Sans Mono',fontsize=7.6,color=RETAIN)
ax.axvline(23.5,color='#5F646B',lw=.7,linestyle='--')
ax.set_yticks([2,1,0],[r['construct'] for r in seqs]);ax.set_xlim(-.8,30);ax.set_ylim(-.7,2.7)
ax.set_xticks([]);bare(ax)
ax.text(29.8,-.6,'INS25–30',ha='right',va='center',color=RETAIN,fontsize=7)
ax.text(14,-.6,'CHGA9',ha='center',va='center',fontsize=6.5,color=MUT)
ax.text(22,-.6,'CHGA17',ha='center',va='center',fontsize=6.5,color=MUT)
ax=fig.add_subplot(gs[2]);panel(ax,'C','Search stages and common evaluation grid');ax.axis('off')
columns=['Search stage','Position(s)','Substitutions','Controls','Total']
widths=[.38,.17,.18,.12,.15];xs=np.cumsum([0]+widths)
for j,h in enumerate(columns):ax.text(xs[j]+.008,.89,h,fontweight='bold',fontsize=7,transform=ax.transAxes)
for i,r in enumerate(design):
    vals=[r['stage'],r['positions'],r['substitutions'],r['controls'],r['scored_sequences']]
    y=.66-i*.22
    for j,v in enumerate(vals):ax.text(xs[j]+.008,y,v,fontsize=7.2,transform=ax.transAxes)
    ax.plot([0,1],[y-.065,y-.065],transform=ax.transAxes,color='#DDDFE2',lw=.6)
ax.text(0,-.12,'7 HLA alleles × 2 presentation models × 3 fixed cutoffs = 42 criteria',transform=ax.transAxes,fontsize=7.2)
ax.text(0,-.32,'8–11mers; 504 leader/junction peptide–HLA pairs per 104-aa construct',transform=ax.transAxes,fontsize=7.0,color='#555B62')
save(fig,'figure1_architecture')

# FIGURE 2: all explicit double substitutions, never biological replicates.
counts=read('figure2_counts.csv');rob=read('figure2_robustness.csv');totals=read('figure2_totals.csv')
variants=[f'L9{x}_T17{y}' for x in ['C','G'] for y in ['C','G','D','P']]
alleles=['HLA-A*02:01','HLA-A*03:01','HLA-A*24:02','HLA-B*40:01','HLA-B*49:01','HLA-C*03:04','HLA-C*07:01']
cutoffs=[.5,1.,2.]
fig=plt.figure(figsize=(7,7.2));gs=fig.add_gridspec(3,2,width_ratios=[1.03,1],height_ratios=[2.1,2.1,2.7],hspace=.71,wspace=.45,left=.18,right=.925,top=.925,bottom=.105)
maxdelta=max(abs(int(r['delta'])) for r in counts)
negative=[plt.cm.Blues(x) for x in np.linspace(.88,.2,maxdelta)]
positive=[plt.cm.Oranges(x) for x in np.linspace(.3,.88,maxdelta)]
cmap=ListedColormap(negative+['#F6F6F6']+positive);norm=BoundaryNorm(np.arange(-maxdelta-.5,maxdelta+1.5),cmap.N)
for p,(model,letter,title) in enumerate([('NetMHCpan EL','A','NetMHCpan EL: change in candidate count'),('MHCflurry presentation','B','MHCflurry presentation, full flanks: change in count')]):
    ax=fig.add_subplot(gs[p,:]);panel(ax,letter,title,pad=25)
    arr=np.array([[int(next(r['delta'] for r in counts if r['variant']==v and r['model']==model and r['context']=='with_flanks' and r['allele']==a and float(r['cutoff_percent'])==t)) for a in alleles for t in cutoffs] for v in variants])
    im=ax.imshow(arr,cmap=cmap,norm=norm,aspect='auto',interpolation='none')
    for i in range(8):
        for j in range(21):
            z=arr[i,j]
            if z:ax.text(j,i,f'{z:+d}',ha='center',va='center',fontsize=6.3,color='white' if abs(z)>=maxdelta*.62 else DARK)
    ax.set_yticks(range(8),[v.replace('_',' + ') for v in variants]);ax.get_yticklabels()[0].set_fontweight('bold')
    ax.set_xticks(range(21),['0.5','1','2']*7);ax.tick_params(length=0,pad=3)
    for j,a in enumerate(alleles):
        ax.text((j*3+1.5)/21,1.045,a.replace('HLA-',''),ha='center',va='bottom',transform=ax.transAxes,fontsize=7)
        if j:ax.axvline(j*3-.5,color='white',lw=1.3)
    ax.add_patch(Rectangle((-.5,-.5),21,1,fill=False,edgecolor=DARK,lw=.8))
    ax.set_xlabel('Rank cutoff (%)',labelpad=3)
    for s in ax.spines.values():s.set_visible(False)
cbax=fig.add_axes([.943,.535,.012,.34]);cb=fig.colorbar(im,cax=cbax);cb.set_ticks([-maxdelta,0,maxdelta]);cb.ax.tick_params(labelsize=6,length=2)
cb.ax.set_title('Δ',fontsize=7,pad=4)
ax=fig.add_subplot(gs[2,0]);panel(ax,'C','Criterion robustness across contexts',pad=8)
for i,v in enumerate(variants):
    for offset,mode,hatch in [(-.16,'with_flanks',''),(.16,'without_flanks','////')]:
        r=next(r for r in rob if r['variant']==v and r['context']==mode);left=0
        for field,col in [('improved',BLUE),('unchanged',GRAY),('worsened',ORANGE)]:
            value=int(r[field]);ax.barh(i+offset,value,left=left,height=.27,color=col,edgecolor='white',linewidth=.25,hatch=hatch);left+=value
        if int(r['worsened']):ax.text(42.8,i+offset,r['worsened'],va='center',fontsize=6.2,color=ORANGE)
ax.set_ylim(7.6,-.6);ax.set_yticks(range(8),[v.replace('_',' / ') for v in variants]);ax.set_xlim(0,47);ax.set_xticks([0,14,28,42]);ax.set_xlabel('Number of criteria (42 per context)')
ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
ax.legend(handles=[Patch(color=BLUE,label='Improved'),Patch(color=GRAY,label='Unchanged'),Patch(color=ORANGE,label='Worsened')],loc='upper left',bbox_to_anchor=(-.30,-.19),ncol=3,frameon=False,fontsize=6.5,handlelength=1.1,columnspacing=.7)

sub=gs[2,1].subgridspec(3,1,hspace=.75)
for j,(model,mode,title) in enumerate([('NetMHCpan EL','with_flanks','NetMHCpan EL'),('MHCflurry presentation','with_flanks','MHCflurry, full flanks'),('MHCflurry presentation','without_flanks','MHCflurry, no flanks')]):
    ax=fig.add_subplot(sub[j]);
    if j==0:panel(ax,'D','Baseline versus selected variant',pad=20)
    for v,col,marker in [('baseline','#6E747B','o'),('L9C_T17C',BLUE,'s')]:
        values=[int(next(r['count'] for r in totals if r['variant']==v and r['model']==model and r['context']==mode and float(r['cutoff_percent'])==t)) for t in cutoffs]
        ax.plot(range(3),values,color=col,marker=marker,markersize=3.5,label='Baseline' if v=='baseline' else 'L9C + T17C')
        for x,y in enumerate(values):ax.annotate(str(y),(x,y),xytext=(0,4 if v=='baseline' else -5),textcoords='offset points',ha='center',fontsize=6.2,color=col)
    ax.text(.02,1.04,title,transform=ax.transAxes,fontsize=6.8,va='bottom')
    ax.set_xlim(-.25,2.25);ax.set_ylim(-10,38);ax.set_xticks(range(3),['0.5','1','2']);ax.set_yticks([0,15,30]);ax.tick_params(labelsize=6.3,pad=2)
    ax.set_ylabel('Count',fontsize=6.8)
    if j==2:ax.set_xlabel('Rank cutoff (%)',fontsize=7)
    if j==0:fig.legend(handles=[Line2D([],[],color='#6E747B',marker='o',markersize=3.5,label='Baseline'),Line2D([],[],color=BLUE,marker='s',markersize=3.5,label='L9C + T17C')],loc='lower center',bbox_to_anchor=(.77,.032),ncol=2,frameon=False,fontsize=6.3,handlelength=1.6)
save(fig,'figure2_design_comparison')

# FIGURE 3: exact presence and sequence projection, no intensity weighting.
data=read('figure3_donor_projection.csv');peptides=['ALWGPDPAAA','AAAFVNQHL','GSHLVEALY','HLVEALYLV'];donors=['HP20289','HP18101','HP19026','R360','R361','R369']
fig=plt.figure(figsize=(7,5.3));gs=fig.add_gridspec(2,1,height_ratios=[1.5,1.8],hspace=.6,left=.245,right=.96,top=.92,bottom=.19)
ax=fig.add_subplot(gs[0]);panel(ax,'A','Projection of curated INS ligands onto the native precursor')
ax.axvspan(10,24.5,color='#FAECDD',zorder=0);ax.axvspan(24.5,46,color='#E7F0F3',zorder=0);ax.axvspan(31.5,42.5,color='#BDD9E0',alpha=.55,zorder=0)
for i,pep in enumerate(peptides):
    r=next(r for r in data if r['peptide']==pep);start,end=int(r['start']),int(r['end']);color=RETAIN if r['footprint']=='retained' else SIGNAL
    ax.add_patch(Rectangle((start-.45,i-.20),end-start+.9,.4,facecolor=color,edgecolor='none'))
    ax.text(start-.8,i,str(start),ha='right',va='center',fontsize=6.5);ax.text(end+.8,i,str(end),ha='left',va='center',fontsize=6.5)
ax.axvline(24.5,color='#7B726A',ls='--',lw=.8);ax.text(24.5,-.55,'24/25',ha='center',fontsize=6.8)
ax.text(37,-.55,'INS32–42 / B8–18',ha='center',fontsize=7,color=RETAIN)
ax.set_yticks(range(4),peptides);ax.set_xlim(10,46);ax.set_ylim(3.5,-.85);ax.set_xticks([10,15,20,25,30,35,40,45]);ax.set_xlabel('Native INS residue')
ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False)
ax=fig.add_subplot(gs[1]);panel(ax,'B','Observed donor–peptide mapping after published curation')
for i,pep in enumerate(peptides):
    for j,donor in enumerate(donors):
        r=next(r for r in data if r['peptide']==pep and r['donor']==donor)
        color=RETAIN if r['footprint']=='retained' else SIGNAL
        ax.add_patch(Rectangle((j-.43,i-.35),.86,.7,facecolor=color if int(r['detected']) else '#F0F1F3',edgecolor='white',lw=.4))
        if int(r['detected']):ax.plot(j,i,'o',color='white',ms=3.5)
ax.set_xlim(-.5,5.5);ax.set_ylim(3.6,-.6);ax.set_yticks(range(4),peptides)
depths=[int(next(r['total_curated_classI_species'] for r in data if r['donor']==d)) for d in donors]
ax.set_xticks(range(6),[d+'\n'+f'{n:,}' for d,n in zip(donors,depths)]);ax.tick_params(length=0,pad=6);bare(ax)
ax.set_xlabel('Donor preparation\nSecond line: total curated class-I peptide species',labelpad=8)
fig.legend(handles=[Patch(color=SIGNAL,label='Changed INS sequence'),Patch(color=RETAIN,label='Retained INS sequence'),Patch(color='#F0F1F3',label='Not detected in curated table')],loc='lower center',bbox_to_anchor=(.57,.015),ncol=3,frameon=False,fontsize=6.8,handlelength=1.2,columnspacing=1.3)
save(fig,'figure3_donor_projection')

# FIGURE 4: model disagreement and concrete residual/affinity tradeoffs.
junction=read('figure4_junction_affinity.csv');residual=read('figure4_residual_candidates.csv');trade=read('figure4_A02_affinity_tradeoff.csv')
fig=plt.figure(figsize=(7,6.0));gs=fig.add_gridspec(2,2,height_ratios=[1.55,1.1],width_ratios=[1.15,1],hspace=.72,wspace=.82,left=.23,right=.975,top=.92,bottom=.12)
ax=fig.add_subplot(gs[0,0]);panel(ax,'A','C*03:04 junction affinity')
for i,r in enumerate(junction):
    for off,field,col,mark in [(-.10,'net_affinity_nm',BLUE,'o'),(.10,'mhc_affinity_nm',ORANGE,'s')]:ax.plot(float(r[field]),i+off,mark,color=col,ms=4)
ax.set_xscale('log');ax.xaxis.set_minor_formatter(NullFormatter());ax.set_xlim(10,1e4);ax.set_xticks([10,100,1000,10000],['10','100','1,000','10,000']);ax.set_yticks(range(4),[r['peptide']+'\n'+r['label'] for r in junction]);ax.set_ylim(3.5,-.6);ax.set_xlabel('Predicted affinity (nM)');ax.tick_params(axis='y',length=0,pad=6)
ax.grid(axis='x',color='#E3E5E7',lw=.5);ax.spines['left'].set_visible(False)
ax=fig.add_subplot(gs[0,1]);panel(ax,'B','Residual candidates: rank')
for i,r in enumerate(residual):
    ax.plot(float(r['net_el_rank_percent']),i-.10,'o',color=BLUE,ms=4)
    ax.plot(float(r['mhc_presentation_rank_percent']),i+.10,'s',color=ORANGE,ms=4)
ax.axvline(.5,color='#92979D',ls=':',lw=.8);ax.set_xscale('log');ax.xaxis.set_minor_formatter(NullFormatter());ax.set_xlim(.1,2);ax.set_xticks([.1,.5,1,2],['0.1','0.5','1','2']);ax.set_yticks(range(3),[r['peptide'] for r in residual]);ax.set_ylim(2.5,-.5);ax.set_xlabel('Model-specific rank (%)');ax.tick_params(axis='y',length=0,pad=5);ax.spines['left'].set_visible(False)
ax.grid(axis='x',color='#E3E5E7',lw=.5)
ax=fig.add_subplot(gs[1,:]);panel(ax,'C','A*02:01 leader affinity tradeoff')
for i,r in enumerate(trade):
    for off,field,col,mark in [(-.10,'net_affinity_nm',BLUE,'o'),(.10,'mhc_affinity_nm',ORANGE,'s')]:
        value=float(r[field]);ax.plot(value,i+off,mark,color=col,ms=4)
        ax.annotate(f'{value:,.0f}',(value,i+off),xytext=(4,0),textcoords='offset points',va='center',fontsize=6.3,color=col)
labels=['Canonical\nVLALLLCAGQV','L9C + T17C\nVLACLLCAGQV','L9G + T17C\nVLAGLLCAGQV']
ax.set_yticks(range(3),labels);ax.set_ylim(2.5,-.5);ax.set_xscale('log');ax.xaxis.set_minor_formatter(NullFormatter());ax.set_xlim(20,25000);ax.set_xticks([100,1000,10000],['100','1,000','10,000']);ax.set_xlabel('Predicted affinity (nM)');ax.tick_params(axis='y',length=0,pad=6);ax.spines['left'].set_visible(False);ax.grid(axis='x',color='#E3E5E7',lw=.5)
fig.legend(handles=[Line2D([],[],color=BLUE,marker='o',linestyle='none',markersize=4,label='NetMHCpan 4.1'),Line2D([],[],color=ORANGE,marker='s',linestyle='none',markersize=4,label='MHCflurry 2.2.1')],loc='lower center',bbox_to_anchor=(.57,.015),ncol=2,frameon=False,fontsize=7)
save(fig,'figure4_model_comparison')

(HERE/'figure_manifest.json').write_text(json.dumps({'matplotlib_version':matplotlib.__version__,'numpy_version':np.__version__,'outputs':outputs,'rendering':'vector PDF/SVG and300dpi PNG; editable SVG text; embedded PDF TrueType fonts'},indent=2))
print('\n'.join(x['file'] for x in outputs))
