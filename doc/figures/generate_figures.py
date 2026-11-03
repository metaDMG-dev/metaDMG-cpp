"""Generate schematic README figures; no empirical data are used.
Requires matplotlib and numpy. Run with python generate_figures.py.
"""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch
OUT=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
ink='#213c46'; teal='#087e8b'; orange='#c97432'; pale='#eaf4f4'
def save(fig,name):
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight',facecolor='white')
    svg = OUT/(name+'.svg')
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    fig.savefig(OUT/(name+'.png'),dpi=160,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def box(ax,x,y,w,h,title,sub,color=pale):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.015,rounding_size=0.02',facecolor=color,edgecolor='#c7d8dc'))
    ax.text(x+w/2,y+h*.69,title,ha='center',va='center',weight='bold',color=ink,fontsize=12)
    ax.text(x+w/2,y+h*.3,sub,ha='center',va='center',color=ink,fontsize=10)
def arrow(ax,a,b):
    ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':teal,'lw':2})
fig,ax=plt.subplots(figsize=(12,4.4));ax.set(xlim=(-.025,1.025),ylim=(0,1));ax.axis('off')
ax.text(0,1,'From mapped reads to interpretable damage estimates',weight='bold',fontsize=17,color=ink)
box(ax,.01,.43,.22,.32,'1  Prepare alignments','MD/NM tags retained\nSort by read name')
box(ax,.31,.64,.28,.22,'2a  getdamage','Whole sample / reference / taxid')
box(ax,.31,.16,.28,.26,'2b  lca','Read classification\n+ matching taxonomy')
box(ax,.68,.43,.29,.32,'3  dfit','Fit terminal damage\nfrom .bdamage.gz counts')
arrow(ax,(.24,.65),(.30,.74));arrow(ax,(.24,.50),(.30,.30));arrow(ax,(.60,.74),(.67,.65));arrow(ax,(.60,.30),(.67,.50))
ax.text(.81,.19,'LCA workflow: aggregate',ha='center',color=ink,weight='bold',fontsize=11)
ax.text(.81,.09,'Roll up statistics + attach fits',ha='center',color=ink,fontsize=10)
arrow(ax,(.82,.42),(.82,.25))
ax.text(.01,-.02,'Counting → model fitting → interpretation. Mapping and read preparation happen upstream.',color=ink,fontsize=11)
save(fig,'workflow')
fig,(a,b)=plt.subplots(1,2,figsize=(12,4.3),gridspec_kw={'width_ratios':[1,1]})
for ax in (a,b):ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
a.set_title('LCA: preserve mapping ambiguity',loc='left',weight='bold',color=ink,pad=18)
box(a,.29,.70,.42,.20,'Genus G','Lowest shared ancestor')
box(a,.02,.27,.39,.19,'Species A','Reference A')
box(a,.59,.27,.39,.19,'Species B','Reference B')
for x in (.21,.79):a.plot([.5,x],[.70,.46],color=teal,lw=2)
a.text(.5,.10,'One read aligns to both species → genus LCA',ha='center',color=ink,fontsize=10)
a.text(.5,-.02,'With --lca_rank species, this group is not\naccepted for species-level damage/statistics.',ha='center',fontsize=10,color=ink)
b.set_title('Weighting: a separate counting decision',loc='left',weight='bold',color=ink,pad=18)
box(b,.08,.70,.84,.20,'One accepted read group','Two contributing alignments (simple example)')
box(b,.04,.28,.42,.25,'--weight_type 0','1 + 1 = 2\nFull alignment weight',color='#fff1e6')
box(b,.55,.28,.42,.25,'--weight_type 1','½ + ½ = 1\nFractional weight')
b.text(.5,.08,'Specify the weighting explicitly in your command.',ha='center',fontsize=10,color=ink)
b.text(.5,-.02,'These are damage-count contributions,\nnot a definition of every output count column.',ha='center',fontsize=10,color=ink)
fig.subplots_adjust(wspace=.22);save(fig,'lca-weighting')
x=np.arange(15);A=.20;q=.30;c=.02
fig,(a,b)=plt.subplots(1,2,figsize=(11.5,4.5))
for ax in (a,b):
    ax.set(xlabel='Distance from read end (x; terminal base = 0)',ylabel='Expected damage frequency',ylim=(0,.26),xlim=(-.3,14.3))
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1))
    ax.grid(axis='y',alpha=.17)
a.plot(x,A*(1-q)**x+c,color=teal,lw=2.7,label='D(x) = A(1 − q)ˣ + c')
a.axhline(c,color=orange,ls='--',label='Background c = 2%')
a.annotate('',xy=(.3,A+c),xytext=(.3,c),arrowprops={'arrowstyle':'<->','color':ink})
a.text(2.3,.14,'Amplitude A = 20%',color=ink,fontsize=10)
a.text(5,.23,'First-base frequency: A + c = 22%',fontsize=10,color=ink)
a.set_title('Amplitude and background',loc='left',weight='bold',color=ink)
a.legend(frameon=False,loc='upper right',bbox_to_anchor=(1,.84),fontsize=9)
for qi,col in ((.15,orange),(.30,teal),(.60,'#647286')):
    b.plot(x,A*(1-qi)**x+c,lw=2.3,label=f'q = {qi:.2f}',color=col)
b.set_title('Larger q → faster terminal decay',loc='left',weight='bold',color=ink)
b.legend(frameon=False)
fig.text(.5,-.02,'Illustrative model curves only • chosen parameters, no measured data • phi controls dispersion, not this mean curve',ha='center',fontsize=10,color=ink)
fig.tight_layout(w_pad=2);save(fig,'damage-model')
