import os, pickle, warnings
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as FA
import mpl_toolkits.axisartist.grid_finder as GF
from matplotlib.projections import PolarAxes
from scipy.stats import pearsonr
warnings.filterwarnings('ignore')

# Analysis working tree = parent of the repo; this file lives at
# <repo>/extremes_validation/ (workspace copy: <workspace>/WRF-CanESM2/extremes_validation/).
REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.environ.get('WRF_CANESM2_BASE', os.path.dirname(REPO_DIR))
MODELS = ['D3','D15','D75','CanRCM4','CanESM2']
LABELS = {'D3':'CanESM2-WRF D03','D15':'CanESM2-WRF D02','D75':'CanESM2-WRF D01',
          'CanESM2':'CanESM2','CanRCM4':'CanRCM4'}
COLORS = {'D3':'#1f77b4','D15':'#ff7f0e','D75':'#2ca02c','CanESM2':'#d62728','CanRCM4':'#9467bd'}
MARKERS = {'D3':'s','D15':'^','D75':'o','CanESM2':'P','CanRCM4':'p'}
FONT=9; FONT_S=7.5
plt.rcParams.update({'font.size':FONT,'axes.labelsize':FONT_S,'xtick.labelsize':FONT_S,
    'ytick.labelsize':FONT_S,'axes.titlesize':FONT,'legend.fontsize':FONT_S,
    'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
    'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':600})
CORR_TICKS=(0.0,0.2,0.4,0.6,0.8,0.9,0.95,1.0)

class TaylorDiagram:
    def __init__(self, refstd, fig, rect, max_std, unit='unit'):
        self.refstd=refstd
        tr=PolarAxes.PolarTransform()
        tlocs=np.arccos(np.array(CORR_TICKS,dtype=float))
        gl1=GF.FixedLocator(tlocs); tf1=GF.DictFormatter(dict(zip(tlocs,map(str,CORR_TICKS))))
        self.smin=0.0; self.smax=max_std
        std_ticks=np.linspace(0,max_std,6)
        gl2=GF.FixedLocator(std_ticks)
        ghelper=FA.GridHelperCurveLinear(tr,extremes=(0,np.pi/2,0,max_std),
            grid_locator1=gl1,tick_formatter1=tf1,grid_locator2=gl2)
        ax=FA.FloatingSubplot(fig,rect,grid_helper=ghelper); fig.add_axes(ax)
        ax.axis['top'].set_axis_direction('bottom'); ax.axis['top'].toggle(ticklabels=False,label=False)
        ax.axis['left'].set_axis_direction('bottom'); ax.axis['left'].label.set_text('Standard deviation (%s)'%unit); ax.axis['left'].label.set_pad(2)
        ax.axis['right'].set_axis_direction('top'); ax.axis['right'].toggle(ticklabels=True)
        ax.axis['right'].major_ticklabels.set_axis_direction('left')
        ax.axis['bottom'].set_visible(False)
        self._ax=ax; self.ax=ax.get_aux_axes(tr)
        self.ax.plot([0.001],refstd,'k*',ls='',ms=13,zorder=1000,clip_on=False,label='observations')
        t=np.linspace(0,np.pi/2); r=np.zeros_like(t)+refstd
        self.ax.plot(t,r,'k--',lw=1.0,zorder=1)
    def add_sample(self,stddev,corr,**kw):
        self.ax.plot(np.arccos(np.clip(corr,-1,1)),stddev,**kw)
    def add_contours(self,levels,**kw):
        rs,ts=np.meshgrid(np.linspace(0,self.smax,120),np.linspace(0,np.pi/2,120))
        rms=np.sqrt(self.refstd**2+rs**2-2*self.refstd*rs*np.cos(ts))
        return self.ax.contour(ts,rs,rms,levels,**kw)
    def label_contours(self,ct,fmt='%.1f',corr=0.55):
        r_mid=0.5*self.smax
        for lvl in ct.levels:
            a=np.arccos(np.clip(corr,-1,1))
            if lvl<self.refstd: a=max(a,np.arcsin(lvl/self.refstd)+0.08)
            a=min(a,0.97*np.pi/2)
            b=-2.0*self.refstd*np.cos(a); disc=b*b-4.0*(self.refstd**2-lvl**2)
            if disc<0: continue
            roots=[(-b-np.sqrt(disc))/2.0,(-b+np.sqrt(disc))/2.0]
            roots=[r for r in roots if 1e-9<=r<=self.smax-1e-9]
            if not roots: continue
            rr=min(roots,key=lambda r:abs(r-r_mid))
            self.ax.text(a,rr,fmt%lvl,ha='center',va='center',fontsize=FONT_S,color='0.35',zorder=500,
                bbox=dict(boxstyle='round,pad=0.12',fc='white',ec='none',alpha=0.8,lw=0))
    def _draw_corr_labels(self):
        for c in CORR_TICKS:
            a=np.arccos(c); dx,dy=np.cos(a),np.sin(a)
            ha='left' if dx>0.5 else 'center'; va='bottom' if dy>0.5 else 'center'
            self.ax.annotate(str(float(c)),xy=(a,self.smax),xytext=(5*dx,5*dy),textcoords='offset points',
                ha=ha,va=va,fontsize=FONT_S,annotation_clip=False)
        a=np.pi/4
        self.ax.annotate('Correlation',xy=(a,self.smax),xytext=(40*np.cos(a),40*np.sin(a)),textcoords='offset points',
            ha='center',va='center',rotation=-45,rotation_mode='anchor',fontsize=FONT_S,annotation_clip=False)
    def _cosmetics(self):
        for loc in self._ax.axis.values():
            try:
                loc.major_ticklabels.set_fontsize(FONT_S); loc.major_ticklabels.set_pad(3)
            except Exception: pass
        self._ax.axis['left'].label.set_fontsize(FONT_S)
        self._draw_corr_labels()

def nice_ceiling(x):
    if not np.isfinite(x) or x<=0: return 1.0
    e=np.floor(np.log10(x)); frac=x/10**e
    for m in (1,2,2.5,5,10):
        if frac<=m: return m*10**e
    return 10*10**e
def round_axis_ceiling(x,n=5):
    if not np.isfinite(x) or x<=0: return n,1.0
    e=int(np.floor(np.log10(x)))
    cands=sorted(m*10**d for d in range(e-2,e+3) for m in (1,2,5))
    step=next((s for s in cands if s*n>=x),cands[-1])
    return step*n,step

with open(os.path.join(BASE,'pickles','extreme_p95_daily.pkl'),'rb') as f:
    daily=pickle.load(f)['daily']
with open(os.path.join(BASE,'pickles','extreme_p95.pkl'),'rb') as f:
    ids=pickle.load(f)['station_ids']

def per_year_p95(sid,kind,tag):
    oc=daily[sid+'|_obs_|'+kind]; od=oc['daily']; yrs=oc['years']
    mc=daily.get(sid+'|'+tag+'|'+kind)
    if mc is None: return None
    o_yr,m_yr=[],[]
    for y in yrs:
        om=od[od.index.year==y].dropna(); mm=mc[mc.index.year==y].dropna()
        if len(om)<5 or len(mm)<5: continue
        o_yr.append(np.percentile(om.values,95)); m_yr.append(np.percentile(mm.values,95))
    return (np.array(o_yr),np.array(m_yr))

obs_ann={}
for sid in ids:
    for kind in ('t','pr'):
        oc=daily[sid+'|_obs_|'+kind]; od=oc['daily']; yrs=oc['years']
        v=[np.percentile(od[od.index.year==y].dropna().values,95) for y in yrs if len(od[od.index.year==y].dropna())>=5]
        obs_ann[(sid,kind)]=np.array(v)

FIG_W=7.087; FIG_H=3.6
fig=plt.figure(figsize=(FIG_W,FIG_H))
x0=0.065; panel_w=0.30; panel_gap=0.10
panel_h=panel_w*FIG_W/FIG_H; panel_y=0.24
VARS=[('t','(a) T p95','°C'),('pr','(b) P p95','mm/day')]
for i,(kind,title,unit) in enumerate(VARS):
    stats={}
    for tag in MODELS:
        o_a=np.concatenate([per_year_p95(s,kind,tag)[0] for s in ids if per_year_p95(s,kind,tag)])
        m_a=np.concatenate([per_year_p95(s,kind,tag)[1] for s in ids if per_year_p95(s,kind,tag)])
        stats[tag]=(m_a.std(), pearsonr(o_a,m_a)[0])
    refstd=float(np.std(np.concatenate([obs_ann[(s,kind)] for s in ids])))
    max_std,step=round_axis_ceiling(max([refstd]+[s for s,_ in stats.values()])*1.15)
    fmt='%.0f' if step>=1 else ('%.1f' if step>=0.1 else '%.2f')
    levels=[nice_ceiling(max_std/5)*(j+1) for j in range(4)]
    rect=[x0+i*(panel_w+panel_gap),panel_y,panel_w,panel_h]
    dia=TaylorDiagram(refstd,fig,rect,max_std=max_std,unit=unit)
    ct=dia.add_contours(levels,colors='0.55',linewidths=0.9,zorder=0)
    dia.label_contours(ct,fmt=fmt,corr=0.55)
    for tag in MODELS:
        std,corr=stats[tag]
        dia.add_sample(std,corr,marker=MARKERS[tag],ms=7,color=COLORS[tag],mec=COLORS[tag],mew=0.8,zorder=3)
    dia._cosmetics()
    fig.text(rect[0]+rect[2]/2,panel_y+panel_h+0.135,title,ha='center',va='bottom',fontsize=FONT)
handles=[plt.Line2D([],[],marker='*',ls='',ms=12,color='black',mec='black',mew=0.8,label='observations')]+[
    plt.Line2D([],[],marker=MARKERS[m],ls='',ms=7,color=COLORS[m],mec=COLORS[m],mew=0.8,label=LABELS[m]) for m in MODELS]
fig.legend(handles=handles,loc='lower center',ncol=len(handles),bbox_to_anchor=(0.5,0.0),frameon=False,
    handlelength=1.0,handletextpad=0.4,columnspacing=1.2,fontsize=FONT_S)
out=os.path.join(BASE,'figures','extremes_validation','p95_taylor_test.png')
fig.savefig(out,dpi=600); print('saved',out)
