import sys, collections
# Replay accepted flips on the initial face list to get the final connectivity,
# then count extraordinary corners / faces among the energy faces (copy == 0).
def run(d, upto=None):
    faces=[]
    for line in open(d+'/inputmesh_faces.csv'):
        if line.startswith('#'): continue
        v0,v1,v2,c=map(int,line.split(',')); faces.append([v0,v1,v2,c])
    bad=0; n=0
    for line in open(d+'/inputEdgeFlips.csv'):
        if line.startswith('iteration'): continue
        it,e,v0,v1,t0,t1,dE,acc=line.split(',')
        if int(acc)!=1: continue
        if upto is not None and int(it)>=upto: break
        v0,v1,t0,t1=map(int,(v0,v1,t0,t1)); n+=1
        hits=[f for f in faces if v0 in f[:3] and v1 in f[:3]]
        if len(hits)!=2: bad+=1; continue
        fa=[f for f in hits if t0 in f[:3]]; fb=[f for f in hits if t1 in f[:3]]
        if len(fa)!=1 or len(fb)!=1: bad+=1; continue
        fa=fa[0]; fb=fb[0]
        # keep orientation: replace v1 in fa by t1, v0 in fb by t0
        fa[fa.index(v1)]=t1; fb[fb.index(v0)]=t0
    val=collections.Counter()
    for f in faces:
        for v in f[:3]: val[v]+=1
    vh=collections.Counter(); 
    for f in faces:
        if f[3]==0:
            pass
    # vertex valence histogram over vertices of energy faces
    everts=set(v for f in faces if f[3]==0 for v in f[:3])
    for v in everts: vh[val[v]]+=1
    nE=reg=ev=multi=0; corners=collections.Counter(); kcalls={}
    for f in faces:
        if f[3]: continue
        nE+=1
        vs=[val[v] for v in f[:3]]
        k=[x for x in vs if x!=6]
        if not k: reg+=1
        else:
            ev+=1
            if len(k)>1: multi+=1
            for x in k: corners[x]+=1
    return dict(flips=n, bad=bad, energyFaces=nE, regular=reg, EV=ev, multiEV=multi,
                EVcorners=sum(corners.values()), corners=dict(sorted(corners.items())),
                valenceHist=dict(sorted(vh.items())), nverts=len(everts))
if __name__=='__main__':
    for d in sys.argv[1:]:
        print(d, run(d))

def kernel_calls(d, depths):
    """Force-pass kernel calls per step over the energy faces, for depths {N: D}."""
    import collections
    faces=[]
    for line in open(d+'/inputmesh_faces.csv'):
        if line.startswith('#'): continue
        v0,v1,v2,c=map(int,line.split(',')); faces.append([v0,v1,v2,c])
    for line in open(d+'/inputEdgeFlips.csv'):
        if line.startswith('iteration'): continue
        it,e,v0,v1,t0,t1,dE,acc=line.split(',')
        if int(acc)!=1: continue
        v0,v1,t0,t1=map(int,(v0,v1,t0,t1))
        hits=[f for f in faces if v0 in f[:3] and v1 in f[:3]]
        fa=[f for f in hits if t0 in f[:3]][0]; fb=[f for f in hits if t1 in f[:3]][0]
        fa[fa.index(v1)]=t1; fb[fb.index(v0)]=t0
    val=collections.Counter()
    for f in faces:
        for v in f[:3]: val[v]+=1
    K=0; const=0; corners=0
    for f in faces:
        if f[3]: continue
        vs=[val[v] for v in f[:3]]; ev=[x for x in vs if x!=6]
        if not ev: K+=1; const+=1
        elif len(ev)==1: K+=3*depths[ev[0]]; corners+=1
        else:
            K+=1+(3-len(ev)); const+=1+(3-len(ev))
            for x in ev: K+=3*depths[x]
            corners+=len(ev)
    return K, const, corners
