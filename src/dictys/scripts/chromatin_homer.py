#!/usr/bin/env python3
# Nikolaos Trasanidis, Lingfei Wang, 2022. All rights reserved.

import numpy as np
import pandas as pd
d=pd.read_csv("15-tf.bed",header=0,index_col=None,sep='\t')
d=d[d['MotifScore']>0].copy()
d['chr']=d['PositionID'].apply(lambda x:x.split(':')[1])
d['start']=d['PositionID'].apply(lambda x:int(x.split(':')[2]))
d['end']=d['PositionID'].apply(lambda x:int(x.split(':')[3]))
d.sort_values(['chr','start','end','Motif Name','MotifScore'],inplace=True)
#Keep site with strongest Homer score
d.drop_duplicates(subset=['chr','start','end','Motif Name'],keep='last',inplace=True)
d['w']=d['PositionID'].apply(lambda x:float(x.split(':')[4])).abs()
d['name']=d['PositionID'].apply(lambda x:':'.join(x.split(':')[1:-1]))
d[['chr','start','end','name','w','Motif Name','MotifScore']].to_csv('16-long.bed',header=False,index=False,sep='\t')

namep=sorted(list(set(d['name'].tolist())))
namem=sorted(list(set(d['Motif Name'].tolist())))
dp=dict(zip(namep,range(len(namep))))
dm=dict(zip(namem,range(len(namem))))
idp=[dp[x] for x in d['name']]
idm=[dm[x] for x in d['Motif Name']]
answ=np.ones((len(dp),len(dm)),dtype='f8')*np.nan
answ[idp,idm]=d['w'].values
answ=pd.DataFrame(answ,index=namep,columns=namem)
answ.to_csv('19-w.tsv',header=True,index=True,sep='\t')
ansh=np.ones((len(dp),len(dm)),dtype='f8')*np.nan
ansh[idp,idm]=d['MotifScore'].values
ansh=pd.DataFrame(ansh,index=namep,columns=namem)
ansh.to_csv('19-h.tsv',header=True,index=True,sep='\t')
