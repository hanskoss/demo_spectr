# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 21:42:47 2022

@author: Hans
"""
from __future__ import division
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(6,8))
ax = Axes3D(fig)
global trans
#%%

""" Stretching 3D plot
To stretch the 3D plot, some fiddling is required.
See the discussion at
https://stackoverflow.com/questions/30223161/
The following lines are taken from that discussion.
Definition of some global variables, so this might be considered "messy"

Set the aspect ratio for the 3D plots here. Careful, correction parameters
will have to be determined (set further below in the script where needed).
The problem is that the operation not only stretches but also translates
the the graph. At this stage is easiest to adjust that tranlation by
empirically determining the required x and y translations for the graph,
setting the xtransl and ytransl parameters below. In order to find a
"lost" graph, zoom can be set to a large number such as 100.
"""
x_scale=1;y_scale=1;z_scale=2
scale=np.diag([x_scale, y_scale, z_scale, 1.0])
scale=scale*(1.0/scale.max())
scale[3,3]=1.0

def short_proj():
    """
    Definition necessary for stretching the 3D plot (see above)
    """    
    return np.dot(Axes3D.get_proj(ax), scale)+transl

def plt3d(axobj,azim,elev,xtransl,ytransl,zoom):
    """3d plot helper function, setting the viewing angles,
    setting the axis direction right and determining trans
    which has to be assigned to the variable transl.
    Has to do with the figure stretching (see above)
    """
    axobj.azim=azim
    axobj.elev=elev
    #axobj.roll=roll #expected in matplotlib 3.6
    trans=np.array([[0,0,0,xtransl],[0,0,0,ytransl],[0,0,0,0],[0,0,0,zoom]])
    axobj.set_xlim(np.sort(list(ax.get_xlim()))[::-1])
    axobj.set_ylim(np.sort(list(ax.get_ylim()))[::1])
    axobj.set_zlim(np.sort(list(ax.get_zlim()))[::1])
    return trans

""" Read in assignment data and sequence
The specific columns= choice might have to be corrected depending on the specific
file format / file.
The file header was cut off in this example.
Columns in this example script are:
    col. 5: residue number
    col. 7: residue name
    col. 8: specific nucleus (e.g., CA)
    col. 11: chemical shift
    col. 9: nucleus type (e.g., H)
The chemical shifts are obtained in the format [CS_res1,CS_res2,CS_res3,...]
with CS_res1 = [chemical shift 1 for res. 1, chemical shift 2 for res.1, ...
                chemical shift 3 for res 1] and corresponding for CS_res2 etc.
each chemical shift n for residue m is defined as:
[specific nucleus (string), chemical shift (float), type (str)]
"""

path='C:\\Users\\Hans\\Desktop\\job applications\\interviews\\edinburgh'
os.chdir(path)
file=open('27095_cut.txt')
csvreader=csv.reader(file, delimiter=' ')
cshifts=[[]];seq=[];init=1
for row in csvreader:
    rowshort=[i for i in row if i != '']
    if init == 1:
        j=rowshort[5]; init=0
        labn=rowshort[7]
    if int(rowshort[5]) > j:
        seq.append([j,labn])
        j=int(rowshort[5])
        labn=rowshort[7]
        if cshifts != [[]]:
            cshifts.append([])
    cshifts[-1].append([rowshort[8],float(rowshort[11]),rowshort[9]])
seq.append([j,rowshort[7]])

"""
The chemical shifts (cshifts) list is inconvenient, we prefer a 3D array.
This is a simple loop to reorder this. For each chemical shift group
"""
points=[]
hsqc=[]
for ii,i in enumerate(cshifts):
    pointgroup=[]
    nseq=[]
    hurdl=0
    for j in i:
        if j[2] == 'H' and j[0] == 'H':
            hsav=j[1]
            hurdl=1
        if hurdl > 0 and j[2] == 'H':
            pointgroup.append([hsav,j[1]])
            hurdl=2
        if hurdl == 2 and j[2] == 'N' and j[0] == 'N':
            nseq.append(j[2])
            for point in pointgroup:
                point.append(j[1])
                points.append(point)
            hsqc.append([seq[ii][0],seq[ii][1],hsav,j[1]])

"""creates about 2000 random cross peaks"""
newpoints=[]
for i in points:
    for j in points:
        if i == j:
            if i[0] == i[1]:
                newpoints.append([i[0],i[1],i[2],1])
            else:
                newpoints.append([i[0],i[1],i[2],(1-np.random.rand(1)/4)[0]])
        else:
            if np.random.randint(1000) == 44:
                newpoints.append([j[0],i[1],j[2],np.random.rand(1)[0]])
    

"""sorts the coordinates according to the order provided in the nuclei list"""
points=np.array(newpoints).T

#%%

ax.clear()
ax.set_proj_type('ortho')
x=points[0]
y=points[1]
z=points[2]
intensity=points[3]
z2=np.ones(shape=x.shape)*max(z)
y2=np.ones(shape=x.shape)*max(y)
dimch=2
#x2,y2,z2=[np.array([l if n != dimch else np.max(points[n]) for l in m]) for n,m in enumerate(points)]



xtransl=[-2.5,0,-2.5]
ytransl=[-0.65,1.4,-1.5]
zoom=[-7.5,2.5,-5]
azim=[-90-0.1,-25,-90-0.1]
elev=[-90-0.1,-50,180-0.1]

#frac2=-0.1
#frac2=0.001
#frac2+=0.1
#frac2+=0.1
frac2=0.5
print(frac2)
sel1=frac2-1
sel2=frac2
if frac2 <=1:
    sel1=0
    sel2=1
    frac=frac2
    xnew=x#z    
    znew=y
    ynewr=z
    ynew=np.array([np.min(z) for i in z])#x
else:
    sel1=1
    sel2=2
    frac=frac2-1
    xnew=x#z    
    znew=y
    ynew=z#x

xtransl=[-2.5-0.1,-1.7-0.1,-2.5-0.1]
ytransl=[0.2-0.6,-4.7-0.6,-4.7-0.6]
zoom=[-8,-8,-8] #-8
azim=[-90,-80,-90]
elev=[-180,-110,-90]
roll=[0,100,0]
ax.xaxis._axinfo['juggled']=(2,0,1)
ax.yaxis._axinfo['juggled']=(2,1,0)
#%%
frac3=0
frac=0.0001
sel1=0
sel2=1

#def update(frac,frac3):
fnlist=[]
for fn,progress in enumerate(np.linspace(0,1,1000)):
#for fn in [999]:
    fno=str(fn)+'.png'
    x0=0
    if 100 <= fn < 200:
        frac=(fn-100)/100
        if frac == 0:
            frac=0.0001
        x0=0
    if fn == 200:
        frac=1
    if 700 <= fn < 800:
        frac3=(fn-700)/100
        x0=0
    if 900 <= fn < 1000:
        sel1=1
        sel2=2
        xnew=x#z    
        znew=y
        znewr=np.min(y)
   #     ynew=z
        frac=(fn-900)/100
        frac3=(fn-900)/100
        x0=0
    if fn == 0 or fn == 800 or x0 == 1 or fn == 1000:
        ax.clear()
        ax.set_xlabel('$\delta(^1HN-H)$ / ppm')
        ax.set_ylabel('$\delta(^{15}N)$ / ppm')
        ax.set_zlabel('$\delta(^1H)$ / ppm')
        if fn < 110: 
            ax.set_yticks([])
            ax.set_ylabel('')
        elif fn < 150:
            ax.set_ylabel('')
        elif fn > 975:
            ax.set_zticks([])
            ax.set_zlabel('')
        elif fn > 950:
            ax.set_zlabel('')

        if fn < 900:
            sc = ax.scatter(xnew, ynewr, znew, s=intensity*0, marker='o',c=intensity, alpha=1, cmap='OrRd',facecolors='none')
            sc = ax.scatter(xnew, ynew+(ynewr-ynew)*frac3, znew, s=intensity*10, marker='o',c=(ynewr-ynew)*(0.0001+frac3)/np.max(ynewr-ynew),norm=plt.Normalize(vmin=0,vmax=0.8), alpha=1, cmap='plasma',facecolors='none')
        else:
            sc = ax.scatter(xnew, ynewr, znew, s=intensity*0, marker='o',c=intensity, alpha=1, cmap='OrRd',facecolors='none')
            sc = ax.scatter(xnew, ynewr, znew-(znew-znewr)*frac3, s=intensity*10, marker='o',c=(ynewr-ynew)*(1-frac3)/np.max(ynewr-ynew),norm=plt.Normalize(vmin=0,vmax=0.8), alpha=1, cmap='plasma',facecolors='none')
        
        transl=plt3d(ax,azim[sel1]*(1-frac)+azim[sel2]*frac,elev[sel1]*(1-frac)+elev[sel2]*frac,xtransl[sel1]*(1-frac)+xtransl[sel2]*frac,ytransl[sel1]*(1-frac)+ytransl[sel2]*frac,zoom[sel1]*(1-frac)+zoom[sel2]*frac)
        ax.get_proj=short_proj
        plt.savefig(fno, dpi=500)
        plt.draw()
    
    
    fnlist.append(fno)

#%%
import ffmpeg
print('hallo')
#os.system("ffmpeg -i 1%2d.png -r 20 video1xx7.avi")
os.system("ffmpeg -i 7%2d.png -r 20 video7xx2.avi")
os.system("ffmpeg -i 9%2d.png -r 20 video9xx2.avi")
        
#%%
from adjustText import adjust_text

sl=['A','S','G','P','N','D','A','Q','E','Y','F','W','R','K','H','L','I','V','C','M','T']
tl=['ALA','SER','GLY','PRO','ASN','ASP','ALA','GLN','GLU','TYR','PHE','TRP','ARG','LYS','HIS','LEU','ILE','VAL','CYS','MET','THR']
def three2one(nam):
    return sl[tl.index(nam)]
fig2=plt.figure(4)
plt.clf()
plt.scatter(xnew,ynewr,color='red',s=3)
texts=[]
for x,y,s,s2 in zip(list(np.array(hsqc).T[2]),list(np.array(hsqc).T[3]),list(np.array(hsqc).T[0]),list(np.array(hsqc).T[1])):
    #print(x,y,s,s2)
    texts.append(plt.text(float(x),float(y),str(s)+str(three2one(s2)),fontsize=8))
plt.xlim(plt.xlim()[::-1])
plt.ylim(plt.ylim()[::-1])
plt.xlabel('$\delta(^1HN-H)$ / ppm')
plt.ylabel('$\delta(^{15}N)$ / ppm')
plt.tight_layout()
adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5))
plt.draw()
#%%


#%%
yl=list(ax.get_ylim())
ynew2=np.array([i for i in z])#x
plt.clf()
zoom=[-8,1000,-8]
transl=plt3d(ax,azim[sel1]*(1-frac)+azim[sel2]*frac,elev[sel1]*(1-frac)+elev[sel2]*frac,xtransl[sel1]*(1-frac)+xtransl[sel2]*frac,ytransl[sel1]*(1-frac)+ytransl[sel2]*frac,zoom[sel1]*(1-frac)+zoom[sel2]*frac)
sc = ax.scatter(xnew, ynew2, znew, s=intensity*10, marker='o',c='black', alpha=1, cmap='coolwarm',facecolors='none')
#ax.set_ylim(set(yl))
plt.draw()

#%%
import imageio
with imageio.get_writer('mygif.gif', mode='I') as writer:
    for filename in fnlist:
        image = imageio.imread(filename)
        writer.append_data(image)
#%%
frac+=1/nsteps


#%%



#%%

for ii in np.arange(0,360,1):
    ax.view_init(elev=10., azim=ii)
    savefig("movie%d.png" % ii)