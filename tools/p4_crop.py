import json, numpy as np, cv2, os
import sys
R, O = sys.argv[1:3]  # render dir (A/ and Cm/, each with ink/ and metric/), output dir
S={a:json.load(open(f'{R}/{a}/metric/metrics.json'))['strips'][0]['slices'] for a in ('A','Cm')}
rows=[]
for sa,sc in zip(S['A'],S['Cm']):
    rows.append((sa['x0'],sa['x1'],sa['fg_pixels'],sc['fg_pixels'],sa['line_score'],sc['line_score']))
with open(f'{O}/slice_table.tsv','w') as f:
    f.write('x0\tx1\tfg_A\tfg_Cm\tline_A\tline_Cm\n')
    for r in rows: f.write('\t'.join(f'{v:.3f}' if isinstance(v,float) else str(v) for v in r)+'\n')
act=[r for r in rows if r[2]+r[3]>20000]
print('slices', len(rows), 'with ink (>20k px either arm):', len(act), 'Cm>A in', sum(r[3]>r[2] for r in act))
b=max(rows,key=lambda r:r[2]); x0,x1=b[0],b[1]
print('crop slice', b)
def strip(a,kind):
    t0,t1=x0//16384,(x1-1)//16384; parts=[]
    for t in range(t0,t1+1):
        p=f'{R}/{a}/ink/w010-135_flat.{t:03d}.jpg' if kind=='ink' else f'{R}/{a}/metric/predictions/w010-135_flat_mask.{t:03d}.png'
        parts.append(cv2.imread(p,cv2.IMREAD_GRAYSCALE))
    im=np.concatenate(parts,1); off=t0*16384
    return im[:,x0-off:x1-off]
ims={}
for a in ('A','Cm'):
    ink=strip(a,'ink'); m=strip(a,'mask')
    if m.max()==1: m=m*255
    ov=cv2.cvtColor(ink,cv2.COLOR_GRAY2BGR); ov[m>127]=(0.4*ov[m>127]+0.6*np.array([40,40,230])).astype(np.uint8)
    cv2.imwrite(f'{O}/crop_{a}_ink.png',ink); ims[a]=(ink,ov)
    print(a,'crop fg px',int((m>127).sum()))
def lab(im,t):
    im=im.copy(); cv2.putText(im,t,(20,60),cv2.FONT_HERSHEY_SIMPLEX,2,(255,255,255) if im.ndim==2 else (0,255,255),4); return im
sep=np.full((12,x1-x0),128,np.uint8)
fig=np.concatenate([lab(ims['A'][0],'A raw gen-19'),sep,lab(ims['Cm'][0],'Cm snap')],0)
cv2.imwrite(f'{O}/crop_A_vs_Cm.jpg',fig,[cv2.IMWRITE_JPEG_QUALITY,92])
sep3=np.full((12,x1-x0,3),128,np.uint8)
fo=np.concatenate([lab(ims['A'][1],'A + ink mask'),sep3,lab(ims['Cm'][1],'Cm + ink mask')],0)
cv2.imwrite(f'{O}/crop_A_vs_Cm_mask.jpg',fo,[cv2.IMWRITE_JPEG_QUALITY,92])
# half-width previews for viewing
for n in ('crop_A_vs_Cm','crop_A_vs_Cm_mask'):
    im=cv2.imread(f'{O}/{n}.jpg'); h,w=im.shape[:2]
    for i in range(4): cv2.imwrite(f'{O}/{n}_part{i}.jpg', im[:, i*w//4:(i+1)*w//4])
