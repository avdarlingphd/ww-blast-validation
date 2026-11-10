#!/usr/bin/env python3
"""Render docs/master_methods_figure.png - the one-figure overview of the BLAST
false-positive validation pipeline (build -> score -> validate reads -> re-aggregate).
Read-level validation v4 (workflow/scripts/false_positive_detection.py --taxonomy; spec docs/FINAL_SPEC_v4.md):
hits are matched to the target by NCBI taxid rolled up to species (names + synonym table only as a fallback).
Each read is EXCLUDED (no named hit), CONFIRMED (target within the near-tie margin of the top named bitscore,
%id>=80; uncultured ignored), EXCLUDED_COMPLEX (every top named hit is a sibling in the target's NCBI complex),
or MISASSIGNED. Per (sample, taxon): validated fraction = CONFIRMED / all reads with a hit; corrected abundance =
validated fraction x full Kraken count; supported if CONFIRMED > 2.
Schematic/illustrative numbers.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

fig, ax = plt.subplots(figsize=(14.5, 17.9))
ax.set_xlim(0, 14.5); ax.set_ylim(0, 17.8); ax.axis("off")
PAE="#AED4F0"; CHO="#A9DCCF"; KPN="#F6CDA0"; GRY="#DBDBDB"; DEC="#DCE6F5"; GRN="#2e7d32"; RED="#c0392b"
def rb(x,y,w,h,fc,ec="#555",lw=1):
    ax.add_patch(FancyBboxPatch((x,y-h),w,h,boxstyle="round,pad=0.02,rounding_size=0.06",lw=lw,edgecolor=ec,facecolor=fc))
def arrow(x0,y0,x1,y1,label=None,lx=0,ly=0,head=True,lw=1.8):
    ax.annotate("",xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle="-|>" if head else "-",lw=lw,color="#333"))
    if label: ax.text((x0+x1)/2+lx,(y0+y1)/2+ly,label,fontsize=12,fontweight="bold",ha="center",va="center",
                      bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none"))

# ===================== PANEL A : build =====================
ax.add_patch(FancyBboxPatch((0.3,12.9),6.5,4.7,boxstyle="round,pad=0.04,rounding_size=0.1",lw=1.2,edgecolor="#999",facecolor="#FBFBFB"))
ax.text(3.55,17.22,"A.  Build the detection (per sample)",ha="center",fontsize=13,fontweight="bold")
ax.text(3.55,16.82,"Kraken assigns each read one taxon",ha="center",fontsize=9.0,style="italic",color="#555")
for i,(rd,org,col) in enumerate([("read r1","P. aeruginosa (287)",PAE),("read r2","V. cholerae (666)",CHO),("read r3","K. pneumoniae (573)",KPN)]):
    yy=16.46-i*0.46
    rb(0.55,yy,1.0,0.36,"#EFEFEF"); ax.text(1.05,yy-0.18,rd,ha="center",va="center",fontsize=8.1)
    ax.annotate("",xy=(1.85,yy-0.18),xytext=(1.6,yy-0.18),arrowprops=dict(arrowstyle="-|>",lw=0.9,color="#777"))
    rb(1.9,yy,2.2,0.36,col); ax.text(2.03,yy-0.18,org,ha="left",va="center",fontsize=8.4)
ax.text(4.35,16.0,"Each tag =\ntaxon | pair end | read-id",ha="left",va="top",fontsize=8.2,color="#444")
ax.text(2.6,15.02,"↓  subsample each read set to ≤ 500 reads (seqtk -s 42)",ha="center",fontsize=7.3,color="#333",fontweight="bold")
ax.text(2.6,14.78,"↓  tag each read by its taxon",ha="center",fontsize=8.2,color="#333",fontweight="bold")
x=0.6
for t,col in [(">287|R1|r1",PAE),(">666|R1|r2",CHO),(">573|R2|r3",KPN)]:
    rb(x,14.5,1.5,0.36,col); ax.text(x+0.75,14.32,t,ha="center",va="center",fontsize=7.9); x+=1.6
arrow(2.6,14.14,2.6,13.9,head=True,lw=1.4)
for k in (2,1): ax.add_patch(FancyBboxPatch((1.0+0.13*k,12.95+0.08*k),3.9,0.62,boxstyle="round,pad=0.02,rounding_size=0.06",lw=1,edgecolor="#bbb",facecolor="#EDE7F2"))
ax.add_patch(FancyBboxPatch((1.0,12.95),3.9,0.62,boxstyle="round,pad=0.02,rounding_size=0.06",lw=1.4,edgecolor="#5b4a7a",facecolor="#E6DAF0"))
ax.text(2.95,13.26,"ONE megablast vs nt\nper sample",ha="center",va="center",fontsize=9.4,fontweight="bold")
ax.text(5.85,13.26,"× 250\nsamples\nin parallel",ha="center",va="center",fontsize=8.4,color="#5b4a7a",fontweight="bold")

# ===================== PANEL B : score =====================
ax.add_patch(FancyBboxPatch((7.0,12.9),7.2,4.7,boxstyle="round,pad=0.04,rounding_size=0.1",lw=1.2,edgecolor="#999",facecolor="#FBFBFB"))
ax.text(10.6,17.22,"B.  Score each read  (by NCBI taxid)",ha="center",fontsize=13,fontweight="bold")
rb(7.3,16.72,2.2,0.4,CHO); ax.text(8.4,16.52,">666|R1|r2  (150 bp)",ha="center",va="center",fontsize=8.4)
cols=["subject organism","NCBI taxid","% id","bitscore"]; cw=[2.3,1.05,0.75,0.9]
data=[("Vibrio cholerae","666","100","267"),("Uncultured bacterium","77133","98","141")]
tx0,ty0=7.3,16.14; rh=0.44; cx=tx0
for j,c in enumerate(cols):
    ax.add_patch(Rectangle((cx,ty0-rh),cw[j],rh,facecolor="#E8E8E8",edgecolor="#999",lw=0.8)); ax.text(cx+cw[j]/2,ty0-rh/2,c,ha="center",va="center",fontsize=8.3,fontweight="bold"); cx+=cw[j]
for i,row in enumerate(data):
    cy=ty0-rh*(i+1); cx=tx0; best=(i==0)
    for j,val in enumerate(row):
        ax.add_patch(Rectangle((cx,cy-rh),cw[j],rh,facecolor="#CDEAD2" if best else "white",edgecolor="#999",lw=0.8)); ax.text(cx+cw[j]/2,cy-rh/2,val,ha="center",va="center",fontsize=8.4,fontweight="bold" if best else "normal"); cx+=cw[j]
ax.text(12.35,15.5,"best hit =\nhighest\nbitscore",ha="left",va="center",fontsize=8.3,color="#333")
ax.text(10.6,14.5,"Target = any hit whose taxid rolls up to the target species (renames such as\nRaoultella → Klebsiella cannot break it); a hit with no taxid is matched by the synonym table",
        ha="center",va="center",fontsize=8.0,color="#333")
ax.text(10.6,14.02,"Near-tie:  target within the margin of the top named score (2, 5 or 10%: pending)",ha="center",fontsize=8.3,fontweight="bold",color=GRN)
rb(7.3,13.78,2.75,0.4,CHO); ax.text(8.68,13.58,"V. cholerae (141)  target",ha="center",va="center",fontsize=8.0)
ax.text(10.35,13.58,"within 2% of",ha="center",va="center",fontsize=8.3,fontweight="bold",color="#333")
rb(11.3,13.78,2.75,0.4,"#FBF0D8"); ax.text(12.68,13.58,"V. mimicus (143)  top named",ha="center",va="center",fontsize=8.0)

arrow(7.0,15.3,6.8,15.3,lw=1.9)   # A -> B

# ===================== PANEL C : validate each read =====================
ax.text(0.3,12.6,"C.  Validate each read, then re-aggregate",fontsize=13,fontweight="bold")
SX,SW=0.6,5.6; LX,LW=8.2,5.6
LEAF={"CONFIRMED":"#CDEAD2","MISASSIGNED":"#F5CBC4","EXCLUDED":"#E6E6E6","EXCLUDED_COMPLEX":"#E1E7F0"}
def dec(y,h,txt,fs=11.0): rb(SX,y+h/2,SW,h,DEC); ax.text(SX+SW/2,y,txt,ha="center",va="center",fontsize=fs)
def leaf(y,name,sub): rb(LX,y+0.5,LW,1.0,LEAF[name]); ax.text(LX+LW/2,y+0.18,name,ha="center",va="center",fontsize=13.0,fontweight="bold"); ax.text(LX+LW/2,y-0.22,sub,ha="center",va="center",fontsize=8.3,style="italic",color="#444")
d1,d2,d3=11.4,9.7,7.85; hb=1.0; h2=1.5; h3=1.3; mis=6.3
ax.text(3.4,12.32,"each read Kraken assigned to the target taxon",ha="center",fontsize=8.8,style="italic",color="#555")
arrow(3.4,12.1,3.4,d1+hb/2,head=True)
dec(d1,hb,"Any hit to a NAMED organism in nt?\n(no nt hit, or uncultured-only  =  excluded)")
dec(d2,h2,"Is a TARGET-species hit (by taxid) within the margin\nof the top NAMED bitscore?  (uncultured ignored)\nAND     %id ≥ 80%",fs=10.4)
dec(d3,h3,"Are ALL top-scoring named hits sibling species\nin the target's NCBI complex?  (e.g. K. pneumoniae complex)",fs=10.0)
leaf(d1,"EXCLUDED","no named hit\n(no nt match, or uncultured-only)")
leaf(d2,"CONFIRMED","counts toward the taxon")
leaf(d3,"EXCLUDED_COMPLEX","a same-complex sibling wins:\nnot evidence against the taxon")
leaf(mis,"MISASSIGNED","an organism outside the complex outscores\nthe target   (false-positive read)")
cx=SX+SW/2
arrow(SX+SW,d1,LX,d1,"No")
arrow(cx,d1-hb/2,cx,d2+h2/2,"Yes",lx=0.32)
arrow(SX+SW,d2,LX,d2,"Yes")
arrow(cx,d2-h2/2,cx,d3+h3/2,"No",lx=0.32)
arrow(SX+SW,d3,LX,d3,"Yes")
arrow(cx,d3-h3/2,cx,mis,"No",lx=0.32,head=False); arrow(cx,mis,LX,mis)
# ----- re-aggregate: validated fraction -> corrected abundance -----
ax.add_patch(FancyBboxPatch((0.6,3.85),13.3,1.75,boxstyle="round,pad=0.03,rounding_size=0.08",lw=1.3,edgecolor="#5b4a7a",facecolor="#F3EFF8"))
ax.text(0.9,5.36,"Re-aggregate per (sample, taxon):",fontsize=10.6,fontweight="bold",color="#4a3a66")
for yy,txt in [
 (5.07,"validated fraction  =  CONFIRMED / all reads with a BLAST hit  (incl. EXCLUDED, EXCLUDED_COMPLEX)     (≤ 500-read subsample)"),
 (4.82,"     alternative (choice pending):  CONFIRMED / (CONFIRMED + MISASSIGNED)   (EXCLUDED reads do not count against the taxon)"),
 (4.57,"corrected abundance  =  validated fraction  ×  reads_clade        (full Kraken clade read count)"),
 (4.32,"relative abundance  =  corrected abundance  /  total QC'd + dehosted reads        (sample library size)"),
 (4.07,"supported detection  if  CONFIRMED reads > 2"),
]:
    ax.text(1.05,yy,txt,fontsize=9.2,color="#222")

# ===================== worked examples =====================
# ---- Example 1: one read ----
ax.add_patch(FancyBboxPatch((0.4,0.25),6.55,3.4,boxstyle="round,pad=0.03,rounding_size=0.08",lw=1.2,edgecolor=GRN,facecolor="#F6FBF7"))
ax.text(0.62,3.42,"Example 1:   one read   >666|R1|r2",fontsize=10.6,fontweight="bold")
e1hx=[0.55,2.85,3.65,4.35,5.20]; e1hw=[2.30,0.80,0.70,0.85,1.40]
for j,h in enumerate(["BLAST hit organism","taxid","% id","bitscore","used?"]):
    ax.add_patch(Rectangle((e1hx[j],2.56),e1hw[j],0.36,facecolor="#E8E8E8",edgecolor="#999",lw=0.8)); ax.text(e1hx[j]+e1hw[j]/2,2.74,h,ha="center",va="center",fontsize=7.8,fontweight="bold")
e1rows=[("Uncultured bacterium","77133","99","145","ignored","#EEEEEE"),
        ("Vibrio mimicus  (relative)","674","98","143","top named","#FBF0D8"),
        ("Vibrio cholerae  (target)","666","98","141","within 2%","#D7EFDD")]
for i,(org,tid,pid,bit,note,bg) in enumerate(e1rows):
    ry=2.56-0.40*(i+1); cyc=ry+0.20
    for j in range(5):
        ax.add_patch(Rectangle((e1hx[j],ry),e1hw[j],0.40,facecolor=bg,edgecolor="#D5E3D8",lw=0.7))
    ax.text(e1hx[0]+0.07,cyc,org,ha="left",va="center",fontsize=7.6,fontweight="bold" if i==2 else "normal")
    ax.text(e1hx[1]+e1hw[1]/2,cyc,tid,ha="center",va="center",fontsize=7.8)
    ax.text(e1hx[2]+e1hw[2]/2,cyc,pid,ha="center",va="center",fontsize=7.8)
    ax.text(e1hx[3]+e1hw[3]/2,cyc,bit,ha="center",va="center",fontsize=7.8)
    ax.text(e1hx[4]+e1hw[4]/2,cyc,note,ha="center",va="center",fontsize=7.3,style="italic",color="#555")
ax.text(0.62,1.04,"Ignore uncultured.  Top NAMED hit = V. mimicus (143).",fontsize=7.9,color="#222")
ax.text(0.62,0.76,"Target (taxid 666 = V. cholerae) scores 141: within 2% of 143,  %id ≥ 80.",fontsize=7.9,color="#222")
ax.text(0.62,0.44,"→  CONFIRMED   (counts toward V. cholerae)",fontsize=8.7,fontweight="bold",color=GRN)
# ---- Example 2: is the (sample, taxon) supported? ----
ax.add_patch(FancyBboxPatch((7.15,0.25),6.95,3.4,boxstyle="round,pad=0.03,rounding_size=0.08",lw=1.2,edgecolor="#5b4a7a",facecolor="#F6F4FB"))
ax.text(7.4,3.42,"Example 2:   is (sample S12, V. cholerae) supported?",fontsize=10.2,fontweight="bold")
for yy,txt in [
 (2.92,"reads in subsample:    CONFIRMED 418   |   MISASSIGNED 60   |   EXCLUDED 22"),
 (2.54,"validated fraction  =  418 / (418 + 60 + 22)  =  83.6%"),
 (2.16,"reads_clade (Kraken)  =  4,200"),
 (1.78,"corrected abundance  =  0.836 × 4,200  =  3,511 reads"),
 (1.40,"total QC'd + dehosted reads  =  18,540,000"),
 (1.02,"relative abundance  =  3,511 / 18,540,000  =  0.0189%"),
 (0.56,"supported?    418 CONFIRMED > 2   →   yes   ✓"),
]:
    ax.text(7.4,yy,txt,fontsize=8.5,color="#222")

plt.savefig("docs/master_methods_figure.png",dpi=185,bbox_inches="tight",facecolor="white")
print("saved")
