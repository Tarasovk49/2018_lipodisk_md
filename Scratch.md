# This is scratch document for algorithms, plots, logs and everything else that can be useful.
### Prepare lipodisk with sensory rhodopsin in DMPC.
1. Prepare topology for DMPC. Cut down the DPPC by two atoms on each chain - C215, C216, C315, C316. Delete all bonded interactions from *lipids.rtp*.
2. Prepare topology for Retinal-lysine (residue REK).
* Upload retinal without C15 linkage to lysine to [TTPMKTOP](http://erg.biophys.msu.ru/erg/tpp/). Manually change some atom names according to alkane and alkene atomname specifications.
* Manually modify *ffbonded.itp*. Add specific angles and dihedrals which are absent in normal opls-aa forcefield. The parameters for such angles and dihedrals were chosen corresponding to alkane and alkene topology specifications. The lines added can be seen below.

**Angles**
```
; Retinal+lysine
  C=      C     NC      1   116.000    585.760   ; wlj
  CT     NC      C      1   118.600    585.760   ; 
  C=     C=     C=      1   124.000    585.760   ; wlj
  C=     CT     CT      1   116.000    585.760   ; wlj
  C=     C=      C      1   124.000    585.760   ; wlj
```
**Dihedrals**
```
; Retinal+lysine
  HC     CT     CT     NC      3      5.77183  -2.67148   0.95814  -4.05848   0.00000   0.00000 ; 
  CT     CT     NC     C       3      3.80117  -6.95172  -1.01671   4.16726   0.00000   0.00000 ; 
  HC     CT     NC     C       3      3.80117  -6.95172  -1.01671   4.16726   0.00000   0.00000 ; JUST COPIED FROM C-C-N-C!
  CT     NC      C     C=      3     58.57600   0.00000 -58.57600   0.00000   0.00000   0.00000 ; imine 
  CT     NC      C     HC      3      3.80117  -6.95172  -1.01671   4.16726   0.00000   0.00000 ; JUST COPIED FROM C-C-N-C!
  CT     CM     C=     C=      3      2.92880  -1.46440   0.20920  -1.67360   0.00000   0.00000 ; hydrocarbon all-atom
  C=     CM     CT     HC      3      0.62760   1.88280   0.00000  -2.51040   0.00000   0.00000 ;
  CM     C=     C=     C=      3     58.57600   0.00000 -58.57600   0.00000   0.00000   0.00000 ; alkene
  CT     C=     C=     C=      3     58.57600   0.00000 -58.57600   0.00000   0.00000   0.00000 ; alkene
  CM     C=     CT     CT      3      2.92880  -1.46440   0.20920  -1.67360   0.00000   0.00000 ; hydrocarbon all-atom
  C=     C=     CT     CT      3      2.92880  -1.46440   0.20920  -1.67360   0.00000   0.00000 ; hydrocarbon all-atom
  C=     C=     C=     HC      3     58.57600   0.00000 -58.57600   0.00000   0.00000   0.00000 ; alkene
  C=     C=     C=     C=      3     58.57600   0.00000 -58.57600   0.00000   0.00000   0.00000 ; alkene
  C=     CT     CT     HC      3      0.62760   1.88280   0.00000  -2.51040   0.00000   0.00000 ;
  C=     CT     CT     CT      3      2.92880  -1.46440   0.20920  -1.67360   0.00000   0.00000 ; hydrocarbon all-atom
  C=     C=     C=     C       3     58.57600   0.00000 -58.57600   0.00000   0.00000   0.00000 ; alkene
  CT     C=     C=     C       3     58.57600   0.00000 -58.57600   0.00000   0.00000   0.00000 ; alkene
  C=     C=     CT     HC      3      0.62760   1.88280   0.00000  -2.51040   0.00000   0.00000 ;
  C=     C=     C      NC      3      5.77183  -2.67148   0.95814  -4.05848   0.00000   0.00000 ; 
  C=     C=     C      HC      3      0.62760   1.88280   0.00000  -2.51040   0.00000   0.00000 ;
  HC     C=     C      HC      3      0.62760   1.88280   0.00000  -2.51040   0.00000   0.00000 ; hydrocarbon *new* 11/99
  HC     C=     C      NC      3      5.77183  -2.67148   0.95814  -4.05848   0.00000   0.00000 ; 
```
