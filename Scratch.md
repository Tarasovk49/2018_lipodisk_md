# This is scratch document for algorithms, plots, logs and everything else that can be useful.
### Prepare lipodisk with sensory rhodopsin in DMPC.
1. Prepare topology for DMPC. Cut down the DPPC by two atoms on each chain - C215, C216, C315, C316. Delete all bonded interactions from *lipids.rtp*.
2. Prepare topology for Retinal-lysine (residue REK).
* Upload retinal without C15 linkage to lysine to [TTPMKTOP](http://erg.biophys.msu.ru/erg/tpp/).
* Modify *ffbonded.itp* and *ffnonbonded.itp*.

**ffbonded.itp**
`; Retinal+lysine
  C=      C     NC      1   116.000    585.760   ; wlj
  CT     NC      C      1   118.600    585.760   ; 
  C=     C=     C=      1   124.000    585.760   ; wlj
  C=     CT     CT      1   116.000    585.760   ; wlj
  C=     C=      C      1   124.000    585.760   ; wlj`
