import MDAnalysis as mda

lipodisk = 'lipodisk_npt.gro'
lipodisk = mda.Universe(lipodisk)
all_lipids = lipodisk.select_atoms("not resid 111 106 and name C211 C212 C213 and resname DMPC")
edge = lipodisk.select_atoms("not resid 111 106 and name C211 C212 C213 and resname DMPC and same resid as (around 10 resname MAL MAR ST1 ST2 DB1 DB2 MAD MA2)")
center = lipodisk.select_atoms("not resid 111 106 and name C211 C212 C213 and resname DMPC and not same resid as (around 10 resname MAL MAR ST1 ST2 DB1 DB2 MAD MA2)")

#center.write("center.pdb")
#edge.write("edge.pdb")

all_lipids.write("index_all_lipids.ndx", name="ALL_LIPIDS")
edge.write("index_edge.ndx", name="EDGE_LIPIDS")
center.write("index_center.ndx", name="CENTER_LIPIDS")
