import pandas as pd
REGIONS = [("Cerebral cortex", ["Cerebral cortex (Cx)", "Perirhinal cortex (area 35) (A35)"]),
           ("Paleocortex and claustrum", ["Paleocortex (PalCx)", "Claustrum"]),
           ("Hippocampus", ["Head of hippocampus (HiH)", "Body of hippocampus (HiB)", "Tail of Hippocampus (HiT)"]),
           ("Amygdala", ["Amygdaloid complex (AMY)", "Extended amygdala (EXA)"]),
           ("Basal nuclei and basal forebrain", ["Basal nuclei (BN)", "Basal forebrain (BF)"]),
           ("Thalamus and epithalamus", ["Thalamus (THM)", "Epithalamus"]),
           ("Hypothalamus", ["Hypothalamus (HTH)"]),
           ("Midbrain", ["Midbrain (M)", "Midbrain (RN)"]),
           ("Pons", ["Pons (Pn)"]),
           ("Medulla", ["Myelencephalon (medulla oblongata) (Mo)"]),
           ("Cerebellum", ["Cerebellum (CB)"]),
           ("Spinal cord", ["Spinal cord"])]
SCS = sorted(pd.read_parquet('sc_map.parquet').sc.unique())
