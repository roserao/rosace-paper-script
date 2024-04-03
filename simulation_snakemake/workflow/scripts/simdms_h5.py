import pandas as pd
import os

file_path = "results/simdms/binding_simulation_depth_200/binding.h5"
effect = pd.read_hdf(file_path, key = "effects")
expected = pd.read_hdf(file_path, key = "expected")
popcounts = pd.read_hdf(file_path, key = "popcounts")
seqcounts = pd.read_hdf(file_path, key = "seqcounts")

dir_path = "results/simdms/binding_simulation_depth_200/Tsv"
os.mkdir(dir_path)

effect.to_csv(os.path.join(dir_path, "effects.tsv"), sep='\t', index=False)
expected.to_csv(os.path.join(dir_path, "expected.tsv"), sep='\t', index=False)
popcounts.to_csv(os.path.join(dir_path, "popcounts.tsv"), sep='\t', index=False)
seqcounts.to_csv(os.path.join(dir_path, "seqcounts.tsv"), sep='\t', index=False)


file_path = "results/simdms/growth_simulation_depth_200/growth.h5"
effect = pd.read_hdf(file_path, key = "effects")
expected = pd.read_hdf(file_path, key = "expected")
popcounts = pd.read_hdf(file_path, key = "popcounts")
seqcounts = pd.read_hdf(file_path, key = "seqcounts")

dir_path = "results/simdms/growth_simulation_depth_200/Tsv"
os.mkdir(dir_path)

effect.to_csv(os.path.join(dir_path, "effects.tsv"), sep='\t', index=False)
expected.to_csv(os.path.join(dir_path, "expected.tsv"), sep='\t', index=False)
popcounts.to_csv(os.path.join(dir_path, "popcounts.tsv"), sep='\t', index=False)
seqcounts.to_csv(os.path.join(dir_path, "seqcounts.tsv"), sep='\t', index=False)
