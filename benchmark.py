import timeit

setup_code = """
from src.amplifyp.dna import Nucleotides
row_val = 'A'
col_val = 'T'
"""

stmt_list = "Nucleotides.GAP in [row_val, col_val]"
stmt_tuple = "Nucleotides.GAP in (row_val, col_val)"
stmt_or = "Nucleotides.GAP == row_val or Nucleotides.GAP == col_val"
stmt_or_rev = "row_val == Nucleotides.GAP or col_val == Nucleotides.GAP"

print("list:", timeit.timeit(stmt_list, setup=setup_code, number=10000000))
print("tuple:", timeit.timeit(stmt_tuple, setup=setup_code, number=10000000))
print("or:", timeit.timeit(stmt_or, setup=setup_code, number=10000000))
print("or_rev:", timeit.timeit(stmt_or_rev, setup=setup_code, number=10000000))
