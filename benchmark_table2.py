import timeit

setup_code = """
from src.amplifyp.dna import Nucleotides
from src.amplifyp.settings import DEFAULT_BASE_PAIR_WEIGHTS

class BasePairWeightsTbl_OLD:
    def __init__(self, row, col, weight):
        self.__row = row
        self.__col = col
        self.__weight = {}
        self.__matrix = [[0] * len(col) for _ in range(len(row))]
        self.__row_max = {}
        for i, row_val in enumerate(self.__row):
            if row_val != Nucleotides.GAP:
                self.__row_max[row_val] = max(weight[i])
            for j, col_val in enumerate(self.__col):
                if Nucleotides.GAP in [row_val, col_val]:
                    val = 0
                else:
                    val = weight[i][j]
                self.__weight[row_val, col_val] = val
                self.__matrix[i][j] = val

class BasePairWeightsTbl_NEW_OR:
    def __init__(self, row, col, weight):
        self.__row = row
        self.__col = col
        self.__weight = {}
        self.__matrix = [[0] * len(col) for _ in range(len(row))]
        self.__row_max = {}
        for i, row_val in enumerate(self.__row):
            if row_val != Nucleotides.GAP:
                self.__row_max[row_val] = max(weight[i])
            for j, col_val in enumerate(self.__col):
                if row_val == Nucleotides.GAP or col_val == Nucleotides.GAP:
                    val = 0
                else:
                    val = weight[i][j]
                self.__weight[row_val, col_val] = val
                self.__matrix[i][j] = val

class BasePairWeightsTbl_NEW_VALUE:
    def __init__(self, row, col, weight):
        self.__row = row
        self.__col = col
        self.__weight = {}
        self.__matrix = [[0] * len(col) for _ in range(len(row))]
        self.__row_max = {}
        gap = Nucleotides.GAP.value
        for i, row_val in enumerate(self.__row):
            if row_val != gap:
                self.__row_max[row_val] = max(weight[i])
            for j, col_val in enumerate(self.__col):
                if row_val == gap or col_val == gap:
                    val = 0
                else:
                    val = weight[i][j]
                self.__weight[row_val, col_val] = val
                self.__matrix[i][j] = val

row = Nucleotides.PRIMER
col = Nucleotides.TEMPLATE
weight = DEFAULT_BASE_PAIR_WEIGHTS._BasePairWeightsTbl__matrix
"""

print("Old INIT:", timeit.timeit("BasePairWeightsTbl_OLD(row, col, weight)", setup=setup_code, number=100000))
print("New INIT OR:", timeit.timeit("BasePairWeightsTbl_NEW_OR(row, col, weight)", setup=setup_code, number=100000))
print("New INIT VALUE:", timeit.timeit("BasePairWeightsTbl_NEW_VALUE(row, col, weight)", setup=setup_code, number=100000))
