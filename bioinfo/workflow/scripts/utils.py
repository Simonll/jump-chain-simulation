from typing import List
from typing import Tuple

import ete3 as phylo
import numpy as np


def quantile(
    arr: List[float], quantile_inf: float = 0.025, quantile_sup: float = 0.975
) -> Tuple[float, float]:
    quantile_inf_v: float = np.quantile(a=arr, q=quantile_inf)
    quantile_sup_v: float = np.quantile(a=arr, q=quantile_sup)
    return (round(quantile_inf_v, 3), round(quantile_sup_v, 3))


def tl_compute(list_of_trees) -> List[float]:
    list_of_tl: List[float] = []
    for t in list_of_trees:
        sum_ = 0
        t_ = phylo.Tree(t)
        for b in t_.traverse("postorder"):
            sum_ += b.dist
        list_of_tl += [sum_]
    return list_of_tl
