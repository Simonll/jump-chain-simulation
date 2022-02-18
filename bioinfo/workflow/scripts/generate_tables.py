# import argparse
# import glob
# import os
# from pathlib import Path
# from typing import Any
# from typing import Dict
# from typing import List
# from typing import Optional
# from typing import Union

# import numpy as np
# import pandas as pd
# from bintools.phylobayes.mcmc_parser import posterior_M0_GTR
# from generate_figure_1 import LIST_OF_GENES

# from generate_M0GTR_pb_mpi_sup import generate_supplements_M0GTR_pb_mpi
# from generate_M0HKY_pb_modified_sup import generate_supplements_M0HKY_pb_modified
# from generate_cabc_sup import generate_supplements_CABC

# def generate_tables(
#     input_M0GTR_pb_mpi_dir: str,
#     pattern_M0GTR_pb_mpi: str,
#     metadata_M0GTR_pb_mpi: Dict[str, Union[str, float]],
#     input_M0HKY_pb_modified_dir: str,
#     pattern_M0HKY_pb_modified: str,
#     metadata_M0HKY_pb_modified: Dict[str, Union[str, float]],
#     input_M0HKY_codeml_dir: str,
#     pattern_M0HKY_codeml: str,
#     metadata_M0HKY_codeml: Dict[str, Union[str, float]],
#     input_M7HKY_codeml_dir: str,
#     pattern_M7HKY_codeml: str,
#     metadata_M7HKY_codeml: Dict[str, Union[str, float]],
#     input_M8HKY_codeml_dir: str,
#     pattern_M8HKY_codeml: str,
#     metadata_M8HKY_codeml: Dict[str, Union[str, float]],
#     output_M0GTR_pb_mpi: str,
#     output_M0HKY_pb_modified: str,
#     output_M0HKY_codeml: str,
#     output_M7HKY_codeml: str,
#     output_M8HKY_codeml: str,
# ) -> bool:

#     try:
#         generate_supplements_M0GTR_pb_mpi(
#             input_dir=input_M0GTR_pb_mpi_dir,
#             pattern=pattern_M0GTR_pb_mpi,
#             metadata=metadata_M0GTR_pb_mpi,
#             output=output_M0GTR_pb_mpi,
#         )
#     except Exception as e:
#         print("something wrong with %s, %s" % (input_M0GTR_pb_mpi_dir, str(e)))
#         return False
#     try:
#         generate_supplements_M0HKY_pb_modified(
#             input_dir=input_M0HKY_pb_modified_dir,
#             pattern=pattern_M0HKY_pb_modified,
#             metadata=metadata_M0HKY_pb_modified,
#             output=output_M0HKY_pb_modified,
#         )
#     except Exception as e:
#         print("something wrong with %s, %s" % (input_M0HKY_pb_modified_dir, str(e)))
#         return False
#     try:
#         generate_supplements_M0HKY_codeml(
#             input_dir=input_M0HKY_codeml_dir,
#             metadata=metadata_M0HKY_codeml,
#             output=output_M0HKY_codeml,
#         )
#     except Exception as e:
#         print("something wrong with %s, %s" % (input_M0HKY_codeml_dir, str(e)))
#         return False
#     try:
#         generate_supplements_M0HKY_codeml(
#             input_dir=input_M7HKY_codeml_dir,
#             metadata=metadata_M7HKY_codeml,
#             output=output_M7HKY_codeml,
#         )
#     except Exception as e:
#         print("something wrong with %s, %s" % (input_M7HKY_codeml_dir, str(e)))
#         return False
#     try:
#         generate_supplements_M7HKY_codeml(
#             input_dir=input_M8HKY_codeml_dir,
#             metadata=metadata_M8HKY_codeml,
#             output=output_M8HKY_codeml,
#         )
#     except Exception as e:
#         print("something wrong with %s, %s" % (input_M8HKY_codeml_dir, str(e)))
#         return False
#     return True
