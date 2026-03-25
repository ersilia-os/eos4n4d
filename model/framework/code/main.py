# imports
import os
import sys
import numpy as np
from ersilia_pack_utils.core import read_smiles, write_out

from predict_weighted_accumulation import (
    get_clean_smiles,
    predict_accumulation_autogluon,
    calculate_accumulation_twinbooster,
    calculate_pairwise_similarity,
    weighted_prediction,
)
import pandas as pd
import twinbooster

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# checkpoints
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))
predictor_path = os.path.join(checkpoints_dir, "autogluon_predictor")
entry_dataset_path = os.path.join(checkpoints_dir, "data", "merged_cleaned_dataset.csv")
tb_model_path = os.path.join(checkpoints_dir, "twinbooster", "scripts", "barlow_twins", "best_model")
tb_lgbm_path = os.path.join(checkpoints_dir, "twinbooster", "scripts", "lgbm", "best_model", "lgbm_model.joblib")

# constants
ECFP_BITS = 1024
TOP_K = 5
DESCRIPTION = (
    "Accumulation of drugs in Gram-negative bacteria using LC–MS/MS as described in provided protocol."
)


def predict_weighted_accumulation(smiles_list):
    molecules = pd.DataFrame({"smiles": smiles_list})
    molecules["smiles"] = molecules["smiles"].apply(get_clean_smiles)

    tb_model = twinbooster.TwinBooster(
        model_path=tb_model_path,
        lgbm_path=tb_lgbm_path,
    )

    accumulation_autogluon = predict_accumulation_autogluon(
        molecules,
        predictor_path=predictor_path,
        ecfp_bits=ECFP_BITS,
    )

    accumulation_twinbooster, _ = calculate_accumulation_twinbooster(
        molecules["smiles"].tolist(),
        tb_model,
        DESCRIPTION,
    )

    similarity_weight = calculate_pairwise_similarity(
        molecules,
        entry_dataset_path,
        ECFP_BITS,
        TOP_K,
    )

    weighted_accumulation = weighted_prediction(
        accumulation_autogluon,
        accumulation_twinbooster,
        similarity_weight,
    )

    return np.array(weighted_accumulation, dtype=np.float32).reshape(-1, 1)


# read SMILES from .csv file, assuming one column with header
_, smiles_list = read_smiles(input_file)

# run model
outputs = predict_weighted_accumulation(smiles_list)

#check input and output have the same length
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

header = ["weighted_accumulation"]

# write output in a .csv file
write_out(outputs, header, output_file, np.float32)
