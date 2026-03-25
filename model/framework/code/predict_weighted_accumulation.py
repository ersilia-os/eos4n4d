#!/usr/bin/env python3
"""
CLI tool to predict weighted antibiotic accumulation scores using
a pre-trained AutoGluon model and TwinBooster.
"""
import argparse
import pandas as pd
import numpy as np

from chembl_structure_pipeline import standardizer
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from autogluon.tabular import TabularPredictor
import twinbooster

SMILES_CACHE = {}


def get_ecfp4_fingerprints(smiles, n_bits=1024):
    """Convert SMILES strings into ECFP4 RDKit fingerprints."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=n_bits)


def fp_to_array(fp, n_bits=1024):
    """Convert RDKit fingerprint to numpy array."""
    arr = np.zeros((n_bits,), dtype=int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def get_clean_smiles(smiles):
    """Standardize and parentize SMILES using ChEMBL structure pipeline."""
    if smiles in SMILES_CACHE:
        return SMILES_CACHE[smiles]
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            SMILES_CACHE[smiles] = None
            return None
        molblock = Chem.MolToMolBlock(mol)
        std_molblock = standardizer.standardize_molblock(molblock)
        parent_molblock, _ = standardizer.get_parent_molblock(std_molblock)
        parent_mol = Chem.MolFromMolBlock(parent_molblock)
        clean_smiles = Chem.MolToSmiles(parent_mol)
        SMILES_CACHE[smiles] = clean_smiles
        return clean_smiles
    except Exception:
        SMILES_CACHE[smiles] = None
        return None


def weighted_prediction(pred1, pred2, weight):
    """Combine two predictions with similarity-based weighting."""
    return weight * pred1 + (1 - weight) * pred2


def pairwise_similarity(fp, fp_list, top_k=5):
    """Mean Tanimoto similarity to the top-k most similar fingerprints."""
    similarities = [DataStructs.TanimotoSimilarity(fp, fp2) for fp2 in fp_list]
    return np.mean(sorted(similarities, reverse=True)[:top_k])


def predict_accumulation_autogluon(df, predictor_path, ecfp_bits):
    """Predict accumulation using a saved AutoGluon predictor."""
    predictor = TabularPredictor.load(predictor_path)

    pred_pairs = []
    for idx, s in df["smiles"].items():
        fp = get_ecfp4_fingerprints(s, n_bits=ecfp_bits)
        if fp is not None:
            pred_pairs.append((idx, fp_to_array(fp, ecfp_bits)))

    result = pd.Series(np.nan, index=df.index, dtype=float)

    if pred_pairs:
        pred_df = pd.DataFrame([x for _, x in pred_pairs])
        probs = predictor.predict_proba(pred_df)[1].values
        for (idx, _), p in zip(pred_pairs, probs):
            result.loc[idx] = p

    return result


def calculate_accumulation_twinbooster(smiles_list, tb_model, description):
    """Predict accumulation using TwinBooster."""
    predictions, confidences = tb_model.predict(
        smiles_list, description, get_confidence=True
    )
    return predictions, confidences


def calculate_pairwise_similarity(df, entry_dataset_path, ecfp_bits, top_k):
    """Compute similarity weights against the entry dataset."""
    entry = pd.read_csv(entry_dataset_path)
    entry["smiles"] = entry["smiles"].apply(get_clean_smiles)
    entry.dropna(subset=["smiles"], inplace=True)

    entry_fps = [get_ecfp4_fingerprints(s, n_bits=ecfp_bits) for s in entry["smiles"]]
    entry_fps = [fp for fp in entry_fps if fp is not None]

    result = pd.Series(np.nan, index=df.index, dtype=float)
    for idx, s in df["smiles"].items():
        fp = get_ecfp4_fingerprints(s, n_bits=ecfp_bits)
        if fp is not None:
            result.loc[idx] = pairwise_similarity(fp, entry_fps, top_k=top_k)

    return result


def parse_args():
    parser = argparse.ArgumentParser(
        description="Predict weighted antibiotic accumulation scores for a SMILES dataset."
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True,
        help="Input CSV file containing SMILES."
    )
    parser.add_argument(
        "--smiles-column", type=str, default="Smiles",
        help="Column name for SMILES in input file."
    )
    parser.add_argument(
        "-e", "--entry-dataset", type=str,
        default="data/entry_dataset/merged_cleaned_dataset.csv",
        help="Entry dataset CSV for similarity calculation."
    )
    parser.add_argument(
        "--predictor-path", type=str, required=True,
        help="Path to saved AutoGluon predictor directory."
    )
    parser.add_argument(
        "-o", "--output", type=str,
        default="ranking/weighted_accumulation_predictions.csv",
        help="Output CSV file for predictions."
    )
    parser.add_argument(
        "--ecfp-bits", type=int, default=1024,
        help="Number of bits for ECFP4 fingerprints."
    )
    parser.add_argument(
        "--top-k", type=int, default=5,
        help="Top K for similarity calculations."
    )
    parser.add_argument(
        "--tb-model-path", type=str, required=True,
        help="Path to TwinBooster model directory."
    )
    parser.add_argument(
        "--tb-lgbm-path", type=str, required=True,
        help="Path to TwinBooster LGBM model."
    )
    parser.add_argument(
        "--description", type=str,
        default="Accumulation of drugs in Gram-negative bacteria using LC–MS/MS as described in provided protocol.",
        help="Description for TwinBooster predictions."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    molecules = pd.read_csv(args.input)
    molecules.rename(columns={args.smiles_column: "smiles"}, inplace=True)
    molecules["smiles"] = molecules["smiles"].apply(get_clean_smiles)

    tb_model = twinbooster.TwinBooster(
        model_path=args.tb_model_path,
        lgbm_path=args.tb_lgbm_path
    )

    molecules["accumulation_autogluon"] = predict_accumulation_autogluon(
        molecules,
        predictor_path=args.predictor_path,
        ecfp_bits=args.ecfp_bits,
    )

    tb_preds, tb_confs = calculate_accumulation_twinbooster(
        molecules["smiles"].tolist(),
        tb_model,
        args.description,
    )
    molecules["accumulation_twinbooster"] = tb_preds
    molecules["confidence_twinbooster"] = tb_confs

    molecules["similarity_weight"] = calculate_pairwise_similarity(
        molecules,
        args.entry_dataset,
        args.ecfp_bits,
        args.top_k,
    )

    molecules["weighted_accumulation"] = weighted_prediction(
        molecules["accumulation_autogluon"],
        molecules["accumulation_twinbooster"],
        molecules["similarity_weight"],
    )

    pd.DataFrame({"weighted_accumulation": molecules["weighted_accumulation"]}).to_csv(
    args.output, index=False
    )
    print(f"Predictions saved to {args.output}")


if __name__ == "__main__":
    main()