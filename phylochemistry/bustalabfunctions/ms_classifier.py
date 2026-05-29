#!/usr/bin/env python3
"""
ms_classifier.py — GC-MS mass spectral classifier for wax compound annotation.

Single-quad EI-GCMS (70eV). Trained on m/z abundance vectors (40–1000, normalized 0–100).
27 compound classes: alkanes, fatty acids, primary alcohols, triterpenoids, ketones, etc.

Usage:
    python3 ms_classifier.py train
    python3 ms_classifier.py predict ms_wide_for_prediction.csv predictions.csv
"""

import sys
import os
import pickle
import argparse

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
from sklearn.feature_selection import f_classif
from sklearn.metrics import accuracy_score, classification_report

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TRAIN_CSV = os.path.join(LIB_DIR, "data_for_forest_building.csv")
NOISY_CSVS = [
    os.path.join(LIB_DIR, f"noisy_data_for_forest_building_{i}.csv")
    for i in range(1, 11)
]
MODEL_PKL = os.path.join(LIB_DIR, "busta_lab_xgb_model_v1.pkl")

N_TOP_FSCORE = 200       # top m/z features by ANOVA F-score
N_TOP_ABUNDANCE = 150    # top m/z features by summed abundance
CONFIDENCE_THRESHOLD = 0.40
N_CV_FOLDS = 5
RANDOM_STATE = 42

# Canonical EI neutral losses for wax compounds (fatty acids, alkanes, alcohols,
# triterpenoids). Excludes 1, 2, 13 (isotope artifacts).
CANONICAL_NEUTRAL_LOSSES = [
    14, 15, 16, 18, 28, 42, 44, 56, 57, 58,
    72, 73, 74, 88, 102, 134, 140, 146, 147,
    148, 150, 154, 155, 162, 164, 168, 170,
    182, 183, 184, 185, 196, 197, 198,
]


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_training_data():
    """Load original and all available noisy training CSVs.

    Returns
    -------
    X_orig : ndarray, shape (415, n_mz)
    y_str  : ndarray of str, shape (415,)
    mz_ints : list of int, length n_mz  — the m/z integer values in column order
    noisy_Xs : list of ndarray, same shape as X_orig
    """
    df = pd.read_csv(TRAIN_CSV)
    y_str = df["compound_common_name"].values
    mz_cols = [c for c in df.columns if c not in ("compound_common_name", "FileName")]
    X_orig = df[mz_cols].values.astype(np.float64)
    mz_ints = [int(c) for c in mz_cols]  # [40, 41, ..., 1000]

    noisy_Xs = []
    for path in NOISY_CSVS:
        if os.path.exists(path):
            nd = pd.read_csv(path)
            noisy_Xs.append(nd[mz_cols].values.astype(np.float64))
        else:
            print(f"  Warning: noisy file not found, skipping: {path}")

    return X_orig, y_str, mz_ints, noisy_Xs


# ---------------------------------------------------------------------------
# Feature configuration (computed once on training data, frozen for prediction)
# ---------------------------------------------------------------------------

def compute_feature_config(X_orig, y_str, mz_ints):
    """Select features from training data.

    Uses union of:
      - top N_TOP_FSCORE m/z by ANOVA F-score (log1p-transformed) — captures
        discriminative molecular ion / high-mass fragment region
      - top N_TOP_ABUNDANCE m/z by summed abundance — captures high-signal
        low-mass fragment series (common hydrocarbon ions at 57, 71, 85, etc.)

    Also stores neutral loss list and confidence threshold.

    Returns a config dict that is saved inside the model pkl and used identically
    at train time and predict time.
    """
    le = LabelEncoder()
    y_enc = le.fit_transform(y_str)

    col_sums = X_orig.sum(axis=0)  # shape (n_mz,)

    # F-score on log1p-transformed non-zero columns only
    nonzero_mask = col_sums > 0
    X_nz = X_orig[:, nonzero_mask]
    mz_nz = [mz_ints[i] for i, keep in enumerate(nonzero_mask) if keep]
    X_log = np.log1p(X_nz)
    f_scores, _ = f_classif(X_log, y_enc)
    f_scores = np.nan_to_num(f_scores, nan=0.0)
    top_f_idx = np.argsort(f_scores)[::-1][:N_TOP_FSCORE]
    top_f_mzs = {mz_nz[i] for i in top_f_idx}

    # Top by summed abundance (all columns including zeros)
    top_a_idx = np.argsort(col_sums)[::-1][:N_TOP_ABUNDANCE]
    top_a_mzs = {mz_ints[i] for i in top_a_idx}

    selected_mzs = sorted(top_f_mzs | top_a_mzs)

    return {
        "selected_mzs": selected_mzs,
        "neutral_losses": CANONICAL_NEUTRAL_LOSSES,
        "confidence_threshold": CONFIDENCE_THRESHOLD,
        "class_names": list(le.classes_),
        "all_train_mzs": mz_ints,  # for reference; not used at predict time
        "n_features": len(selected_mzs) + len(CANONICAL_NEUTRAL_LOSSES),
    }


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

def extract_features(X_raw, mz_index_map, config):
    """Build feature matrix from raw m/z abundance vectors.

    Works identically for training data (m/z 40–1000) and prediction data
    from R (m/z 0–1000). Column alignment is by m/z value, not position.

    Parameters
    ----------
    X_raw : ndarray, shape (n_samples, n_cols)
    mz_index_map : dict  {mz_int: column_index_in_X_raw}
    config : dict from compute_feature_config()

    Returns
    -------
    X_feat : ndarray, shape (n_samples, n_features)
        [log1p(selected m/z abundances)] + [neutral loss features]
    """
    selected_mzs = config["selected_mzs"]
    neutral_losses = config["neutral_losses"]
    n = X_raw.shape[0]

    # Part 1: log1p of selected m/z abundances
    mz_feat = np.zeros((n, len(selected_mzs)), dtype=np.float64)
    for j, mz in enumerate(selected_mzs):
        if mz in mz_index_map:
            mz_feat[:, j] = X_raw[:, mz_index_map[mz]]
    mz_feat_log = np.log1p(mz_feat)

    # Part 2: neutral loss features
    # For each spectrum: find base peak (highest m/z with >5% relative abundance),
    # then record log1p abundance at base_mz - loss for each canonical loss.
    nl_feat = np.zeros((n, len(neutral_losses)), dtype=np.float64)
    all_mzs_sorted = sorted(mz_index_map.keys())  # sorted list of available m/z

    for i in range(n):
        row = X_raw[i, :]
        max_val = row.max()
        if max_val == 0:
            continue  # all-zero spectrum — leave neutral loss features as zero

        rel_abund = (row / max_val) * 100.0

        # Significant ions: >5% relative abundance
        sig_mzs = [mz for mz in all_mzs_sorted if rel_abund[mz_index_map[mz]] > 5.0]
        if not sig_mzs:
            continue

        # Base peak: highest m/z among significant ions (favors molecular ion region)
        base_mz = max(sig_mzs)

        for k, loss in enumerate(neutral_losses):
            target_mz = base_mz - loss
            if target_mz in mz_index_map:
                nl_feat[i, k] = np.log1p(row[mz_index_map[target_mz]])

    return np.hstack([mz_feat_log, nl_feat])


# ---------------------------------------------------------------------------
# Classifier factory
# ---------------------------------------------------------------------------

def get_classifier():
    """Return XGBoost classifier if available, else HistGradientBoosting."""
    try:
        import xgboost as xgb
        clf = xgb.XGBClassifier(
            n_estimators=300,
            max_depth=4,
            learning_rate=0.05,
            subsample=0.8,
            colsample_bytree=0.8,
            min_child_weight=1,
            gamma=0,
            reg_alpha=0.1,
            reg_lambda=1.0,
            eval_metric="mlogloss",
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )
        print("  Using XGBoost classifier")
        return clf
    except ImportError:
        from sklearn.ensemble import HistGradientBoostingClassifier
        clf = HistGradientBoostingClassifier(
            max_iter=300,
            max_depth=4,
            learning_rate=0.05,
            min_samples_leaf=5,
            random_state=RANDOM_STATE,
        )
        print("  XGBoost not found — using sklearn HistGradientBoostingClassifier")
        return clf


# ---------------------------------------------------------------------------
# Cross-validation (leak-free)
# ---------------------------------------------------------------------------

def run_cv(X_orig, y_str, mz_ints, noisy_Xs, config):
    """5-fold stratified CV with augmentation applied to training folds only.

    Validation folds use original (unaugmented) spectra only, so accuracy
    estimates are honest and not inflated by noise copies of training data.
    """
    le = LabelEncoder()
    le.classes_ = np.array(config["class_names"])
    y = le.transform(y_str)

    mz_index_map = {mz: i for i, mz in enumerate(mz_ints)}

    skf = StratifiedKFold(n_splits=N_CV_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    fold_accuracies = []
    all_val_true, all_val_pred = [], []

    for fold, (train_idx, val_idx) in enumerate(skf.split(X_orig, y)):
        # Validation: original spectra only
        X_val_raw = X_orig[val_idx]
        y_val = y[val_idx]

        # Training: original fold + all noisy copies at same row indices
        X_train_parts = [X_orig[train_idx]]
        y_train_parts = [y[train_idx]]
        for noisy_X in noisy_Xs:
            X_train_parts.append(noisy_X[train_idx])
            y_train_parts.append(y[train_idx])

        X_train_raw = np.vstack(X_train_parts)
        y_train = np.concatenate(y_train_parts)

        X_train_feat = extract_features(X_train_raw, mz_index_map, config)
        X_val_feat = extract_features(X_val_raw, mz_index_map, config)

        clf = get_classifier()
        clf.fit(X_train_feat, y_train)

        y_pred = clf.predict(X_val_feat)
        acc = accuracy_score(y_val, y_pred)
        fold_accuracies.append(acc)
        all_val_true.extend(y_val)
        all_val_pred.extend(y_pred)
        print(f"  Fold {fold + 1}/{N_CV_FOLDS}: accuracy = {acc:.4f}")

    overall_acc = accuracy_score(all_val_true, all_val_pred)
    print(f"\nCV Accuracy (mean ± std): {np.mean(fold_accuracies):.4f} ± {np.std(fold_accuracies):.4f}")
    print(f"Pooled CV Accuracy:       {overall_acc:.4f}")
    print("\nPer-class report (pooled CV predictions on original spectra only):")
    print(
        classification_report(
            all_val_true,
            all_val_pred,
            target_names=config["class_names"],
            zero_division=0,
        )
    )
    return fold_accuracies


# ---------------------------------------------------------------------------
# Train command
# ---------------------------------------------------------------------------

def cmd_train():
    print("=" * 60)
    print("Loading training data...")
    X_orig, y_str, mz_ints, noisy_Xs = load_training_data()
    print(f"  {X_orig.shape[0]} original spectra, {len(np.unique(y_str))} classes")
    print(f"  {len(noisy_Xs)} noisy augmentation files loaded")

    print("\nComputing feature configuration...")
    config = compute_feature_config(X_orig, y_str, mz_ints)
    print(f"  Selected m/z features:   {len(config['selected_mzs'])}")
    print(f"  Neutral loss features:   {len(config['neutral_losses'])}")
    print(f"  Total feature dimensions: {config['n_features']}")

    print(f"\nRunning {N_CV_FOLDS}-fold stratified CV")
    print("  (augmentation applied to training folds only — honest accuracy estimate)")
    run_cv(X_orig, y_str, mz_ints, noisy_Xs, config)

    print("\nTraining final model on ALL data (original + all noisy copies)...")
    le = LabelEncoder()
    y = le.fit_transform(y_str)
    mz_index_map = {mz: i for i, mz in enumerate(mz_ints)}

    X_all_parts = [X_orig] + noisy_Xs
    y_all = np.concatenate([y] + [y] * len(noisy_Xs))
    X_all_raw = np.vstack(X_all_parts)

    print(f"  Total training rows: {X_all_raw.shape[0]}")
    X_all_feat = extract_features(X_all_raw, mz_index_map, config)

    clf = get_classifier()
    clf.fit(X_all_feat, y_all)

    bundle = {"model": clf, "config": config}
    with open(MODEL_PKL, "wb") as f:
        pickle.dump(bundle, f)
    print(f"\nModel saved to: {MODEL_PKL}")

    # Confidence threshold analysis on original training spectra (in-sample)
    print("\n--- Confidence threshold analysis (in-sample, original spectra) ---")
    X_orig_feat = extract_features(X_orig, mz_index_map, config)
    proba = clf.predict_proba(X_orig_feat)
    max_proba = proba.max(axis=1)
    for thresh in [0.30, 0.40, 0.50, 0.60, 0.70, 0.80]:
        n_unknown = (max_proba < thresh).sum()
        pct = n_unknown / len(max_proba) * 100
        print(f"  threshold={thresh:.2f}: {n_unknown:3d}/{len(max_proba)} flagged as 'unknown' ({pct:.1f}%)")
    print("=" * 60)


# ---------------------------------------------------------------------------
# Predict command
# ---------------------------------------------------------------------------

def cmd_predict(input_csv, output_csv):
    # Load model
    if not os.path.exists(MODEL_PKL):
        print(f"ERROR: Model not found at {MODEL_PKL}")
        print("Run 'python3 ms_classifier.py train' first.")
        sys.exit(1)

    with open(MODEL_PKL, "rb") as f:
        bundle = pickle.load(f)
    clf = bundle["model"]
    config = bundle["config"]

    # Read input from R (columns: peak_unique_id, 0, 1, ..., 1000)
    df = pd.read_csv(input_csv)

    # Separate peak_unique_id from m/z columns
    # m/z columns are those whose names are integers (possibly as strings)
    mz_col_names = []
    non_mz_cols = []
    for c in df.columns:
        try:
            int(c)
            mz_col_names.append(c)
        except ValueError:
            non_mz_cols.append(c)

    mz_ints_in_file = [int(c) for c in mz_col_names]
    X_raw = df[mz_col_names].values.astype(np.float64)
    mz_index_map = {mz: i for i, mz in enumerate(mz_ints_in_file)}

    X_feat = extract_features(X_raw, mz_index_map, config)

    proba = clf.predict_proba(X_feat)
    class_names = config["class_names"]
    threshold = config["confidence_threshold"]

    predicted = []
    confidence = []
    for i in range(len(X_feat)):
        max_idx = int(np.argmax(proba[i]))
        max_prob = float(proba[i, max_idx])
        if max_prob >= threshold:
            predicted.append(class_names[max_idx])
        else:
            predicted.append("unknown")
        confidence.append(round(max_prob, 4))

    out_df = pd.DataFrame({"predicted_compound": predicted, "confidence": confidence})
    out_df.to_csv(output_csv, index=False)

    n_unknown = sum(1 for p in predicted if p == "unknown")
    print(f"Predictions written to {output_csv}")
    print(f"  {len(predicted)} peaks processed")
    print(f"  {n_unknown} peaks below confidence threshold ({threshold}) → 'unknown'")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="GC-MS mass spectral classifier for wax compound annotation"
    )
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("train", help="Train model on data_for_forest_building.csv and save pkl")

    predict_p = sub.add_parser("predict", help="Predict compound identities from ms_wide CSV")
    predict_p.add_argument("input_csv", help="Path to ms_wide_for_prediction.csv (from R app)")
    predict_p.add_argument("output_csv", help="Path to write predictions.csv")

    args = parser.parse_args()

    if args.command == "train":
        cmd_train()
    elif args.command == "predict":
        cmd_predict(args.input_csv, args.output_csv)
    else:
        parser.print_help()
        sys.exit(1)
