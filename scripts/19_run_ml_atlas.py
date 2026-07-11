#!/usr/bin/env python3
"""PoxHostAtlas interpretable machine-learning layer.

What the below code will be doing is creating a harmonized cross-study, sample-level expression matrix (gene-symbol
space) and this will be done from the three Tier-A studies that expose the raw counts (GSE278320,
GSE287860, GSE288000) and also what is going to be done is that this will ask whether host translation/helicase features are predicting being a
poxvirus infected across *unseen* studies.

The validation for this is on purpose going to be leave-one-DATASET-out (and leave-one-virus-out) and this will never
be a random sample split. This is because of random splits leak study-specific batch signal.

Outputs:
  results/ml/sample_feature_matrix.csv
  results/ml/leave_dataset_out_performance.csv
  results/ml/leave_virus_out_performance.csv
  results/ml/feature_importance_consensus.csv
  results/ml/ablation_results.csv
  results/ml/negative_control_results.csv
"""

from __future__ import annotations

import gzip
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, average_precision_score, balanced_accuracy_score,
    f1_score, matthews_corrcoef, roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[1]
PROCESSED_DIR = REPO_ROOT / "data" / "processed"
EXTERNAL_RAW = REPO_ROOT / "data" / "external" / "expansion_raw"
META_DIR = REPO_ROOT / "results" / "meta_analysis"
OUT_DIR = REPO_ROOT / "results" / "ml"

TARGET_PATTERN = re.compile(r"^(DHX|DDX|EIF|RPS|RPL)", re.IGNORECASE)
DHX_DDX = re.compile(r"^(DHX|DDX)", re.IGNORECASE)
EIF_RIBO = re.compile(r"^(EIF|RPS|RPL)", re.IGNORECASE)
ISG = ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "OAS3",
       "IFI6", "IFI27", "IFI44", "RSAD2", "STAT1", "STAT2", "IRF7", "IRF9", "DDX58", "IFIH1"]
HOUSEKEEPING = ["ACTB", "GAPDH", "TBP", "PGK1", "PPIA", "B2M", "HPRT1", "GUSB",
                "TUBB", "RPLP0", "YWHAZ", "SDHA", "UBC", "PPIB", "GAPDHS"]
RNG = np.random.default_rng(42)


def gene_symbol_map() -> dict[str, str]:
    counts = pd.read_csv(PROCESSED_DIR / "counts.csv", usecols=["gene_id", "gene_symbol"], dtype="string")
    mapping = {}
    for gid, sym in zip(counts["gene_id"], counts["gene_symbol"]):
        gid = str(gid); sym = str(sym)
        mapping[gid] = sym
        mapping[gid.split(".", 1)[0]] = sym
    return mapping


def collapse_to_symbols(counts_genes_by_samples: pd.DataFrame) -> pd.DataFrame:
    """Input rows=genes(symbol), cols=samples -> collapse duplicate symbols by sum."""
    counts_genes_by_samples.index = counts_genes_by_samples.index.astype(str).str.upper()
    return counts_genes_by_samples.groupby(level=0).sum()


def log_cpm(samples_by_genes: pd.DataFrame) -> pd.DataFrame:
    lib = samples_by_genes.sum(axis=1).replace(0, np.nan)
    cpm = samples_by_genes.div(lib, axis=0) * 1e6
    return np.log2(cpm.fillna(0) + 1.0)


def load_gse278320() -> tuple[pd.DataFrame, pd.Series]:
    counts = pd.read_csv(PROCESSED_DIR / "model_ready_counts.csv", index_col=0)  # samples x gene_id
    meta = pd.read_csv(PROCESSED_DIR / "model_ready_metadata.csv", index_col=0)
    mapping = gene_symbol_map()
    sym = pd.Series(counts.columns.map(lambda g: mapping.get(g, mapping.get(str(g).split(".")[0], g))))
    genes_by_samples = counts.T
    genes_by_samples.index = sym.values
    collapsed = collapse_to_symbols(genes_by_samples)  # symbol x samples
    sbg = collapsed.T  # samples x symbol
    labels = (meta.loc[sbg.index, "infection"].astype(str).str.upper() != "MOCK").astype(int)
    return sbg, labels


def load_raw_matrix(path: Path, control_tokens, infected_tokens, mapping: dict[str, str]) -> tuple[pd.DataFrame, pd.Series]:
    with gzip.open(path, "rt") as h:
        frame = pd.read_csv(h, sep="\t", index_col=0)
    frame.index = frame.index.astype(str).str.strip('"')
    frame.columns = frame.columns.astype(str).str.strip('"')
    # map versioned/base Ensembl IDs to HGNC symbols
    sym = frame.index.map(lambda g: mapping.get(g, mapping.get(str(g).split(".", 1)[0], g)))
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0)
    frame.index = sym
    collapsed = collapse_to_symbols(frame)
    sbg = collapsed.T
    label = {}
    for c in sbg.index:
        cl = c.lower()
        if any(t in cl for t in infected_tokens):
            label[c] = 1
        elif any(t in cl for t in control_tokens):
            label[c] = 0
    keep = [c for c in sbg.index if c in label]
    sbg = sbg.loc[keep]
    return sbg, pd.Series({c: label[c] for c in keep})


def build_atlas() -> tuple[pd.DataFrame, pd.Series, pd.Series, pd.Series]:
    mapping = gene_symbol_map()
    s1, y1 = load_gse278320()
    s2, y2 = load_raw_matrix(EXTERNAL_RAW / "GSE287860_raw_feature_counts.txt.gz",
                             ["cherry"], ["m3", "m003"], mapping)
    s3, y3 = load_raw_matrix(EXTERNAL_RAW / "GSE288000_raw_feature_counts2.txt.gz",
                             ["cherry"], ["m003", "m3"], mapping)
    common = s1.columns.intersection(s2.columns).intersection(s3.columns)
    s1, s2, s3 = s1[common], s2[common], s3[common]
    X1, X2, X3 = log_cpm(s1), log_cpm(s2), log_cpm(s3)
    X = pd.concat([X1, X2, X3])
    y = pd.concat([y1, y2, y3])
    study = pd.Series(
        ["GSE278320"] * len(X1) + ["GSE287860"] * len(X2) + ["GSE288000"] * len(X3),
        index=X.index,
    )
    virus = study.map({"GSE278320": "Vaccinia", "GSE287860": "Myxoma", "GSE288000": "Myxoma"})
    X.index = [f"{study[i]}::{i}" for i in X.index]
    y.index = X.index; study.index = X.index; virus.index = X.index
    return X, y, study, virus


def models() -> dict:
    return {
        "elasticnet_logreg": Pipeline([("sc", StandardScaler()),
            ("clf", LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5,
                                       C=0.5, max_iter=5000, class_weight="balanced"))]),
        "random_forest": RandomForestClassifier(n_estimators=400, max_depth=6,
                                                 class_weight="balanced", random_state=42),
        "gradient_boosting": GradientBoostingClassifier(n_estimators=200, max_depth=3, random_state=42),
        "linear_svm": Pipeline([("sc", StandardScaler()),
            ("clf", SVC(kernel="linear", C=0.5, probability=True, class_weight="balanced", random_state=42))]),
    }


def score(y_true, y_pred, y_prob) -> dict:
    out = {
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "f1": f1_score(y_true, y_pred, zero_division=0),
        "mcc": matthews_corrcoef(y_true, y_pred) if len(set(y_true)) > 1 else np.nan,
    }
    if len(set(y_true)) > 1:
        try:
            out["roc_auc"] = roc_auc_score(y_true, y_prob)
            out["pr_auc"] = average_precision_score(y_true, y_prob)
        except ValueError:
            out["roc_auc"] = out["pr_auc"] = np.nan
    else:
        out["roc_auc"] = out["pr_auc"] = np.nan
    return out


def feature_sets(X: pd.DataFrame) -> dict:
    cols = X.columns
    tf = [c for c in cols if TARGET_PATTERN.match(c)]
    meta_sig = []
    top = META_DIR / "top_conserved_translation_factors.csv"
    if top.exists():
        meta_sig = [g for g in pd.read_csv(top, index_col=0).index.astype(str).str.upper() if g in cols]
    fs = {
        "all_genes": list(cols),
        "translation_factors": tf,
        "dhx_ddx_helicases": [c for c in cols if DHX_DDX.match(c)],
        "eif_rps_rpl": [c for c in cols if EIF_RIBO.match(c)],
        "meta_signature": meta_sig if meta_sig else tf[:40],
    }
    n = len(fs["meta_signature"])
    fs["random_matched"] = list(RNG.choice(cols, size=min(n, len(cols)), replace=False))
    return {k: v for k, v in fs.items() if v}


def predict_holdout(model, X_tr, y_tr, X_te) -> tuple[np.ndarray, np.ndarray]:
    model.fit(X_tr, y_tr)
    pred = model.predict(X_te)
    if hasattr(model, "predict_proba"):
        prob = model.predict_proba(X_te)[:, 1]
    elif hasattr(model, "decision_function"):
        d = model.decision_function(X_te); prob = (d - d.min()) / (np.ptp(d) + 1e-9)
    else:
        prob = pred.astype(float)
    return pred, prob


def leave_group_out(X, y, group, feats, label_col, group_name) -> pd.DataFrame:
    rows = []
    for fs_name, cols in feats.items():
        Xf = X[cols]
        for m_name, model in models().items():
            for held in group.unique():
                tr = group != held; te = group == held
                if y[tr].nunique() < 2 or y[te].nunique() < 1:
                    continue
                pred, prob = predict_holdout(model, Xf[tr], y[tr], Xf[te])
                s = score(y[te].to_numpy(), pred, prob)
                s.update({"feature_set": fs_name, "model": m_name, label_col: held,
                          "n_train": int(tr.sum()), "n_test": int(te.sum()), "n_features": len(cols)})
                rows.append(s)
    df = pd.DataFrame(rows)
    fname = "leave_dataset_out_performance.csv" if group_name == "study" else "leave_virus_out_performance.csv"
    df.to_csv(OUT_DIR / fname, index=False)
    return df


def importance_consensus(X, y, study, feats) -> pd.DataFrame:
    cols = feats["translation_factors"]
    Xf = X[cols]
    agg = pd.DataFrame(index=cols)
    # elastic-net coefficient stability across leave-dataset-out folds
    coef_runs = []
    perm_runs = []
    for held in study.unique():
        tr = study != held; te = study == held
        if y[tr].nunique() < 2:
            continue
        en = Pipeline([("sc", StandardScaler()),
            ("clf", LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5,
                                       C=0.5, max_iter=5000, class_weight="balanced"))])
        en.fit(Xf[tr], y[tr])
        coef_runs.append(pd.Series(np.abs(en.named_steps["clf"].coef_[0]), index=cols))
        rf = RandomForestClassifier(n_estimators=400, max_depth=6, class_weight="balanced", random_state=42)
        rf.fit(Xf[tr], y[tr])
        if y[te].nunique() > 1:
            pi = permutation_importance(rf, Xf[te], y[te], n_repeats=10, random_state=42, scoring="balanced_accuracy")
            perm_runs.append(pd.Series(pi.importances_mean, index=cols))
    agg["mean_abs_enet_coef"] = pd.concat(coef_runs, axis=1).mean(axis=1) if coef_runs else 0.0
    agg["enet_selection_freq"] = (pd.concat(coef_runs, axis=1) > 1e-6).mean(axis=1) if coef_runs else 0.0
    agg["mean_perm_importance"] = pd.concat(perm_runs, axis=1).mean(axis=1) if perm_runs else 0.0
    for c in ["mean_abs_enet_coef", "mean_perm_importance"]:
        rng = agg[c].max() - agg[c].min()
        agg[c + "_norm"] = (agg[c] - agg[c].min()) / rng if rng > 0 else 0.0
    agg["ml_importance_score"] = (agg["mean_abs_enet_coef_norm"] + agg["mean_perm_importance_norm"]
                                  + agg["enet_selection_freq"]) / 3
    agg = agg.sort_values("ml_importance_score", ascending=False)
    agg.index.name = "gene"
    agg.to_csv(OUT_DIR / "feature_importance_consensus.csv")
    return agg


def _svm():
    return Pipeline([("sc", StandardScaler()),
        ("clf", SVC(kernel="linear", C=0.5, probability=False, class_weight="balanced", random_state=42))])


def ablation(X, y, study, feats, importance: pd.DataFrame) -> pd.DataFrame:
    """Ablation within the predictive translation-factor feature space.

    Uses the linear SVM (the model that generalizes cross-study) and asks whether
    removing the conserved meta/ML signature degrades leave-dataset-out prediction
    more than removing random matched translation factors.
    """
    base_cols = list(feats["translation_factors"])
    ranked = [g for g in importance.index.tolist() if g in base_cols]
    rows = []

    def loso_bacc(cols):
        accs = []
        for held in study.unique():
            tr = study != held; te = study == held
            if y[tr].nunique() < 2 or y[te].nunique() < 1:
                continue
            pred, _ = predict_holdout(_svm(), X[cols][tr], y[tr], X[cols][te])
            accs.append(balanced_accuracy_score(y[te], pred))
        return float(np.mean(accs)) if accs else np.nan

    baseline = loso_bacc(base_cols)
    rows.append({"condition": "baseline_translation_factors", "n_removed": 0,
                 "loso_balanced_accuracy": baseline, "drop": 0.0})

    for n in [10, 25, 50]:
        remove = set(ranked[:n])
        bacc = loso_bacc([c for c in base_cols if c not in remove])
        rows.append({"condition": f"remove_top{n}_signature", "n_removed": len(remove),
                     "loso_balanced_accuracy": bacc, "drop": baseline - bacc})
        rand_baccs = []
        for _ in range(10):
            rand_remove = set(RNG.choice(base_cols, size=len(remove), replace=False))
            rand_baccs.append(loso_bacc([c for c in base_cols if c not in rand_remove]))
        mean_rand = float(np.mean(rand_baccs))
        rows.append({"condition": f"remove_random{n}_matched", "n_removed": len(remove),
                     "loso_balanced_accuracy": mean_rand, "drop": baseline - mean_rand})

    all_cols = feats["all_genes"]
    for name, genes in [("remove_DHX_DDX", feats.get("dhx_ddx_helicases", [])),
                        ("remove_ISG", [g for g in ISG if g in all_cols])]:
        keep = [c for c in base_cols if c not in set(genes)]
        if not keep:
            continue
        bacc = loso_bacc(keep)
        rows.append({"condition": name, "n_removed": len(set(genes) & set(base_cols)),
                     "loso_balanced_accuracy": bacc, "drop": baseline - bacc})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "ablation_results.csv", index=False)
    return df


def negative_controls(X, y, study, feats) -> pd.DataFrame:
    rows = []

    def loso(model_fn, Xf, target):
        accs = []
        for held in study.unique():
            tr = study != held; te = study == held
            if target[tr].nunique() < 2 or target[te].nunique() < 1:
                continue
            pred, _ = predict_holdout(model_fn(), Xf[tr], target[tr], Xf[te])
            accs.append(balanced_accuracy_score(target[te], pred))
        return float(np.mean(accs)) if accs else np.nan

    rf = _svm
    tf_cols = feats["translation_factors"]

    rows.append({"control": "true_labels_translation_factors", "loso_balanced_accuracy": loso(rf, X[tf_cols], y)})
    # label permutation (within study to respect class balance)
    perm_accs = []
    for _ in range(10):
        yp = y.copy()
        for s in study.unique():
            idx = study[study == s].index
            yp.loc[idx] = RNG.permutation(yp.loc[idx].to_numpy())
        perm_accs.append(loso(rf, X[tf_cols], yp))
    rows.append({"control": "label_permutation_mean", "loso_balanced_accuracy": float(np.nanmean(perm_accs))})
    # housekeeping genes
    hk = [g for g in HOUSEKEEPING if g in X.columns]
    if hk:
        rows.append({"control": "housekeeping_genes", "loso_balanced_accuracy": loso(rf, X[hk], y)})
    # interferon only
    isg = [g for g in ISG if g in X.columns]
    if isg:
        rows.append({"control": "interferon_ISG_only", "loso_balanced_accuracy": loso(rf, X[isg], y)})
    # random matched
    rand_accs = [loso(rf, X[list(RNG.choice(X.columns, size=len(tf_cols), replace=False))], y) for _ in range(5)]
    rows.append({"control": "random_matched_genes_mean", "loso_balanced_accuracy": float(np.nanmean(rand_accs))})
    # study-label prediction (should be HIGH -> proves batch exists, motivates LODO)
    study_codes = study.astype("category").cat.codes
    from sklearn.model_selection import cross_val_score
    rfc = RandomForestClassifier(n_estimators=300, max_depth=6, random_state=42)
    try:
        sc = cross_val_score(rfc, X[tf_cols], study_codes, cv=3, scoring="accuracy").mean()
    except Exception:
        sc = np.nan
    rows.append({"control": "study_identity_prediction_cv3_accuracy", "loso_balanced_accuracy": sc})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "negative_control_results.csv", index=False)
    return df


def _loso_predictions(X, y, study, cols, model_fn):
    """Return out-of-fold predictions: each sample predicted when its study is held out."""
    yt, yp, ypr = [], [], []
    for held in study.unique():
        tr = study != held; te = study == held
        if y[tr].nunique() < 2 or y[te].nunique() < 1:
            continue
        m = model_fn()
        m.fit(X[cols][tr], y[tr])
        pred = m.predict(X[cols][te])
        if hasattr(m, "predict_proba"):
            prob = m.predict_proba(X[cols][te])[:, 1]
        elif hasattr(m, "decision_function"):
            d = m.decision_function(X[cols][te]); prob = (d - d.min()) / (np.ptp(d) + 1e-9)
        else:
            prob = pred.astype(float)
        yt.extend(y[te].tolist()); yp.extend(pred.tolist()); ypr.extend(np.asarray(prob).tolist())
    return np.array(yt), np.array(yp), np.array(ypr)


def _loso_bacc(X, y, study, cols, model_fn):
    accs = []
    for held in study.unique():
        tr = study != held; te = study == held
        if y[tr].nunique() < 2 or y[te].nunique() < 1:
            continue
        m = model_fn()
        m.fit(X[cols][tr], y[tr])
        accs.append(balanced_accuracy_score(y[te], m.predict(X[cols][te])))
    return float(np.mean(accs)) if accs else np.nan


def expression_variance_bins(X: pd.DataFrame, n_bins: int = 10) -> pd.Series:
    mean = X.mean(axis=0); var = X.var(axis=0)
    mb = pd.qcut(mean.rank(method="first"), n_bins, labels=False)
    vb = pd.qcut(var.rank(method="first"), n_bins, labels=False)
    return (mb.astype(str) + "_" + vb.astype(str)).rename("bin")


def null_distribution(X, y, study, feats, n_draws=1000) -> dict:
    """Expression/variance-matched random gene-set null vs the conserved signature."""
    sig = list(feats["meta_signature"])
    bins = expression_variance_bins(X)
    bin_to_genes = {b: bins.index[bins == b].tolist() for b in bins.unique()}
    sig_bacc = _loso_bacc(X, y, study, sig, _svm)
    draws = []
    for _ in range(n_draws):
        chosen = []
        for g in sig:
            pool = bin_to_genes.get(bins.get(g), X.columns.tolist())
            chosen.append(RNG.choice(pool))
        draws.append(_loso_bacc(X, y, study, list(dict.fromkeys(chosen)), _svm))
    draws = np.array([d for d in draws if np.isfinite(d)])
    pd.DataFrame({"random_set_balanced_accuracy": draws}).to_csv(
        OUT_DIR / "null_distribution_random_sets.csv", index=False)
    # empirical p: P(random >= signature)
    emp_p = (np.sum(draws >= sig_bacc) + 1) / (len(draws) + 1)
    pct = float((draws < sig_bacc).mean() * 100)
    summary = pd.DataFrame([{
        "signature_balanced_accuracy": sig_bacc,
        "n_random_sets": len(draws),
        "null_mean": float(draws.mean()), "null_sd": float(draws.std()),
        "percentile_vs_null": pct, "empirical_p": float(emp_p),
        "signature_size": len(sig),
    }])
    summary.to_csv(OUT_DIR / "signature_vs_null_summary.csv", index=False)
    print(f"Signature outperformed {pct:.1f}% of {len(draws)} matched random sets (empirical p={emp_p:.3g})")
    return {"sig_bacc": sig_bacc, "draws": draws, "percentile": pct, "emp_p": float(emp_p)}


def bootstrap_cis(X, y, study, feats, n_boot=2000) -> pd.DataFrame:
    """Bootstrap CIs for pooled out-of-fold LODO predictions (linear SVM, translation factors)."""
    cols = feats["translation_factors"]
    svm = lambda: Pipeline([("sc", StandardScaler()),
        ("clf", SVC(kernel="linear", C=0.5, probability=True, class_weight="balanced", random_state=42))])
    yt, yp, ypr = _loso_predictions(X, y, study, cols, svm)
    rows = []

    def ci(fn, name):
        vals = []
        n = len(yt)
        for _ in range(n_boot):
            idx = RNG.integers(0, n, n)
            if len(set(yt[idx])) < 2:
                continue
            try:
                vals.append(fn(yt[idx], yp[idx], ypr[idx]))
            except ValueError:
                continue
        vals = np.array(vals)
        rows.append({"metric": name, "point_estimate": fn(yt, yp, ypr),
                     "ci_lo": float(np.percentile(vals, 2.5)), "ci_hi": float(np.percentile(vals, 97.5))})

    ci(lambda t, p, pr: balanced_accuracy_score(t, p), "balanced_accuracy")
    ci(lambda t, p, pr: roc_auc_score(t, pr), "roc_auc")
    ci(lambda t, p, pr: average_precision_score(t, pr), "pr_auc")
    ci(lambda t, p, pr: f1_score(t, p, zero_division=0), "f1")
    ci(lambda t, p, pr: matthews_corrcoef(t, p), "mcc")
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "bootstrap_confidence_intervals.csv", index=False)
    print("\nBootstrap 95% CIs (LODO, linear SVM, translation factors):")
    print(df.to_string(index=False))
    return df


def per_study_performance(X, y, study, feats) -> pd.DataFrame:
    cols = feats["translation_factors"]
    svm = lambda: Pipeline([("sc", StandardScaler()),
        ("clf", SVC(kernel="linear", C=0.5, probability=True, class_weight="balanced", random_state=42))])
    rows = []
    virus_map = {"GSE278320": "Vaccinia", "GSE287860": "Myxoma", "GSE288000": "Myxoma"}
    for held in study.unique():
        tr = study != held; te = study == held
        if y[tr].nunique() < 2 or y[te].nunique() < 1:
            continue
        m = svm(); m.fit(X[cols][tr], y[tr])
        pred = m.predict(X[cols][te]); prob = m.predict_proba(X[cols][te])[:, 1]
        s = score(y[te].to_numpy(), pred, prob)
        rows.append({"held_out_study": held, "virus": virus_map.get(held, "?"),
                     "n_test": int(te.sum()), "n_infected": int(y[te].sum()),
                     "balanced_accuracy": s["balanced_accuracy"], "roc_auc": s["roc_auc"],
                     "pr_auc": s["pr_auc"], "mcc": s["mcc"]})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "per_held_out_study_performance.csv", index=False)
    print("\nPer-held-out-study performance:")
    print(df.to_string(index=False))
    return df


def ablation_significance(X, y, study, feats, importance, n_perm=1000) -> pd.DataFrame:
    """Empirical significance of signature ablation vs matched random ablation."""
    base_cols = list(feats["translation_factors"])
    ranked = [g for g in importance.index.tolist() if g in base_cols]
    baseline = _loso_bacc(X, y, study, base_cols, _svm)
    rows = []
    for n in [10, 25, 50]:
        sig_drop = baseline - _loso_bacc(X, y, study, [c for c in base_cols if c not in set(ranked[:n])], _svm)
        rand_drops = []
        for _ in range(n_perm):
            rr = set(RNG.choice(base_cols, size=n, replace=False))
            rand_drops.append(baseline - _loso_bacc(X, y, study, [c for c in base_cols if c not in rr], _svm))
        rand_drops = np.array(rand_drops)
        emp_p = (np.sum(rand_drops >= sig_drop) + 1) / (len(rand_drops) + 1)
        pct = float((rand_drops < sig_drop).mean() * 100)
        rows.append({"top_n_removed": n, "signature_drop": sig_drop,
                     "random_drop_mean": float(rand_drops.mean()),
                     "random_drop_sd": float(rand_drops.std()),
                     "signature_worse_than_pct_random": pct, "empirical_p": float(emp_p),
                     "n_permutations": len(rand_drops)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "ablation_significance.csv", index=False)
    print("\nAblation significance (signature vs matched random removal):")
    print(df.to_string(index=False))
    return df


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    X, y, study, virus = build_atlas()
    meta = pd.concat([y.rename("infected"), study.rename("study"), virus.rename("virus")], axis=1)
    meta.to_csv(OUT_DIR / "sample_metadata.csv")
    X.to_csv(OUT_DIR / "sample_feature_matrix.csv")
    print(f"Atlas: {X.shape[0]} samples x {X.shape[1]} genes; "
          f"{int(y.sum())} infected / {int((y==0).sum())} control across {study.nunique()} studies")

    feats = feature_sets(X)
    print("Feature sets:", {k: len(v) for k, v in feats.items()})

    lodo = leave_group_out(X, y, study, feats, "held_out_study", "study")
    leave_group_out(X, y, virus, feats, "held_out_virus", "virus")

    imp = importance_consensus(X, y, study, feats)
    abl = ablation(X, y, study, feats, imp)
    nc = negative_controls(X, y, study, feats)
    null_distribution(X, y, study, feats, n_draws=1000)
    bootstrap_cis(X, y, study, feats, n_boot=2000)
    per_study_performance(X, y, study, feats)
    ablation_significance(X, y, study, feats, imp, n_perm=1000)

    best = lodo[lodo["feature_set"] == "translation_factors"].groupby("model")["balanced_accuracy"].mean().sort_values(ascending=False)
    print("\nLeave-dataset-out balanced accuracy (translation factors):")
    print(best.to_string())
    print("\nAblation (drop in LOSO balanced accuracy):")
    print(abl[["condition", "loso_balanced_accuracy", "drop"]].to_string(index=False))
    print("\nNegative controls:")
    print(nc.to_string(index=False))
    print("\nTop ML-important translation factors:")
    print(imp.head(10)[["ml_importance_score"]].to_string())


if __name__ == "__main__":
    main()
