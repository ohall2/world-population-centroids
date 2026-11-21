#!/usr/bin/env python3
"""Skriver docs/extended/A3_detection.md ur A3:s CSV-filer (planen §5).

Alla tal hämtas ur docs/extended/detection_*.csv. Ändra aldrig talen för hand;
kör om tools/detection_model.py och sedan detta skript.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DS = {"census": "USA", "inegi": "Mexiko", "geostat": "Europa"}
MODEL = {"snr": "log(\\|t\\|/σ)", "abs_t": "log\\|t\\|"}


def f(x, nd=2) -> str:
    if pd.isna(x):
        return "–"
    return f"{x:.{nd}f}".replace(".", ",").replace("-", "−")


def period(ds, fr, to) -> str:
    return f"{DS.get(ds, ds)} {int(fr)}–{int(to)}"


def main(argv=None):
    p = argparse.ArgumentParser(description="Skriv A3-resultatdokumentet ur CSV.")
    p.add_argument("--dir", type=Path, default=Path("docs/extended"))
    a = p.parse_args(argv)
    d = a.dir

    rho = pd.read_csv(d / "detection_rho.csv")
    models = pd.read_csv(d / "detection_models.csv")
    test = pd.read_csv(d / "detection_test.csv")
    cal = pd.read_csv(d / "detection_calibration.csv")
    mm = pd.read_csv(d / "detection_min_move.csv")

    L = ["# A3 – Detektionsmodell (utforskande)", "",
         "Frågan (docs/ANALYSIS_PLAN.md §5): går det att ange per enhet hur stor "
         "en förflyttning måste vara för att synas?", "",
         "Brusmåttet σ per enhet tas ur utgåvan: `noise_dmsp_km` inom DMSP, "
         "√(`noise_dmsp_km`² + `noise_viirs_km`²) över skarven. e är "
         "ändpunktsskattningen (utgåvans ljuscentroider y1→y2); |t| den verkliga "
         "förflyttningen i facit. Allt på nivå 2. Måtten och facitläsningen är de "
         "frysta i `tools/change_detection.py`; logistisk regression och AUC görs "
         "i numpy. Alla tal är skrivna ur `docs/extended/detection_*.csv` med "
         "`tools/write_detection_doc.py`.", ""]

    # Steg 1
    L += ["## Steg 1: σ mot fel (Spearmans ρ, nivå 2)", "",
          "| Period | n | ρ(σ, \\|e−t\\|) |", "| --- | --- | --- |"]
    for _, r in rho.iterrows():
        L += [f"| {period(r.dataset, r['from'], r.to)} | {int(r.n)} | "
              f"{f(r.rho)} ({f(r.rho_lo)}–{f(r.rho_hi)}) |"]
    h3a = bool((rho.rho > 0).all() and (rho.rho_lo > 0).all())
    L += ["", f"**H3a** (ρ > 0 i alla facit): {'stöds' if h3a else 'stöds inte'}. "
          f"Större brus hör ihop med större fel i alla facit och perioder.", ""]

    # Steg 2
    L += ["## Steg 2: logistisk modell, anpassad på USA counties 2000–2020", "",
          "Rätt riktning ~ β0 + β1 · prediktor.", "",
          "| Modell | prediktor | n | β0 | β1 |", "| --- | --- | --- | --- | --- |"]
    for _, r in models.iterrows():
        L += [f"| {r.model} | {MODEL.get(r.model, r.predictor)} | {int(r.n)} | "
              f"{f(r.beta0, 3)} | {f(r.beta1, 3)} |"]
    L += [""]

    # Steg 3: prövning
    L += ["## Steg 3: prövning (AUC och Brier)", "",
          "| Prövningsfacit | modell | n | basandel | AUC | Brier | Brier (konstant) |",
          "| --- | --- | --- | --- | --- | --- | --- |"]
    for _, r in test.iterrows():
        L += [f"| {period(r.test_dataset, r['from'], r.to)} | {MODEL.get(r.model, r.model)} | "
              f"{int(r.n)} | {f(r.base_rate)} | {f(r.auc, 3)} | {f(r.brier, 3)} | "
              f"{f(r.brier_const, 3)} |"]
    mx = test[test.test_dataset == "inegi"].set_index("model")
    h3b = bool(mx.loc["snr", "auc"] > mx.loc["abs_t", "auc"])
    L += ["", f"**H3b** (σ-modellen högre AUC än \\|t\\|-modellen på Mexiko): "
          f"{'stöds' if h3b else 'stöds inte'}. σ-modellens AUC "
          f"{f(mx.loc['snr','auc'],3)} mot {f(mx.loc['abs_t','auc'],3)} för \\|t\\|-modellen.", ""]

    # Kalibrering (σ-modellen) på Mexiko
    mxcal = cal[cal.test_dataset == "inegi"]
    L += ["### Kalibrering på Mexiko (σ-modellen, tio lika stora klasser)", "",
          "| Klass | n | förutsagd | observerad | skillnad |", "| --- | --- | --- | --- | --- |"]
    for _, r in mxcal.iterrows():
        L += [f"| {int(r.bin)} | {int(r.n)} | {f(r.pred_mean)} | {f(r.obs_prop)} | {f(r.abs_gap)} |"]
    within = int((mxcal.abs_gap <= 0.05).sum())
    h3c = within >= 8
    L += ["", f"**H3c** (förutsagd och observerad skiljer ≤ 0,05 i minst åtta av tio "
          f"klasser på Mexiko): {'stöds' if h3c else 'stöds inte'}. {within} av tio "
          f"klasser inom 0,05.", ""]

    # Steg 4
    L += ["## Steg 4: minsta upptäckbara förflyttning (km)", "",
          "Det |t| där σ-modellen (anpassad på USA 2000–2020) ger 0,8 sannolikhet "
          "för rätt riktning, per enhet. Fördelning över enheterna i varje facit.", "",
          "| Period | n | σ median | min. flytt median | q25 | q75 |",
          "| --- | --- | --- | --- | --- | --- |"]
    for _, r in mm.iterrows():
        L += [f"| {period(r.dataset, r['from'], r.to)} | {int(r.n)} | {f(r.sigma_med)} | "
              f"{f(r.min_move_med)} | {f(r.min_move_q25)} | {f(r.min_move_q75)} |"]
    L += [""]

    # Sammanfattning
    L += ["## Hypoteser", "",
          "| Hypotes | Utfall |", "| --- | --- |",
          f"| H3a: ρ > 0 i alla facit | {'stöds' if h3a else 'stöds inte'} |",
          f"| H3b: σ-modellen högre AUC på Mexiko | {'stöds' if h3b else 'stöds inte'} |",
          f"| H3c: kalibrering ≤ 0,05 i ≥ 8/10 klasser på Mexiko | {'stöds' if h3c else 'stöds inte'} |",
          ""]

    out = d / "A3_detection.md"
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"skrev {out}  (H3a {'ja' if h3a else 'nej'}, H3b {'ja' if h3b else 'nej'}, "
          f"H3c {'ja' if h3c else 'nej'})")


if __name__ == "__main__":
    main()
