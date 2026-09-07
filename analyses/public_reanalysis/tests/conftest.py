from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

MODULE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_ROOT))


@pytest.fixture
def synthetic_cells() -> pd.DataFrame:
    rows = []
    values = {
        "P1": {"IP": [0, 0, 1, 0], "D7": [1, 1, 0, 1]},
        "P2": {"IP": [0, 1, 0, 0], "D7": [1, 0, 0, 0]},
        "P3": {"IP": [1, 0, 0, 0], "D7": [0, 0, 0, 0]},
    }
    for patient, stages in values.items():
        for stage, counts in stages.items():
            for index, count in enumerate(counts):
                rows.append(
                    {
                        "dataset": "TEST",
                        "patient_id": patient,
                        "product": "product-a",
                        "stage": stage,
                        "biological_sample_id": f"{patient}_{stage}",
                        "cell_id": f"cell{index}",
                        "total_umi": 100 + count,
                        "cxcr6_umi": count,
                    }
                )
    return pd.DataFrame(rows)


@pytest.fixture
def synthetic_contrasts() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "dataset": "TEST",
                "contrast": "D7_minus_IP",
                "baseline_stage": "IP",
                "followup_stage": "D7",
                "min_cells": 1,
                "role": "test",
                "family_id": "none",
            }
        ]
    )
