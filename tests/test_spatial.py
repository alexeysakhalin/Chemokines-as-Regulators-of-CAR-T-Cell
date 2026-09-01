"""Tests for hierarchical and point-based spatial analyses."""

from __future__ import annotations

import unittest

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from chemokine_cart.spatial import (
    distance_to_landmarks,
    neighborhood_enrichment,
    plot_conceptual_spatial_panel,
    plot_spatial_panel,
    stratified_permutation_test,
    summarize_patient_distances,
    validate_spatial_table,
)


def make_spatial_table() -> pd.DataFrame:
    records = [
        # Patient 1, section 1: exact one-dimensional distances are easy to audit.
        ("P1", "S1", "v1", 0.0, 0.0, "endothelial", "vessel", False, False),
        ("P1", "S1", "t1", 10.0, 0.0, "tumor", "tumor_nest", False, False),
        ("P1", "S1", "s1", 20.0, 0.0, "fibroblast", "stroma", False, False),
        ("P1", "S1", "l1", 30.0, 0.0, "lymphatic", "lymphatic", False, False),
        ("P1", "S1", "q_lig", 1.0, 0.0, "myeloid", None, True, False),
        ("P1", "S1", "q_cart", 9.0, 0.0, "car_t", None, False, True),
        # Patient 1, section 2 intentionally lacks a lymphatic landmark.
        ("P1", "S2", "v1", 0.0, 5.0, "endothelial", "vessel", False, False),
        ("P1", "S2", "t1", 10.0, 5.0, "tumor", "tumor_nest", False, False),
        ("P1", "S2", "s1", 20.0, 5.0, "fibroblast", "stroma", False, False),
        ("P1", "S2", "q_lig", 2.0, 5.0, "myeloid", None, True, False),
        ("P1", "S2", "q_cart", 8.0, 5.0, "car_t", None, False, True),
    ]
    return pd.DataFrame(
        records,
        columns=[
            "patient_id",
            "section_id",
            "cell_id",
            "x",
            "y",
            "cell_type",
            "structure",
            "is_ligand_source",
            "is_car_t",
        ],
    )


class ValidateSpatialTableTests(unittest.TestCase):
    def test_validation_normalizes_coordinates_and_boolean_flags(self) -> None:
        table = make_spatial_table()
        table["x"] = table["x"].astype(str)
        table["is_car_t"] = table["is_car_t"].map({True: "yes", False: "no"})

        validated = validate_spatial_table(table, boolean_columns=("is_car_t",))

        self.assertTrue(pd.api.types.is_float_dtype(validated["x"]))
        self.assertTrue(pd.api.types.is_bool_dtype(validated["is_car_t"]))
        self.assertEqual(int(validated["is_car_t"].sum()), 2)
        self.assertEqual(table["x"].dtype, object, "Input must not be mutated.")

    def test_validation_rejects_duplicate_keys_and_nonfinite_coordinates(self) -> None:
        table = make_spatial_table()
        duplicate = pd.concat([table, table.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            validate_spatial_table(duplicate)

        nonfinite = table.copy()
        nonfinite.loc[0, "x"] = np.inf
        with self.assertRaisesRegex(ValueError, "non-finite"):
            validate_spatial_table(nonfinite)


class LandmarkDistanceTests(unittest.TestCase):
    def test_distances_are_section_specific_and_scaled(self) -> None:
        distances = distance_to_landmarks(make_spatial_table(), coordinate_scale_um=2.0)

        ligand_s1 = distances.loc[
            distances["section_id"].eq("S1")
            & distances["cell_id"].eq("q_lig")
            & distances["query_role"].eq("ligand_source")
        ].set_index("landmark")
        self.assertEqual(ligand_s1.loc["vessel", "distance_um"], 2.0)
        self.assertEqual(ligand_s1.loc["tumor_nest", "distance_um"], 18.0)
        self.assertEqual(ligand_s1.loc["stroma", "distance_um"], 38.0)
        self.assertEqual(ligand_s1.loc["lymphatic", "distance_um"], 58.0)
        self.assertEqual(ligand_s1.loc["vessel", "nearest_landmark_id"], "v1")

        car_s1 = distances.loc[
            distances["section_id"].eq("S1")
            & distances["cell_id"].eq("q_cart")
            & distances["query_role"].eq("car_t")
        ].set_index("landmark")
        self.assertEqual(car_s1.loc["tumor_nest", "distance_um"], 2.0)

        missing_local_landmark = distances.loc[
            distances["section_id"].eq("S2")
            & distances["cell_id"].eq("q_cart")
            & distances["landmark"].eq("lymphatic"),
            "distance_um",
        ]
        self.assertEqual(len(missing_local_landmark), 1)
        self.assertTrue(missing_local_landmark.isna().all())

    def test_patient_summary_weights_sections_equally(self) -> None:
        distances = pd.DataFrame(
            {
                "patient_id": ["P1", "P1", "P1"],
                "section_id": ["S1", "S1", "S2"],
                "query_role": ["car_t"] * 3,
                "landmark": ["vessel"] * 3,
                "distance_um": [1.0, 3.0, 9.0],
            }
        )

        summary = summarize_patient_distances(distances).iloc[0]

        self.assertEqual(summary["n_sections"], 2)
        self.assertEqual(summary["n_query_points"], 3)
        self.assertEqual(summary["pooled_median_distance_um"], 3.0)
        self.assertEqual(summary["mean_section_median_distance_um"], 5.5)


class StratifiedPermutationTests(unittest.TestCase):
    @staticmethod
    def make_comparison_table() -> pd.DataFrame:
        records = []
        for patient in ("P1", "P2"):
            for section in ("S1", "S2"):
                records.extend(
                    [
                        (patient, section, "modified", 1.0),
                        (patient, section, "modified", 2.0),
                        (patient, section, "control", 8.0),
                        (patient, section, "control", 9.0),
                    ]
                )
        return pd.DataFrame(
            records,
            columns=["patient_id", "section_id", "product", "distance_um"],
        )

    def test_permutation_is_hierarchical_and_reproducible(self) -> None:
        kwargs = {
            "value_col": "distance_um",
            "group_col": "product",
            "group_a": "modified",
            "group_b": "control",
            "n_permutations": 199,
            "seed": 17,
        }
        first = stratified_permutation_test(self.make_comparison_table(), **kwargs)
        second = stratified_permutation_test(self.make_comparison_table(), **kwargs)

        self.assertEqual(first["observed_difference"], -7.0)
        self.assertEqual(first["n_patients"], 2)
        self.assertEqual(first["n_sections"], 4)
        self.assertGreater(first["p_value"], 0.0)
        self.assertLessEqual(first["p_value"], 1.0)
        np.testing.assert_array_equal(first["null_distribution"], second["null_distribution"])


class NeighborhoodEnrichmentTests(unittest.TestCase):
    @staticmethod
    def clustered_table() -> pd.DataFrame:
        records = []
        for index, x in enumerate((0.0, 0.4, 0.8)):
            records.append(("P1", "S1", f"A{index}", x, 0.0, "A"))
        for index, x in enumerate((10.0, 10.4, 10.8)):
            records.append(("P1", "S1", f"B{index}", x, 0.0, "B"))
        return pd.DataFrame(
            records,
            columns=["patient_id", "section_id", "cell_id", "x", "y", "cell_type"],
        )

    def test_clustered_types_show_homotypic_enrichment(self) -> None:
        result = neighborhood_enrichment(
            self.clustered_table(), radius=1.0, n_permutations=199, seed=23
        ).set_index(["focal_type", "neighbor_type"])

        self.assertEqual(result.loc[("A", "A"), "observed_edges"], 6)
        self.assertEqual(result.loc[("A", "B"), "observed_edges"], 0)
        self.assertGreater(result.loc[("A", "A"), "log2_enrichment"], 0.0)
        self.assertLess(result.loc[("A", "B"), "log2_enrichment"], 0.0)
        self.assertTrue(result["q_value"].between(0.0, 1.0).all())

        repeated = neighborhood_enrichment(
            self.clustered_table(), radius=1.0, n_permutations=199, seed=23
        )
        pd.testing.assert_frame_equal(result.reset_index(), repeated)

    def test_neighborhood_permutation_is_invariant_to_input_row_order(self) -> None:
        original = self.clustered_table()
        shuffled = original.sample(frac=1.0, random_state=41).reset_index(drop=True)
        kwargs = {"radius": 1.0, "n_permutations": 199, "seed": 23}

        first = neighborhood_enrichment(original, **kwargs)
        second = neighborhood_enrichment(shuffled, **kwargs)

        pd.testing.assert_frame_equal(first, second)

    @staticmethod
    def patient_weighting_table(*, expanded_first_patient: bool) -> pd.DataFrame:
        records: list[tuple[str, str, str, float, float, str]] = []
        pair_count = 20 if expanded_first_patient else 1
        for pair_index in range(pair_count):
            offset = float(pair_index * 10)
            records.extend(
                [
                    ("P1", "S1", f"P1_F_{pair_index}", offset, 0.0, "F"),
                    ("P1", "S1", f"P1_N_{pair_index}", offset + 0.25, 0.0, "N"),
                ]
            )
        records.extend(
            [
                ("P2", "S1", "P2_F", 0.0, 0.0, "F"),
                ("P2", "S1", "P2_X", 0.25, 0.0, "X"),
                ("P2", "S1", "P2_N", 10.0, 0.0, "N"),
            ]
        )
        return pd.DataFrame(
            records,
            columns=["patient_id", "section_id", "cell_id", "x", "y", "cell_type"],
        )

    def test_patient_equal_effect_is_not_dominated_by_large_section(self) -> None:
        kwargs = {
            "radius": 1.0,
            "focal_types": ("F",),
            "neighbor_types": ("N",),
            "n_permutations": 49,
            "seed": 31,
        }
        compact = neighborhood_enrichment(
            self.patient_weighting_table(expanded_first_patient=False), **kwargs
        ).iloc[0]
        expanded = neighborhood_enrichment(
            self.patient_weighting_table(expanded_first_patient=True), **kwargs
        ).iloc[0]

        self.assertEqual(compact["n_patients"], 2)
        self.assertEqual(expanded["n_patients"], 2)
        self.assertEqual(compact["n_patients_with_edges"], 2)
        self.assertEqual(compact["n_patients_total"], 2)
        self.assertEqual(compact["n_sections_total"], 2)
        self.assertEqual(compact["observed_edges"], 1)
        self.assertEqual(expanded["observed_edges"], 20)
        self.assertAlmostEqual(compact["observed_patient_mean_edge_fraction"], 0.25)
        self.assertAlmostEqual(
            expanded["observed_patient_mean_edge_fraction"],
            compact["observed_patient_mean_edge_fraction"],
        )

    @staticmethod
    def patient_section_weighting_table(*, first_patient_sections: int) -> pd.DataFrame:
        records: list[tuple[str, str, str, float, float, str]] = []
        for section_index in range(first_patient_sections):
            section = f"S{section_index + 1}"
            records.extend(
                [
                    ("P1", section, f"P1_F_{section_index}", 0.0, 0.0, "F"),
                    ("P1", section, f"P1_N_{section_index}", 0.25, 0.0, "N"),
                ]
            )
        records.extend(
            [
                ("P2", "S1", "P2_F", 0.0, 0.0, "F"),
                ("P2", "S1", "P2_X", 0.25, 0.0, "X"),
                ("P2", "S1", "P2_N", 10.0, 0.0, "N"),
            ]
        )
        return pd.DataFrame(
            records,
            columns=["patient_id", "section_id", "cell_id", "x", "y", "cell_type"],
        )

    def test_patient_equal_effect_is_not_dominated_by_patient_with_many_sections(self) -> None:
        kwargs = {
            "radius": 1.0,
            "focal_types": ("F",),
            "neighbor_types": ("N",),
            "n_permutations": 49,
            "seed": 37,
        }
        compact = neighborhood_enrichment(
            self.patient_section_weighting_table(first_patient_sections=1), **kwargs
        ).iloc[0]
        many_sections = neighborhood_enrichment(
            self.patient_section_weighting_table(first_patient_sections=20), **kwargs
        ).iloc[0]

        self.assertEqual(compact["n_sections_with_edges"], 2)
        self.assertEqual(many_sections["n_sections_with_edges"], 21)
        self.assertEqual(compact["observed_edges"], 1)
        self.assertEqual(many_sections["observed_edges"], 20)
        self.assertAlmostEqual(compact["observed_patient_mean_edge_fraction"], 0.25)
        self.assertAlmostEqual(
            many_sections["observed_patient_mean_edge_fraction"],
            compact["observed_patient_mean_edge_fraction"],
        )


class PlottingTests(unittest.TestCase):
    def tearDown(self) -> None:
        plt.close("all")

    def test_real_spatial_panel(self) -> None:
        figure, axis = plot_spatial_panel(make_spatial_table(), patient_id="P1", section_id="S1")
        self.assertIs(figure, axis.figure)
        self.assertIn("Patient P1", axis.get_title())
        self.assertEqual(axis.get_aspect(), 1.0)

    def test_conceptual_panel_is_deterministic(self) -> None:
        first_figure, first_axis = plot_conceptual_spatial_panel(seed=11)
        first_offsets = first_axis.collections[0].get_offsets().copy()
        second_figure, second_axis = plot_conceptual_spatial_panel(seed=11)
        second_offsets = second_axis.collections[0].get_offsets().copy()
        self.assertIs(first_figure, first_axis.figure)
        self.assertIs(second_figure, second_axis.figure)
        np.testing.assert_array_equal(first_offsets, second_offsets)


if __name__ == "__main__":
    unittest.main()
