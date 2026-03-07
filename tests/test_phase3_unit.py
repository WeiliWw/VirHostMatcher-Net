import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from predictor import HostPredictor
from src.Variables import REGRESSION_COEFFICIENTS
from src.config import RuntimeConfig
from src.neighborhood import neighborhood_calculator
from src.pipeline import validate_runtime_config


class NeighborhoodTests(unittest.TestCase):
    def test_neighborhood_outputs_expected_shape_and_values(self):
        query_virus = pd.DataFrame(
            [[0.2, 0.8]],
            index=["q1"],
            columns=["b1", "b2"],
        )
        interaction = pd.DataFrame(
            [[1, 0], [0, 1]],
            index=["b1", "b2"],
            columns=["h1", "h2"],
        )

        pos_sv, neg_sv = neighborhood_calculator(query_virus, interaction)
        self.assertEqual(pos_sv.shape, (1, 2))
        self.assertEqual(neg_sv.shape, (1, 2))
        self.assertAlmostEqual(pos_sv.loc["q1", "h1"], 0.2)
        self.assertAlmostEqual(pos_sv.loc["q1", "h2"], 0.8)
        self.assertAlmostEqual(neg_sv.loc["q1", "h1"], 0.8)
        self.assertAlmostEqual(neg_sv.loc["q1", "h2"], 0.2)


class ConfigValidationTests(unittest.TestCase):
    def _make_config(self, root, data_dir):
        return RuntimeConfig(
            query_virus_dir=os.path.join(root, "query"),
            output_dir=os.path.join(root, "out"),
            intermediate_dir=os.path.join(root, "tmp"),
            data_dir=data_dir,
            genome_list=None,
            short_contig=False,
            topN=1,
            num_threads=1,
        )

    def test_validate_runtime_config_accepts_valid_layout(self):
        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "query"))
            data_dir = os.path.join(root, "data")
            os.makedirs(os.path.join(data_dir, "tables"))
            os.makedirs(os.path.join(data_dir, "crispr_db_prefix"))
            os.makedirs(os.path.join(data_dir, "host_wish_model"))

            config = self._make_config(root, data_dir)
            validate_runtime_config(config)

    def test_validate_runtime_config_rejects_missing_required_subdir(self):
        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "query"))
            data_dir = os.path.join(root, "data")
            os.makedirs(os.path.join(data_dir, "tables"))
            os.makedirs(os.path.join(data_dir, "crispr_db_prefix"))

            config = self._make_config(root, data_dir)
            with self.assertRaises(SystemExit):
                validate_runtime_config(config)


class PredictorScoringTests(unittest.TestCase):
    def test_get_scores_matches_manual_logistic_transform(self):
        predictor = HostPredictor.__new__(HostPredictor)
        predictor._short = False
        predictor._virus_index = pd.Index(["q1"])
        predictor._host_index = pd.Index(["h1", "h2"])
        predictor.s2star = pd.DataFrame([[0.2, 0.3]], index=["q1"], columns=["h1", "h2"])
        predictor.posSV = pd.DataFrame([[0.4, 0.1]], index=["q1"], columns=["h1", "h2"])
        predictor.negSV = pd.DataFrame([[0.05, 0.2]], index=["q1"], columns=["h1", "h2"])
        predictor.crispr = pd.DataFrame([[0.0, 0.6]], index=["q1"], columns=["h1", "h2"])

        predictor.get_scores()

        intercept, s2c, posc, negc, crc = REGRESSION_COEFFICIENTS
        linear = intercept + s2c * predictor.s2star + posc * predictor.posSV + negc * predictor.negSV + crc * predictor.crispr
        expected = 1 - 1 / (np.exp(linear) + 1)

        np.testing.assert_allclose(predictor.score.values, expected.values, rtol=1e-10, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
