"""Deterministic, no-GE-solve checks for source-equation audit."""
import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np


S = Path(__file__).resolve().parents[1] / "scripts" / "stage6_author_depreciation_audit.py"
spec = importlib.util.spec_from_file_location("stage6_audit", S)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class Parameter:
    k_initial = 10
    q_initial = 1
    k_const_hat = 0.25


class TestStage6Audit(unittest.TestCase):
    def test_gross_portfolio_and_revaluation_identity(self):
        grid = np.array([-1., 1., 3.])
        mu = np.array([[0.2], [0.4], [0.4]])
        k, b = audit.source_portfolio(grid, mu, Parameter())
        self.assertTrue(np.all(k >= 0))
        self.assertTrue(np.allclose(k + b, grid))
        ratio = 1.2
        wealth = k + ratio * b
        self.assertTrue(np.allclose(wealth, ratio*grid+(1-ratio)*k))

    def test_source_scalar_roots_and_nonzero_minimizer(self):
        # C + E/x - D/x^3 = 0 admits two positive roots for C<0.
        c, e, d = -1., 3., 1.
        f = lambda x: c + e/x - d/x**3
        roots = audit.positive_roots(f)
        self.assertEqual(len(roots), 2)
        for root in roots:
            self.assertLess(abs(f(root)), 1e-10)
        self.assertGreater(abs(f(.8)), 1e-4)

    def test_budget_exact_and_negative_resources(self):
        # The source consumption floor cannot satisfy the borrowing constraint
        # if resources R a + income - bmin <= 0.
        R, a, labor, bmin, m = 1.01, -.5, .1, -.1, 1.5
        resources = R*a + labor - bmin
        self.assertLess(resources, 0)
        c = 1e-5
        next_asset = bmin
        residual = m*c + next_asset - R*a - labor
        self.assertGreater(residual, 0)

    def test_exact_source_goods_equation(self):
        f = lambda x: audit.scalar_goods(
            x, a_pre=10., a_target=9.5, wage_labor=2., profits=.1,
            gross_return=1.02, x_pre=1., B_pre=-2., y=4.,
            investment=.3, adj_cost=.1, d_export=.7,
            theta_export=3., omega=.4, theta=1.)
        x = 1.
        expected = ((10-2)*1.02-9.5+2+.1)/(1+.4/.6)-(4-.3-.1-.7)
        self.assertAlmostEqual(f(x), expected)

    def test_recover_saved_original_price(self):
        root = Path(__file__).resolve().parents[1] / "results" / "author_reference"
        prices = np.loadtxt(root/"leveraged_source_terms.csv")
        meta = np.loadtxt(root/"shock_inherited_metadata.csv")
        self.assertEqual(prices.shape, (200,))
        self.assertAlmostEqual(prices[0], .90985851141107077)
        self.assertAlmostEqual(100*math.log(meta[3]/prices[0]),
                               11.99871760, places=6)


if __name__ == "__main__":
    unittest.main()
