"""Stage 7 source-constrained-policy budget and interpolation unit tests.

Full saved 200-quarter original MATLAB workspace is evaluated in CI script.
"""
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from stage7_constrained_policy_audit import (
    source_constrained_coarse,
    budget_consistent_binding,
    audit_stage7,
)
from stage6_author_depreciation_audit import read_original


class TestStage7ConstrainedPolicy(unittest.TestCase):
    def test_unfloored_affine_source_pchip_is_exact(self):
        coarse=np.array([-3.0,-2.0,-1.0,0.0,1.0,3.0])
        fine=np.linspace(-3,3,161)
        raw,_=source_constrained_coarse(coarse,np.array([5.0]),1.02,-.1,1.6)
        interp=PchipInterpolator(coarse,raw[:,0])(fine)
        exact=(1.02*fine+5+.1)/1.6
        self.assertLess(np.max(np.abs(interp-exact)),1e-12)

    def test_floor_kink_breaks_linear_budget_between_coarse_knots(self):
        coarse=np.array([-3.,-2.,-1.,0.,1.,3.])
        fine=np.linspace(-3,3,161)
        raw,capped=source_constrained_coarse(coarse,np.array([1.5]),1.02,-.1,1.6)
        source=PchipInterpolator(coarse,capped[:,0])(fine)
        exact=(1.02*fine+1.5+.1)/1.6
        self.assertTrue(np.any(raw<1e-5))
        self.assertGreater(np.max(np.abs(source-exact)),1e-3)
        # At an occupied positive-resource state near the coarse kink,
        # a nonzero interpolation error does not imply negative consumption.
        near=(fine>=-1.3)&(fine<=-.9)
        self.assertTrue(np.any(source[near]>0))

    def test_occupied_binding_exact_budget_and_unchanged_nonbinding(self):
        a=np.array([-.5,0.,1.])
        mu=np.array([[.2],[.3],[.5]])
        income=np.array([1.])
        r=1.01; bmin=-.1; m=1.5
        s=np.array([[-.1],[-.1],[.8]])
        old=np.array([[.2],[.7],[.8]])
        corrected,binding,exact=budget_consistent_binding(
            old,s,a,income,r,bmin,m,mu)
        self.assertTrue(np.array_equal(binding[:,0],[True,True,False]))
        np.testing.assert_allclose(
            m*corrected[binding]+s[binding],
            (r*a[:,None]+income[None,:])[binding],atol=1e-14)
        self.assertEqual(corrected[2,0],old[2,0])
        self.assertTrue(np.all(corrected[binding]>0))

    def test_real_author_saved_regression(self):
        path=Path(__file__).resolve().parents[1]/"results"/"author_reference"/"exp_flex_lev16.mat"
        metric, rows = audit_stage7(read_original(path))
        self.assertAlmostEqual(metric["source_max_absolute_budget_error"],
                               .006961346273794,places=10)
        self.assertLess(metric["corrected_max_absolute_budget_error"],1e-10)
        self.assertAlmostEqual(metric["corrected_minus_original_consumption"],
                               .00009567180543346,places=11)
        self.assertAlmostEqual(metric["decomposition_portfolio_revaluation"],
                               .002320487927977,places=10)
        self.assertAlmostEqual(metric["decomposition_individual_budget_error"],
                               -.00009567180543341,places=11)
        self.assertAlmostEqual(metric["binding_max_abs_original_minus_PCHIP_floored"],
                               0.,places=12)
        self.assertEqual(metric["worst_asset_index_1based"],116)
        self.assertEqual(metric["worst_productivity_index_1based"],1)
        self.assertGreater(len(rows),0)

    def test_nonpositive_binding_resources_are_rejected(self):
        with self.assertRaisesRegex(ValueError,"infeasible"):
            budget_consistent_binding(np.array([[.1]]),np.array([[-.1]]),
                                      np.array([-3.]),np.array([1.]),
                                      1.01,-.1,1.5,np.array([[1.]]))


if __name__=="__main__":
    unittest.main()
