"""Stage 7 conservative distribution step unit tests."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from stage7_full_transition_budget_audit import transport_one


class TestConservativeTransport(unittest.TestCase):
    def test_barycenter_and_markov_masses(self):
        fine=np.array([0.,1.,2.,4.])
        mu=np.array([[.25,.0],[.0,.5],[.25,.0],[.0,.0]])
        savings=np.array([[.6,1.],[1.,1.5],[3.,2.],[4.,4.]])
        P=np.array([[.8,.2],[.3,.7]])
        result=transport_one(mu,savings,fine,P)
        self.assertAlmostEqual(result.sum(),1.)
        target=np.sum(mu*savings)
        self.assertAlmostEqual(np.sum(result*fine[:,None]),target)
        self.assertTrue(np.all(result>=0))
        self.assertAlmostEqual(result[:,0].sum(),.5*.8+.5*.3)
        self.assertAlmostEqual(result[:,1].sum(),.5*.2+.5*.7)

    def test_positive_mass_outside_grid_is_error(self):
        fine=np.array([0.,1.,2.])
        mu=np.array([[.4],[.6],[.0]])
        bad=np.array([[.5],[3.],[100.]])
        with self.assertRaisesRegex(ValueError,"outside"):
            transport_one(mu,bad,fine,np.array([[1.]]))

    def test_zero_mass_offgrid_is_ignored(self):
        fine=np.array([0.,1.,2.])
        mu=np.array([[1.],[0.],[0.]])
        sav=np.array([[.5],[3.],[np.nan]])
        out=transport_one(mu,sav,fine,np.array([[1.]]))
        np.testing.assert_allclose(out[:,0],np.array([.5,.5,0.]),atol=1e-14)


if __name__=="__main__":
    unittest.main()
