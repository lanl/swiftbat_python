import unittest

import swiftbat
import numpy as np
from pathlib import Path


class TestXY(unittest.TestCase):
    def getsample(self):
        """
        Test data from BAT data

        from astropy.io import fits
        d = fits.getdata('/tmp/sw01090661000bevshto_uf.evt.gz')
        sample = np.array(d[10:1000])
        np.save(open(file, "wb"), sample)
        """
        file = Path(__file__).parent.joinpath("sampleevents.np")
        sample = np.load(open(file, "rb"))
        return sample

    def test_detid2xy_scalar(self):
        sample = self.getsample()
        for samplerow in sample[0:100]:
            x, y = swiftbat.detid2xy(samplerow["DET_ID"])
            assert (x, y) == (samplerow["DETX"], samplerow["DETY"])

    def test_detid2xy_array(self):
        sample = self.getsample()
        x, y = swiftbat.detid2xy(sample["DET_ID"])
        assert np.allclose(x, sample["DETX"])
        assert np.allclose(y, sample["DETY"])

    def test_xy2detid_scalar(self):
        sample = self.getsample()
        for samplerow in sample[0:100]:
            detid = swiftbat.xy2detid(samplerow["DETX"], samplerow["DETY"])
            assert detid == samplerow["DET_ID"]

    def test_xy2detid_array(self):
        sample = self.getsample()
        detid = swiftbat.xy2detid(sample["DETX"], sample["DETY"])
        assert np.allclose(detid, sample["DET_ID"])


if __name__ == "__main__":
    unittest.main()
