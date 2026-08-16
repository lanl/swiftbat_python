import datetime
import unittest

import swiftbat
from swiftbat.clockinfo import clockErrData
from pathlib import Path
import tempfile
import datetime


class ClockTestCase(unittest.TestCase):
    def test_utcf(self):
        testvector = [
            (swiftbat.string2met("2018-08-29 23:31:01", correct=False), -21.50867)
        ]
        clockerrdata = clockErrData()
        clockerrfile = Path(clockerrdata._clockfile)
        self.assertTrue(clockerrfile.exists())
        for t, utcf_expected in testvector:
            utcf_returned = swiftbat.utcf(t)
            self.assertAlmostEqual(utcf_returned, utcf_expected, places=3)

    def test_clock_download(self):
        clockerrdata = clockErrData()
        with tempfile.TemporaryDirectory() as tempclockdir:
            clockerrdata.updateclockfiles(tempclockdir, ifolderthan_days=0)
            file = next(Path(tempclockdir).glob("swclock*.fits"))
            self.assertTrue(file.exists())


if __name__ == "__main__":
    unittest.main()
