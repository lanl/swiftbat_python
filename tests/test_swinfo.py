import datetime
import unittest

import swiftbat


class SwInfoTestCase(unittest.TestCase):
    def test_source(self):
        s1020 = swiftbat.source(ra=10, dec=20)
        sunnow = swiftbat.sunsource()
        sun = swiftbat.sunsource("2023-03-21 06:00")  # Near the Vernal equinox
        crab = swiftbat.source("Crab")
        cygx1 = swiftbat.source("Cyg_x-1")
        for source in [s1020, sun, crab, cygx1]:
            print(f"{source.name:10s} {source.ra_deg:8.3f} {source.dec_deg:8.3f}")
        pass

    def test_swinfo_cli(self):
        swiftbat.swinfo_main(
            "swinfo_test  2020-05-05T12:34:56 -o -s cyg_X-1 -S".split()
        )
        """
        swinfo  2020-05-05T12:34:56 -o -s cyg_X-1 -S
        Swift Zenith(RA,dec):         232.97, -20.46
        Swift Location(lon,lat,alt):  -179.31 E, -20.53 N, 549 km
        Time(MET + UTCF -> UT):   610374920.889 + -24.888 -> 2020-05-05T12:34:55.999
        YYMMDD_SOD:                   200505_45295.999       DOY=126
        Obs Sequence Number:          00033349058
        Obs Target Name:              SGR 1935+2154
        Obs Date and Times:           2020-05-05   12:32:37 - 13:01:34
        Obs Pointing(ra,dec,roll):    293.724, 21.897, 62.77
        Cyg_X-1_imageloc (boresight_dist, angle): (14, -133)
        Cyg_X-1_exposure (cm^2 cos adjusted): 4230
        Cyg_X-1_altitude: 29 (up)
        Sun_imageloc (boresight_dist, angle): (101, -85)
        Sun_exposure (cm^2 cos adjusted): 0
        Sun_altitude: -57 (down)
        """


if __name__ == "__main__":
    unittest.main()
