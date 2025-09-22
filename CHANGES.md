# Changes for swiftbat_python

## v 0.1.5 2025-05-16

- Moved all conditions for the calculation of the exposure to be inside batExposure
- Made all returned values of batExposure be float values (both the area exposure and the cos(theta))
- Modified simbadlocation to be able to access table RA/Dec values based on lowercase or uppercase headers, which changes with different versions of astroquery
- Use `https` to access heasarc.
- Update TLEs less frequently (10 days is adequate for source visibility).

## v 0.1.4 2023-08-05

- Started CHANGES.md file
- Replace `pyephem` with `skyfield`
  - pyephem is no longer supported
- Pointing history replaced screen-scraping with `swifttool` database access
  - Appeared as DOS on the scraped website when I ran too often
- General decrufting
  - reformatted with black
  - Code still suffers from being the first Python I ever wrote, on python2.6
  - Some (but not much) type hinting)
- Requirements updated:
  - Python >= 3.9
  - astropy >= 5
  - skyfield >= 1.4
