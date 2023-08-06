# Changes for swiftbat_python

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
