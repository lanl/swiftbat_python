"""
Find and manage Swift data files from Heasarc
"""

from astroquery.heasarc import Heasarc
import datetime
from swiftbat.swutil import any2datetime, datetime2mjd


def heasarcquery(cache=True, **queryargs):
    heasarc = Heasarc()
    if False:
        queryargs.setdefault('position','0.0d 0.0d')
        queryargs.setdefault('radius','361 deg')
        payload = heasarc.query_region(get_query_payload=True, **queryargs)
        for field in ['Entry', 'Coordinates', 'Radius']:
            payload.pop(field)
        result = heasarc.query_async(payload, cache=cache)
        result = heasarc._parse_result(result)
    else:
        object_name = queryargs.pop('object_name', None)
        result = heasarc.query_object(object_name=object_name, cache=cache, **queryargs)
    return result
    

def triggerrecords(time=None, timewindow_d=1/24, trignum=None, cache=True, **queryargs):
    """Get trigger records from the 'swifttdrss' table of HEASARC

    Table definition:
    https://heasarc.gsfc.nasa.gov/W3Browse/all/swifttdrss.html
    Do your own Quality Control on the results

    time can be a single value (which will be extended by timewindow_d days 
        on each side) or a list, where the extreme values will be used as a range
    trignum can be a single value, or a list where the extreme values 
        will be used as an inclusive range
    
    time gets converted to an MJD to match the 'TIME' column of the result
    trignum goes to the TARGET_ID column.
    
    Any other arguments are passed to the query
    """
    if not queryargs.pop('force', False):            
        if trignum is None and time is None and len(queryargs) == 0:
            raise RuntimeError("If you want the full table, set force=True")
        
    queryargs.setdefault('mission', 'swifttdrss')
    queryargs.setdefault('fields', 'all')
    # queryargs.setdefault('entry', '')
    if time is not None:
        time = any2datetime(time)
        if not isinstance(time, list):
            time = [time]     # Start and end
        mjds = [datetime2mjd(t) for t in time]
        if len(mjds) == 1:
            mets = [mjds[0] - timewindow_d,
                    mjds[0] + timewindow_d]
        queryargs['TIME'] = f"{min(mjds)} .. {max(mjds)}"
    if trignum is not None:
        try:
            trignum = [tr_ for tr_ in trignum]
        except TypeError:
            trignum = trignum
            pass
        queryargs['target_id'] = (f"{trignum[0]}" if len(trignum) == 1 
                                  else f"{min(trignum)} .. {max(trignum)}")
    return heasarcquery(cache=cache, **queryargs)

def main_test():
    result_trigrange = triggerrecords(trignum=[1104735, 1104873])
    result_timerange = triggerrecords(time=["2022-05-01", "2022-05-02"])
    pass
    
if __name__ == '__main__':
    main_test()
    