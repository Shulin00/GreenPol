#find coordinates of celestial bodies

import sys
##sys.path.append('C:/users/labuser/anaconda/lib/site-packages')
from datetime import datetime
from astropy.coordinates import AltAz, Angle, EarthLocation, ICRS, SkyCoord
from astropy import units as u
import ephem
from numpy import sin,cos,arcsin,arccos,pi

def location_config(LOCATION):
    locations = dict(
          Barcroft  = EarthLocation( lat=Angle(37.5838176, 'deg'),
                                    lon=Angle(-118.2373297, 'deg'), 
                                    height=3800 * u.m),
          Greenland = EarthLocation( lat=Angle(72.5796, 'deg'), 
                                    lon=Angle(-38.4592, 'deg'),
                                    height=3200 * u.m),
          UCSB      = EarthLocation( lat=Angle(34.414, 'deg'),
                                    lon=Angle(-119.843, 'deg'),
                                    height=14 * u.m),
  )
    return locations[LOCATION]

def azalt_to_radec(LOCATION,AZ,ALT):
    Location=location_config(LOCATION)
    time=datetime.utcnow()
    altaz=AltAz(alt=ALT*u.deg,az=AZ*u.deg,location=location,obstime=time)
    frame='icrs'
    frame=astropy.coordinates.frame_transform_graph.lookup_name(frame)()
    radec=altaz.transform_to(frame)
    ra=radec.ra.deg
    dec=radec.dec.deg
    
    return ra,dec
    

def radec_to_azalt(LOCATION,RA,DEC):
    Location=location_config(LOCATION)
    time=datetime.utcnow() 
    radec=SkyCoord(ra=RA*u.deg,dec=DEC*u.deg,frame='icrs',unit='deg')
    altaz=radec.transform_to(AltAz(obstime=time,location=location))
    az=altaz.az.deg
    el=altaz.alt.deg

    return az, el
  
  

def getlocation(LOCATION, CBODY):

  #observation locations
  locations = dict(
          Barcroft  = EarthLocation( lat=Angle(37.5838176, 'deg'),
                                    lon=Angle(-118.2373297, 'deg'), 
                                    height=3800 * u.m),
          Greenland = EarthLocation( lat=Angle(72.5796, 'deg'), 
                                    lon=Angle(-38.4592, 'deg'),
                                    height=3200 * u.m),
          UCSB      = EarthLocation( lat=Angle(34.414, 'deg'),
                                    lon=Angle(-119.843, 'deg'),
                                    height=14 * u.m),
  )

  #celestial bodies
  cbodies = dict(
  	  Sun     = ephem.Sun(),
  	  Moon    = ephem.Moon(),
  	  Mercury = ephem.Mercury(),
      Venus   = ephem.Venus(),
      Mars    = ephem.Mars(),
      Jupiter = ephem.Jupiter(),
      Saturn  = ephem.Saturn(),
      Uranus  = ephem.Uranus(),
      Neptune = ephem.Neptune()
  )

  #observer location
  location = locations[LOCATION]

  #current utc time
  time = str(datetime.utcnow())

  #celestial body of interest
  if type(CBODY)==list:
      az=CBODY[1]
      alt=CBODY[2]
      return az,alt
  if type(CBODY)==str:
      cbody = cbodies[CBODY]
      #this method take 2x as long and produces a slightly different azel, I am not sure which is more accurate
      '''##########################
        #compute radec of body at current time
        cbody.compute(time)

        #create icrs object to convert to az el
        icrs = ICRS(ra = Angle(str(cbody.ra) + 'hours'), dec = Angle(str(cbody.dec) + 'degrees'))
        #print icrs.ra.deg, icrs.dec.deg
        #convert to altaz coordinates
        altaz = icrs.transform_to(AltAz(obstime=t, location=location))
        print altaz.az.deg, altaz.alt.deg
        #############################
      '''
      #set observer location and epoch
      obs = ephem.Observer()
      obs.lon, obs.lat = str(location.longitude.deg), str(location.latitude.deg)
      obs.elevation = float(str(location.height).split()[0])
      obs.date = time
      #compute celestial body az/el coordinates given observer spacetime coordiantes
      cbody.compute(obs)
      print cbody.az
      #convert azimuth to degrees
      az = str(cbody.az).split(':')
      az = [float(i) for i in az]
      az = (az[0] + az[1]/60. + az[2]/60./60.)
      #convert altitude to degrees
      alt = str(cbody.alt).split(':')
      alt = [float(i) for i in alt]
      alt = (alt[0] + alt[1]/60. + alt[2]/60./60.)
      #return azimuth and altitude of celestial body
      return az, alt
if __name__=='__main__':
  getlocation(LOCATION='UCSB',CBODY='Sun',Az=None,El=None)
