"""
Functions for calculating 1) initial values, 
                          2) sun angle, 
                          3) vertical mixing coefficients, 
                          4) surface deposition, 
                          5) vertical mixing
"""

import numpy as np
import math

__all__ = ("init", "zenith", "atk", "sinks", "newc")

def atk(H2O, nlev, theta, t_prof, p_prof, vcp, wind, z):

  """Calculate vertical mixing coefficients                 
  
  Parameters:
  -----------
  H2O: int
    Index for H2O concentration in vcp
  
  nlev: int
    Number of height levels
    
  theta: list
    
  t_prof: list
    Temperature profile (K)
  
  p_prof: list
    Pressure profile (hPa)
 
  vcp: 2-D list
    Concentration of all species
    
  wind: list
    Wind profile
 
  z: list
    Heights (meter)
  	
  Return:
  -----------
  akh: list
  
  theta: list
  
  """

  ak = 0.4
  grav = 9.81
  airmol = 0.1529e-04
  thelam = 30
  brcr = 0.25
  xx0 = 0
  xxh = 0
  kpbl = 0
  
  stable = False
  knov = False
  
  akh = np.zeros(nlev)
  ajnn = np.zeros(nlev)
  akhst = np.zeros(nlev)
  akm = np.zeros(nlev)
  akmst = np.zeros(nlev)
  akm_ysu = np.zeros(nlev)
  thvx = np.zeros(nlev)
  ri_pro = np.zeros(nlev)
  
  for i in range(nlev):
    theta[i] = t_prof[i] * ((100000 / p_prof[i])**0.286)
      
    im = max(i - 1, 0)
    ip = min(i + 1, nlev - 1)
      
    akhst[i] = 0.05 * akh[i] + 0.9 * akh[im] + 0.05 * akh[ip] + airmol
    akmst[i] = 0.05 * akm[i] + 0.9 * akm[im] + 0.05 * akm[ip] + airmol
      
  for i in range(nlev - 1):
    im = max(i - 1, 0)
    ip = min(i + 1, nlev - 1)
      
    dz_xhu = z[i + 1] - z[i]
    strat = (theta[i + 1] - theta[i]) / dz_xhu
      
    if(strat <= 0):
      xx0 = xx0 + (z[ip] - z[i])
      knov = (theta[i] + 0.1) < theta[0]
        
      if(knov):
        if(xxh <= 0):
          xxh = z[i]
          
    else:
      knov = False
    
  stable_xhu = False
    
  for i in range(nlev):
    thvx[i] = theta[i] * (1 + 0.61 * vcp[i, H2O] * 0.622)
      
  for i in range(15, nlev - 1):
    if(not stable_xhu):
      spdk2 = max(wind[i]**2, 1)
      brup = (thvx[i] - thvx[14]) * (9.8 * (z[i] - z[14]) / thvx[14]) / spdk2
      kpbl = i
      stable_xhu = brup > brcr
        
  for i in range(nlev - 1):
    dz_xhu = z[i + 1] - z[i]
    strat = (theta[i + 1] - theta[i]) / dz_xhu
    vzbet = (wind[i + 1] - wind[i]) / dz_xhu
      
    zm = 0.5 * (z[i] + z[i + 1])
    tm = 0.5 * (t_prof[i] + t_prof[i + 1])
      
    vzbet2 = vzbet * vzbet
      
    ri = grav * strat / (tm * vzbet2)
    ri_pro[i] = ri
      
    xlam = max(thelam, 0.1 * xx0, 0.1 * xxh)
    
    if(zm > max(xx0, xxh)):
      xlam = thelam
        
    All = ak * zm / (1 + ak * zm / xlam)
      
    ajnn[i] = (All**2) * math.sqrt(vzbet2)
      
    if(ri < 0):
      akm[i] = ajnn[i] * math.sqrt(1 - 11 * ri) + airmol
      akh[i] = akm[i] * 1.35 * (1 - 5.5 * ri) / (1 - 3 * ri) + airmol
      
    else:
      ri = min(ri, 1)
      akm[i] = ajnn[i] / math.sqrt(1 + 6 * ri) + airmol
        
      Ustar = 0.4
      Z_over_L = 0.1
        
      akm_ysu[i] = (0.4 * Ustar / (1 + 5 * Z_over_L) * z[i] * ((1 - z[i] / z[kpbl])**3))
                    
      if(z[kpbl] < 400 and z[i] < z[kpbl]):
        akm[i] = max(akm[i], akm_ysu[i])
        
      if(stable):
        akh[i] = akm[i] * 1.35 / (1 + 6 * ri) + airmol
        
      else:
        akh[i] = akm[i]
        
    akm[i] = 0.1 * akmst[i] + 0.9 * akm[i]
    akh[i] = 0.1 * akhst[i] + 0.9 * akh[i]
    
  return akh, theta

def declin(decmax, icumdy, iyear, pid180):
  """Calculate solar daclination
  
  Parameters:
  -----------
  decmax: float
  
  icumdy: int
    Starting day of staring year
  
  iyear: int
    Starting year
    
  pid180: float
    pi divided by 180
  
  Return:
  -----------
  sindec: float
  
  cosdec: float
  
  eqtm: float
  
  """

  kday = (iyear - 1977) * 365 + icumdy + 28123
  xm = (-1 + 0.9856 * kday) * pid180
  delnu = 2 * 0.01674 * math.sin(xm) + 1.25 * 0.01674 * 0.01674 * math.sin(2*xm)
  slong = (-79.828 + 0.9856479 * kday) * pid180 + delnu
  decl = math.asin(decmax * math.sin(slong))
  sindec = math.sin(decl)
  cosdec = math.cos(decl)
  eqtm = (9.4564 * math.sin(2 * slong) / cosdec - 4 * delnu / pid180) / 60
  
  return sindec, cosdec, eqtm

def init(Temp, iyear, jday, stdlng, xlong, xlat, nlev, init_values, H2O):
  """Generate initial values                 
  
  Parameters:
  -----------
  Temp: boolean
    True: use nlev-based function to calculate temperature
    False: all nlevs have the same temperature
  
  iyear: int
    Starting year
    
  jday: int
    Starting day of starting year
    
  stdlng: float
    Standard longitutde is 75 plus 15 for daylight saving time
  
  xlong: float
    Longitude
    
  xlat: float
    Latitude
  
  nlev: int
    Number of height levels
    
  init_values: list
    List of initial values for all species
  	
  Return:
  -----------
  dlong: float
  
  sinlat: float
  
  coslat: float
  
  sindec: float
  
  cosdec: float
  
  eqtm: float
  
  t_prof: list
    Temperature profile (K)
  
  theta: list
  
  p_prof: list
    Pressure profile (hPa)
 
  z: list
    Heights (meter)
  
  dz: list
    Difference in heights between heights (meter)
  
  dzf: list
  
  wind: list
    Wind profile
  
  relh: list
    Humidity
  
  """

  pres = 1013.25*100 # Surface presssure, pascal. TODO: Get surface pressure from user input
  pid180 = 3.1415926535/180
  decmax = math.sin(23.44 * pid180)

  z = np.zeros(nlev)
  zf = np.zeros(nlev)
  dz = np.zeros(nlev)
  dzf = np.zeros(nlev)
  t_prof = np.zeros(nlev)
  theta = np.zeros(nlev)
  p_prof = np.zeros(nlev)
  relh = np.zeros(nlev)
  vcp = np.zeros((nlev, len(init_values)))

  icumdy = jday
  
  dlong = (stdlng - xlong) / 15
  sinlat = math.sin(xlat * pid180)
  coslat = math.cos(xlat * pid180)

  # Solar daclination
  sindec, cosdec, eqtm = declin(decmax, icumdy, iyear, pid180) 

  # Altitude calculation
  for i in range(nlev):
  
    if(i < 9):
      z[i] = (i + 1) * 0.1
    elif(i < 47):
      z[i] = 1 + ((i + 1) - 10) * 0.5
    elif(i < 127):
      z[i] = 20 + ((i + 1) - 48)
    elif(i < 188):
      z[i] = 100 + ((i + 1) - 128) * 5
    elif(i < 218):
      z[i] = 400 + ((i + 1) - 188) * 20
    else:
      z[i] = 1000 + ((i + 1) - 218) * 100
  
  zf[nlev - 1] = (z[nlev - 1] + 0.5 * (z[nlev - 1] - z[nlev - 2]))
                    
  for i in range(nlev - 1):
    zf[i] = 0.5 * (z[i + 1] + z[i])
    
    if(i > 0):
      dzf[i] = zf[i] - zf[i - 1]
    
    dz[i] = z[i + 1] - z[i]
  
  dzf[0] = zf[0]
  dzf[nlev - 1] = dzf[nlev - 2]  
  dz[nlev - 1] = dz[nlev - 2]

  # Atmospheric variables
  for i in range(nlev):
    p_prof[i] = pres * (((1 - (z[i]) / 44308))**(1 / 0.19023)) # Pressure
    
    if Temp: # Calculate temperature using a function
      t_prof[i] = 273 - z[i] * 0.0085 #0.0055
      
      if z[i] > 1500: 
        t_prof[i] = (300 -1500 * 0.005) - (z[i] - 1500) * 0.002
      
      if z[i] < 28:
        t_prof[i] = 290.00 + z[i] * 0.20 
    
    else: # Temperature is a constant
      t_prof[i] = 273
    
    theta[i] = t_prof[i] * ((100000 / p_prof[i])**0.286)
    
    wind = wprofil(nlev, z)
  
    # Humidity
    relh[i] = 0.7
    
    if(z[i] < 500):
      relh[i] = 0.7 - z[i] * 0.2/500
  
    # Assign initial values
    for k in range(len(init_values)):
      vcp[i, k] = init_values[k]  
    
    vcp[i, H2O]   = (relh[i] * 610.7 * math.exp(17.1536 * (t_prof[i] - 273.15) / (t_prof[i] - 38.33))) / p_prof[i]
    
  return (dlong, sinlat, coslat, sindec, cosdec, eqtm, 
      t_prof, theta, p_prof, z, dz, dzf, wind, relh, vcp)

def newc(nlev, akh, dzf, dz, deltim, nspec_host, vcp, source, H2O, all_spc):
  
  """Calculate surface deposition
  
  Parameters:
  -----------
  nlev: int
    Number of height levels
  
  akh: list
  
  dzf: list
  
  dz: list
    Difference in heights between heights (meter)
  
  deltim: int
  
  nspec_host: int
    Number of species
  
  vcp: 2-D list
    Concentration of all species
    
  source: 2-D list
    Constant source/sink terms for species
  
  H2O: int
    Index for H2O concentration in vcp
    
  Return:
  -----------
  vcp: 2-D list
    Concentration of all species
  
  """
  
  nr = 1
    
  a0 = np.zeros(nlev)
  b0 = np.zeros(nlev)
  c0 = np.zeros(nlev)
  cs = np.zeros(nlev)
  x  = np.zeros(nlev)
    
  for i in range(1, nlev - 1):
      
    a0[i] = -akh[i - 1] / dzf[i] / dz[i - 1] * deltim
    c0[i] = -akh[i] / dzf[i] / dz[i] * deltim
    b0[i] = 1 + (-a0[i] - c0[i]) 
      
  for i in range(nspec_host):
    
    if i == 51:
      continue
        
    for j in range(nlev):
      cs[j] = vcp[j, i]
      
    for j in range(1, nlev - 1):
      x[j] = cs[j] + source[j, i] * deltim
          
    at = 0
    bt = 1
    dt = cs[nlev - 1]
    ab = 1
    bb = 0
    db = cs[0] + source[0, i] * deltim
      
    if(i == H2O):
      db = cs[0]
    
    if "O3P" in all_spc: 
      if(i != all_spc.index("O3P")):
        cs = solve(ab, bb, db, at, bt, dt, a0, b0, c0, x, cs, nr, nlev)
      
    for j in range(nlev):
      vcp[j, i] = max(cs[j], 0)
        
    vcp[0, i] = vcp[0, i] * 0.1 + vcp[1, i] * 0.9
    
  return vcp
    

def sinks(akh, f0, henry, O3, nlev, nspec_host, relh, source, vcp, dz, z):
  
  """Calculate surface deposition
  
  Parameters:
  -----------
  akh: list
  
  f0: list
  
  henry: list
    Henry's law constant (to be confirmed)
  
  O3: int
    Index for O3 concentration in vcp
  
  nlev: int
    Number of height levels
  
  nspec_host: int
    Number of species
  
  relh: list
    Humidity
  
  source: 2-D list
    Constant source/sink terms for species
  
  vcp: 2-D list
    Concentration of all species
  
  dz: list
    Difference in heights between heights (meter)
  
  z: list
    List of heights
    
  Return:
  -----------
  Vdep: list
  
  source: 2-D list
    Constant source/sink terms for species
  
  """
    
  rcutref = 3000
    
  rgs = 500
  rgo = 200
    
  f0[O3] = 1
    
  if(relh[0] > 0.8):
    rgo = 200 + 1800 * (relh[0] - 0.8) / (0.999 - 0.8)
    f0[O3] = 1 - 0.9 * (relh[0] - 0.8) / (0.999 - 0.8)
      
  Vdep = np.zeros(nspec_host)
    
  for i in range(nspec_host - 9): # only for gas species
    rgi = 1 / (henry[i] * 1e-5 / rgs + f0[i] / rgo)
    ra0 = z[1] / max(1e-4, akh[0])
      
    sh1 = 1 / (ra0 + rgi)
    Vdep[i] = sh1
    si = -sh1 * vcp[0, i] / dz[0]
    source[0, i] = source[0, i] + si
      
  return Vdep, source


def solve(ab, bb, db, at, bt, dt, a, b, c, d, var, nr, nlev):
  
  alpha = np.zeros(nlev)
  beta = np.zeros(nlev)
  xx = np.zeros(nlev)
  
  alpha[nr - 1] = -bb / ab
  beta[nr - 1] = db / ab
    
  m1 = nlev - 1
    
  for j in range(nr, m1):
    xx[j] = (a[j] * alpha[j-1] + b[j])
    alpha[j] = -c[j] / xx[j]
    beta[j] = -(a[j] * beta[j-1] - d[j]) / xx[j]
      
  var[nlev - 1] = (dt - at * beta[m1 - 1]) / (bt + alpha[m1 - 1] * at)
    
  for jj in range(nr - 1, m1):
    j = nlev - 2 - jj
    var[j] = var[j + 1] * alpha[j] + beta[j]
       
    if(var[j] < 0):
      var[j] = 0.0001 * var[j+1]
        
  return var

def wprofil(nlev, z):
  """Calculate wind profile
  
  Parameters:
  -----------
  nlev: int
    Number of height levels
  
  z: list
    List of heights
    
  Return:
  -----------
  wind: list
    Wind profile
  
  """
  
  wind = np.zeros(nlev)
  
  vg = 7.0
  ustara = vg / 30
  wind[0] = 0
  
  for i in range(1, nlev):
    wind[i] = ustara * math.log((z[i]) / 0.01) / 0.4
    
  return wind

def zenith(timloc, eqtm, dlong, sinlat, sindec, coslat, cosdec):
  """Calculate Sun Angle
  
  Parameters:
  -----------
  timloc: int
  
  eqtm: float
  
  dlong: float
  
  sinlat: float
  
  sindec: float
  
  coslat: float
  
  cosdec: float
  
  Return:
  -----------
  coszen: float
  
  zenang: float
  """
  
  crtzen = 0
  pid = 3.1415926535
  pid2 = 3.1415926535/2
  
  timsun = timloc + eqtm + dlong
  hrang = (timsun - 12) * pid2 / 6
  zenang = math.acos(sinlat * sindec + coslat * cosdec * math.cos(hrang))
  sunazm = math.asin(cosdec * math.sin(hrang) / math.sin(zenang))
  
  if(sindec > sinlat):
    crtzen = math.acos(1)
  else:
    crtzen = math.acos(sindec / sinlat)
    
  if(zenang > crtzen):
    sunazm = (pid - abs(sunazm)) * sunazm / abs(sunazm)
  
  sunazm = sunazm + pid
  coszen = math.cos(zenang)
  
  if(timloc == 4.5):
    print("Check coszen: ", coszen)
    
  return coszen
