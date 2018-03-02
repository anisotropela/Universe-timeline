import argparse
import numpy as np
from numpy import sqrt,pi,exp
import scipy as sp
from astropy.cosmology import FlatLambdaCDM, Planck15, WMAP9
from astropy import units as u
from astropy.units import cds
import matplotlib.pyplot as plt
import astropy.constants as cc
from scipy.special import zeta
from scipy.optimize import newton
cosmoP = FlatLambdaCDM(H0       = 67.81,
                       Om0      = 0.308,
                       Ob0      = .0484,
                     # Onu0     = 3.710122469245305e-05,
                     # Ogamma0  = 5.373825182529615e-05,
                       name     ='Planck')
"""
Purpose:
    Calculate properties of the Universe at time t

Usage:
    From the command line:
    > python timeline time unit

Arguments:
    time    Time quantity
    unit    Unit of time. Allowed values are
                s:      Seconds
                min:    Minutes
                h:      Hours
                day:    Days
                yr:     Years
                kyr:    kilo-years
                Myr:    mega-years
                Gyr:    giga-years

Examples:
    > python timeline 1e-32 s    # Properties just after inflation
    > python timeline 13.79 Gyr  # Properties today

Notes:
    To get sensible units, you may need to change manually in the print statements.
"""

def uniProp(t,               #Time with unit
            cosmo  = cosmoP, #Cosmology, astropy-style
            Runit  = u.m,    #Unit to display size of Universe in
            output = False   #Output Evernote string
            ):
    """
    Properties of the Universe at time t.
    Note that z_eq is calculated from the cosmology, and doesn't agree exactly
    with the quoted value in Planck (2016).

    Usage:
    >>> from astropy import units as u
    >>> uniProp(t=1e-32*u.s) # Properties just after inflation
    """


    def photonPressure(T):  #Photon pressure
        return (pi**2 * cc.k_B**4 / (45 * cc.c**3 * cc.h**3) * T**4)


    def nph(T):             #Photon number density
        return 16 * pi * (cc.k_B*T / (cc.h*cc.c))**3 * zeta(3)


    def t_radmat(a,a_eq):   #Exact t-a relation in a radiation+matter universe
        """Ryden (2003): Eq. 6.37 (p. 94)"""
        import decimal
        from decimal import Decimal as D
        decimal.getcontext().prec = 50
        f1    = 4*a_eq**2 / (3*sqrt(cosmo.Onu0 + cosmo.Ogamma0))
        f2    = float(D(1) - (D(1) - D(a)/(2*D(a_eq))) * (D(1) + D(a)/D(a_eq)).sqrt())
        return (f1*f2 / cosmo.H0).decompose()

    assert (isinstance(t,u.Quantity)) and (t.unit.is_equivalent(u.s)), 'Keyword `t` must have units of time.'
    a_eq = (cosmo.Onu0 + cosmo.Ogamma0) / cosmo.Om0
    z_eq = 1/a_eq - 1
    t_eq = (4/3. * (1 - 1/sqrt(2)) * (cosmo.Onu0 + cosmo.Ogamma0)**1.5/cosmo.Om0**2 / cosmo.H0).to(u.yr)
    assert 1e-32*u.s <= t, 't must be > t_endOfInflation'

    R0     = cosmo.comoving_distance(2.7e7)  #Distance in cm to particle horizon today; the 2.7e7 is roughly the highest redshift it can take, but using 1e7, or even 1e6 or 1e4 gives almost the same result
    rho_eq = cosmo.critical_density(z_eq)
    rho0   = cosmo.critical_density0

    t_thres = 20 * u.kyr # Threshold between a(t) schemes
    if t <= t_eq:
        print('Photon epoch')
        if t <= t_thres:
            a = sqrt(2 * sqrt(cosmo.Onu0+cosmo.Ogamma0) * (cosmo.H0 * t).decompose()).value
        else:
            a = (a_eq * sqrt(t/t_eq).decompose()).value
        rho = rho_eq * (a_eq/a)**4
        z   = 1./a - 1
        T   = cosmo.Tcmb0 / a
    else:
        print('Matter epoch')
        from astropy.cosmology import z_at_value
        z   = z_at_value(cosmo.age, t, zmax=3500)
        a   = cosmo.scale_factor(z)
        rho = cosmo.critical_density(z)
        T   = cosmo.Tcmb(z)

    Pph  = photonPressure(T)
    R    = a * R0                          #Radius of Universe at end of inflation
    rhob = cosmo.critical_density0 * cosmo.Ob0 / a**3
    n_ph = nph(T)
    H    = cosmo.H(z)
    X    = .75
    Y    = 1 - X
    mu   = 1 / (X + Y/4)
    nbar = rhob / (mu*cc.m_p)
    Pbar = nbar * cc.k_B * T
    E    = cc.k_B * T

    if T > 4500*u.K:
        xe = 1.
    elif 4500*u.K >= T > 500*u.K:
        if T > 2500*u.K:
            x0 = .999
        elif 2500*u.K >= T > 100*u.K:
            x0 = .5
        elif 100*u.K >= T > 50*u.K:
            x0 = 1e-5
        xe  = newton(Saha,x0,args=(T.value,cosmo))
    elif T <= 500*u.K:
        xe = 1e-10
    mfp   = meanFreePath(rhob,xe,X)
    d_hor = dP(z,cosmo)

    u.c = 2.99792458e10 * u.cm / u.s
  # c = u.def_unit('c', 2.99792458e10 * u.cm / u.s)

    print('Age at rad-mat eq.:       {:.0f}'.format(cosmo.age(z_eq).to(u.yr)))
    print('Redshift at rad-mat eq.:  {:.0f}'.format(z_eq))
    print('Radius of Universe today: {:2.2f}'.format(R0.to(u.Glyr)))
    print('Physical properties at t:')
    print(' - Scale factor:       {:.3g}'.format(a))
    print(' - Redshift:           {:.4g}'.format(z))
    print(' - Hubble parameter:   {:.3g}'.format(H))
    print('   - H(t) in c / m:    {:.3g}'.format(H.to(u.c/u.m).value))
    print(' - e-foldings since t: {:2.1f}'.format(np.log(1/a)))
    print(' - Dist. to par. hor.  {:.3g}'.format(d_hor.to(Runit)))
    print(' - Radius:             {:.3g}'.format(R.to(Runit)))
    print(' - Diameter:           {:.3g}'.format(2*R.to(Runit)))
    print(' - Density:            {:.3g}'.format(rho))
    print(' - Photon no. density: {:.3g}'.format(n_ph.to(u.cm**(-3))))
    print(' - Naryon no. density: {:.3g}'.format(nbar.to(u.cm**(-3))))
    print(' - Ionized fraction    {:.3g}'.format(xe))
    print(' - Photon mfp:         {:.3g}'.format(mfp.to(u.pc)))
    print('   - mfp / d_hor       {:.5g}'.format((mfp / d_hor).decompose()))
    print(' - Photon pressure:    {:.3g}'.format(Pph.to(u.cds.atm)))
    print(' - Baryon pressure:    {:.3g}'.format(Pbar.to(u.cds.atm)))
    print(' - Temperature:        {:.1e}'.format(T))
    print(' - Energy:             {:.3f}'.format(E.to(u.MeV)))

    tt = t_radmat(a,a_eq)
    print('-----------')
    print('Precision:', (tt/t).decompose())

    if output:
        numbers = [t,
                   T,
                   E.to(u.eV),
                   d_hor.to(u.Mpc),
                   R.to(u.Mpc),
                   H.to(u.km/u.s/u.Mpc),
                   rho,
                   a,
                   z,
                   n_ph.to(u.cm**(-3)),
                   mfp.to(u.Gpc),
                   Pph.to(u.cds.atm)]
        print('              {:.3g}     {:.3g}                {:.3g}  {:.3g} {:.3g}  {:.3g}  {:.3g}                                     {:.3g}   {:.3g}   {:.3g}      {:.3g}    {:.3g}'.format(*numbers))


def mue(X):             #Mean molecular mass per electron
    return 2. / (1+X)


def ne(rhob,xe,X):      #Number density of electrons
    return xe * rhob / (mue(X)*cc.m_p)


def meanFreePath(rhob,xe,X):     #Mean free path to Thomson scattering
    return 1 / (ne(rhob,xe,X) * cc.sigma_T)


def Saha(xe,T,cosmo):
    return xe**2/(1-xe) - 5.8e15/(cosmo.Ob0*cosmo.h**2*(T/1e4)**1.5) * exp(-157000/T)


def dP(z,cosmo):
    from astropy.cosmology import z_at_value

    a_eq = (cosmo.Onu0 + cosmo.Ogamma0) / cosmo.Om0
    z_eq = 1/a_eq - 1
    t_eq = (4/3. * (1 - 1/sqrt(2)) * (cosmo.Onu0 + cosmo.Ogamma0)**1.5/cosmo.Om0**2 / cosmo.H0).to(u.s).value

    def inva(t):
        if t <= t_eq:
            a = sqrt(2 * sqrt(cosmo.Onu0+cosmo.Ogamma0) * (cosmo.H0 * t).decompose()).value
        else:
            z = z_at_value(cosmo.age, t*u.s, zmax=1e4)
            a = cosmo.scale_factor(z)
        return 1 / a

    eta,err = sp.integrate.quad(inva,0,cosmo.age(z).to(u.s).value) * u.s

    return cc.c * eta / (1+z)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('time', type=float, help='Time quantity')
    parser.add_argument('unit', type=str,   help='Time unit (s,min,day,yr,kyr,Myr,Gyr)')
    args = parser.parse_args()

    if args.unit == 's':
        time = args.time * u.s
    elif args.unit == 'min':
        time = args.time * u.min
    elif args.unit == 'day':
        time = args.time * u.day
    elif args.unit == 'yr':
        time = args.time * u.yr
    elif args.unit == 'kyr':
        time = args.time * u.kyr
    elif args.unit == 'Myr':
        time = args.time * u.Myr
    elif args.unit == 'Gyr':
        time = args.time * u.Gyr

    uniProp(t=time, Runit=u.Glyr)

if __name__ == '__main__':
    main()
