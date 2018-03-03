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

def uniProp(t,               #Time with unit
            cosmo  = cosmoP, #Cosmology, astropy-style
            Runit  = u.m,    #Unit to display size of Universe in
            output = False   #Output Evernote string
            ):
    """
    Purpose:
        Calculate properties of the Universe at time t.

    History:
        The code was written in conjunction with the (Danish) popular science
        article "Big Bang --- en oejenvidneberetning" ("Big Bang --- an
        eyewitness account"), in order to calculate various properties of the
        Universe. It is meant for being run from the command line, since this
        is more tractable to non-experts, but can also be run from a Python
        environment is a slightly different way

    Usage:
        From the command line:
            > python timeline time unit [-Runit distance_unit]
        From Python environment:
            >>> from astropy import units as u
            >>> import timeline
            >>> timeline.uniProp(t=time*unit [,Runit=distance_unit])

    Arguments (for command line):
        time    Time quantity, i.e. a number
        unit    Unit of time. Allowed values are
                    s:      Seconds
                    min:    Minutes
                    h:      Hours
                    day:    Days
                    yr:     Years
                    kyr:    kilo-years
                    Myr:    mega-years
                    Gyr:    giga-years
    Optional arguments
        -Runit distance_unit  Sensible units for distances are estimated. If you
                              want other units, give the "-Runit" keyword,
                              followed by your desired units, which can be any
                              astropy unit, e.g. angstrom, km,

    Examples:
        > python timeline 1e-32 s            # Properties just after inflation
        > python timeline 13.79 Gyr          # Properties today
        > python timeline 500 Myr -Runit Gpc # Properties 500 million years after
                                             # Big Bang, but use Gpc (giga-parsec)
                                             # for distances

    Same examples from Python environment:
        >>> from astropy import units as u
        >>> import timeline
        >>> timeline.uniProp(1e-32*u.s)
        >>> timeline.uniProp(13.79*u.Gyr)
        >>> timeline.uniProp(500*u.Myr, Runit=u.Gpc)

    Notes:
        z_eq is calculated from the cosmology, and doesn't agree exactly with
        the quoted value in Planck Collaboration et al. (2016).

        For the article, I assumed a Planck(ian) cosmology, but with massless
        neutrinos. The astropy.cosmology.Planck13/15 cosmologies assume massive
        neutrinos, with alters the ratio between photons at late epochs.
        I think this means that matter-radiation equality is calculated
        somewhat differently that "a_eq = (Onu0+Ogamm0)/Om0", but I haven't got
        the time to go into this, and anyway the results are correct to within
        orders of unity, and decreasing with time.
    """


    def photonPressure(T):
        """Photon pressure"""
        return (pi**2 * cc.k_B**4 / (45 * cc.c**3 * cc.h**3) * T**4)


    def nph(T):
        """Photon number density"""
        return 16 * pi * (cc.k_B*T / (cc.h*cc.c))**3 * zeta(3)


    def t_radmat(a,a_eq):
        """Exact t-a relation in a radiation+matter universe.
        Ryden (2003): Eq. 6.37 (p. 94)"""
        import decimal
        from decimal import Decimal as D
        decimal.getcontext().prec = 50
        f1    = 4*a_eq**2 / (3*sqrt(cosmo.Onu0 + cosmo.Ogamma0))
        f2    = float(D(1) - (D(1) - D(a)/(2*D(a_eq))) * (D(1) + D(a)/D(a_eq)).sqrt())
        return (f1*f2 / cosmo.H0).decompose()

    # Calculate raiation-matter equality
    assert (isinstance(t,u.Quantity)) and (t.unit.is_equivalent(u.s)), 'Keyword `t` must have units of time.'
    assert 1e-32*u.s <= t, 't must be > t_endOfInflation'
    a_eq = (cosmo.Onu0 + cosmo.Ogamma0) / cosmo.Om0
    z_eq = 1/a_eq - 1
    t_eq = (4/3. * (1 - 1/sqrt(2)) * (cosmo.Onu0 + cosmo.Ogamma0)**1.5/cosmo.Om0**2 / cosmo.H0).to(u.yr)
    rho_eq = cosmo.critical_density(z_eq)

    # Calculate redshift, density, and temperature
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
        tt  = t_radmat(a,a_eq)
        print(' - Precision:', (tt/t).decompose())
    else:
        print('Matter epoch')
        from astropy.cosmology import z_at_value
        z   = z_at_value(cosmo.age, t, zmax=3500)
        a   = cosmo.scale_factor(z)
        rho = cosmo.critical_density(z)
        T   = cosmo.Tcmb(z)

    # More properties
    H    = cosmo.H(z)                       #Hubble parameter
    R0   = cosmo.comoving_distance(2.7e7)   #Distance in cm to particle horizon today; the 2.7e7 is roughly the highest redshift it can take, but using 1e7, or even 1e6 or 1e4 gives almost the same result
    d_hor = dP(z,cosmo)                     # Horizon distance at t
    d_H   = cc.c / H                        # Hubble radius at t
    Pph  = photonPressure(T)                #Photon pressure
    R    = a * R0                           #Radius of Universe at end of inflation
    rhob = cosmo.critical_density0 * cosmo.Ob0 / a**3   #Baryon density at t
    n_ph = nph(T)                           #Photon number desity 
    X    = .75                              #Primordial hydrogen mass fraction
    Y    = 1 - X                            #Primordial helium
    mu   = 1 / (X + Y/4)                    #Mean molecular mass per electron
    nbar = rhob / (mu*cc.m_p)               #Number density of baryons
    Pbar = nbar * cc.k_B * T                #Baryon pressure
    E    = cc.k_B * T                       #Energy per particle

    # Calculate ionization fraction and resulting mean free path of photons
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

    # Print results
    u.c = 2.99792458e10 * u.cm / u.s
  # c = u.def_unit('c', 2.99792458e10 * u.cm / u.s)
    print('Age at rad-mat eq.:       {:.0f}'.format(cosmo.age(z_eq).to(u.yr)))
    print('Redshift at rad-mat eq.:  {:.0f}'.format(z_eq))
    print('Radius of Universe today: {:2.2f}'.format(R0.to(u.Glyr)))
    print()
    print('Physical properties at t = {:}:'.format(t))
    print(' * Expansion:')
    print('   - Scale factor:                         {:.3g}'.format(a))
    print('   - Redshift:                             {:.4g}'.format(z))
    print('   - Hubble parameter:                     {:.3g}'.format(H))
 #  print('     - H(t) in c / m:                      {:.3g}'.format(H.to(u.c/u.m).value))
    print('   - e-foldings since t:                   {:2.1f}'.format(np.log(1/a)))
    print(' * Size:')
    print('   - Radius of observable Universe at t:   {:.3g}'.format(d_hor.to(Runit)))
    print("   - Radius of today's obs. Universe at t: {:.3g}".format(R.to(Runit)))
    print('   - Hubble distance:                      {:.3g}'.format(d_H.to(Runit)))
    print(' * Gas and radiation:')
    print('   - Temperature:                          {:.1e}'.format(T))
    print('   - Energy:                               {:.3g}'.format(E.to(u.MeV)))
    print('   - Energy density:                       {:.3g}'.format(rho))
    print('   - Ionized fraction:                     {:.3g}'.format(xe))
    print('   - Photon mean free path:                {:.3g}'.format(mfp.to(Runit)))
    print('     - mfp / d_H                           {:.3g}'.format((mfp / d_H).decompose()))
    print('   - Photon no. density:                   {:.3g}'.format(n_ph.to(u.cm**(-3))))
    print('   - Baryon no. density:                   {:.3g}'.format(nbar.to(u.cm**(-3))))
    print('   - Photon pressure:                      {:.3g}'.format(Pph.to(u.cds.atm)))
    print('   - Baryon pressure:                      {:.3g}'.format(Pbar.to(u.cds.atm)))

    # Print one-lines (for timeline table)
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


def mue(X):
    """Mean molecular mass per electron"""
    return 2. / (1+X)


def ne(rhob,xe,X):
    """Number density of electrons"""
    return xe * rhob / (mue(X)*cc.m_p)


def meanFreePath(rhob,xe,X):
    """Mean free path to Thomson scattering"""
    return 1 / (ne(rhob,xe,X) * cc.sigma_T)


def Saha(xe,T,cosmo):
    """Saha equation solver"""
    return xe**2/(1-xe) - 5.8e15/(cosmo.Ob0*cosmo.h**2*(T/1e4)**1.5) * exp(-157000/T)


def dP(z,cosmo):
    """
    Particle horizon
    """
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
    parser.add_argument('time',               type=float, help='Time quantity')
    parser.add_argument('unit',               type=str,   help='Time unit (s,min,day,yr,kyr,Myr,Gyr)')
    parser.add_argument('-Runit', default='', type=str,   help='Unit for distances (not exactly sure why this works)')
    args = parser.parse_args()

    if args.unit == 's':
        time  = args.time * u.s
        Runit = u.cm
    elif args.unit == 'min':
        time  = args.time * u.min
        Runit = u.AU
    elif args.unit == 'day':
        time  = args.time * u.day
        Runit = u.lyr
    elif args.unit == 'yr':
        time  = args.time * u.yr
        Runit = u.lyr
    elif args.unit == 'kyr':
        time  = args.time * u.kyr
        Runit = u.Mlyr
    elif args.unit == 'Myr':
        time  = args.time * u.Myr
        Runit = u.Mlyr
    elif args.unit == 'Gyr':
        time  = args.time * u.Gyr
        Runit = u.Glyr
    else:
        print("Sorry, unit `"+args.unit+"` is not implemented.\nYou're welcome to go ahead and do it yourself,\nand then send me a pull request.")
        exit()

    if args.Runit != '':
        Runit = args.Runit

    uniProp(t=time, Runit=Runit)

if __name__ == '__main__':
    main()
