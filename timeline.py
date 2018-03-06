import argparse
import numpy as np
from numpy import sqrt,pi,exp
import scipy as sp
from astropy.cosmology import FlatLambdaCDM,WMAP5,WMAP7,WMAP9,Planck13,Planck15
from astropy import units as u
from astropy.units import cds
import matplotlib.pyplot as plt
import astropy.constants as cc
from scipy.special import zeta
from scipy.optimize import newton
mycosmo = FlatLambdaCDM(H0       = 67.81,
                        Om0      = 0.308,
                        Ob0      = .0484,
                      # Onu0     = 3.710122469245305e-05,
                      # Ogamma0  = 5.373825182529615e-05,
                        name     ='My cosmology')

def uniProp(t,                  #Time with unit
            cosmo  = Planck15,  #Cosmology, astropy-style
            Runit  = u.m,       #Unit to display size of Universe in
            output = False      #Output Evernote string
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
            > python timeline.py time unit [-Runit my_dist_unit]
        From Python environment:
            >>> from astropy import units as u
            >>> import timeline
            >>> timeline.uniProp(t=time*unit [,Runit=my_dist_unit])

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
    Optional arguments:
        -Runit my_dist_unit  Sensible units for distances are estimated. If you
                             want other units, give the "-Runit" keyword,
                             followed by your desired units, which can be any
                             astropy unit, e.g. angstrom, km, Gpc.
        -cosmo my_cosmology  Set of cosmological parameters, astropy style.
                             Allowed cosmologies are WMAP5, WMAP7, WMAP9,
                             Planck13, and Planck15 (default)

    Examples:
        > python timeline.py 1e-32 s            # Properties just after inflation
        > python timeline.py 13.79 Gyr          # Properties today
        > python timeline.py 500 Myr -Runit Gpc # Properties 500 million years after
                                                # Big Bang, but use Gpc (giga-parsec)
                                                # for distances

    Same examples from Python environment:
        >>> from astropy import units as u
        >>> import timeline
        >>> timeline.uniProp(1e-32*u.s)
        >>> timeline.uniProp(13.79*u.Gyr)
        >>> timeline.uniProp(500*u.Myr, Runit=u.Gpc)
    """


    def photonPressure(T):
        """Photon pressure"""
        return (pi**2 * cc.k_B**4 / (45 * cc.c**3 * cc.h**3) * T**4)


    def nph(T):
        """Photon number density"""
        return 16 * pi * (cc.k_B*T / (cc.h*cc.c))**3 * zeta(3)

    # Calculate radiation-matter equality
    assert (isinstance(t,u.Quantity)) and (t.unit.is_equivalent(u.s)), '\n\nKeyword `t` must have units of time.'
    assert t >= 1e-32*u.s,    '\n\nt must be >= the end of inflation, which is assumed to be at 1e-32 s.'
    assert t <= cosmo.age(0), "\n\nt must be <= the age of the Universe, which for the chosen cosmology is {:}".format(cosmo.age(0))
    a_eq = newton(a_eqSolver,3400.,args=(cosmo,))
    z_eq = 1/a_eq - 1
    t_eq = cosmo.age(z_eq).to(u.yr) # (4/3. * (1 - 1/sqrt(2)) * (cosmo.Onu0 + cosmo.Ogamma0)**1.5/cosmo.Om0**2 / cosmo.H0).to(u.yr)
    rho_eq = cosmo.critical_density(z_eq)

    # Calculate matter-dark energy equality
    a_DE = (cosmo.Om0 / cosmo.Ode0)**.3333333
    z_DE = 1/a_DE - 1
    t_DE = cosmo.age(z_DE).to(u.Gyr)

    # Calculate redshift, density, and temperature
    if t <= t_eq:
        epoch = 'photon epoch'
        a   = (a_eq * sqrt(t/t_eq).decompose()).value
        rho = rho_eq * (a_eq/a)**4
        z   = 1./a - 1
        T   = cosmo.Tcmb0 / a
    else:
        epoch = 'matter epoch' if t<t_DE else 'dark energy epoch'
        from astropy.cosmology import z_at_value
        z   = z_at_value(cosmo.age, t, zmin=0, zmax=5000)
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
    print('Cosmology: ', cosmo)
    print('Redshift and age at radiation-matter equality:   {:.0f} {:.0f}'.format(z_eq, (t_eq).to(u.yr)))
    print('Redshift and age at matter-dark energy equality: {:.2f} {:.2f}'.format(z_DE, (t_DE).to(u.Gyr)))
    print('Radius of Universe today:                        {:2.2f}'.format(R0.to(u.Glyr)))
    print()
    print('Physical properties at t = {:}'.format(t) + ' ('+epoch+'):')
    print(' * Expansion:')
    print('   - Scale factor:                         {:.3g}'.format(a))
    print('   - Redshift:                             {:.4g}'.format(z))
    print('   - Hubble parameter:                     {:.3g}'.format(H))
 #  print('     - H(t) in c / m:                      {:.3g}'.format(H.to(u.c/u.m).value))
 #  print('   - e-foldings since t:                   {:2.1f}'.format(np.log(1/a)))
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

  # a_eq = (cosmo.Onu0 + cosmo.Ogamma0) / cosmo.Om0
    a_eq = newton(a_eqSolver,3400.,args=(cosmo,))
    z_eq = 1/a_eq - 1
    t_eq = cosmo.age(z_eq).to(u.s).value

    def inva(t):
        if t <= t_eq:
          # a = sqrt(2 * sqrt(cosmo.Onu0+cosmo.Ogamma0) * (cosmo.H0 * t).decompose()).value
            a = a_eq * sqrt(t/t_eq)
        else:
            z = z_at_value(cosmo.age, t*u.s, zmax=1e4)
            a = cosmo.scale_factor(z)
        return 1 / a

    eta,err = sp.integrate.quad(inva,0,cosmo.age(z).to(u.s).value) * u.s

    return cc.c * eta / (1+z)


def a_eqSolver(a,cosmo):
    z   = 1./a - 1
    Og0 = cosmo.Ogamma0
    Om0 = cosmo.Om0
    fa  = cosmo.nu_relative_density(z)
    return (1+fa) * Og0/Om0 - a


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('time',                       type=float, help='Time quantity')
    parser.add_argument('unit',                       type=str,   help='Time unit (s,min,day,yr,kyr,Myr,Gyr)')
    parser.add_argument('-Runit', default='',         type=str,   help='Unit for distances (not exactly sure why this works)')
    parser.add_argument('-cosmo', default='Planck15', type=str,   help='Cosmology (WMAP5, WMAP7, WMAP9, Planck13, Planck15). Default is Planck15')
    args = parser.parse_args()

    # Set time unit
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

    # Set distance unit
    if args.Runit != '':
        Runit = args.Runit

    #Set cosmology
    if args.cosmo == 'WMAP5':
        cosmo = WMAP5
    elif args.cosmo == 'WMAP7':
        cosmo = WMAP7
    elif args.cosmo == 'WMAP9':
        cosmo = WMAP9
    elif args.cosmo == 'Planck13':
        cosmo = Planck13
    elif args.cosmo == 'Planck15':
        cosmo = Planck15
    else:
        print("Sorry, cosmology `"+args.cosmo+"` is not implemented.\nYou're welcome to go ahead and do it yourself,\nand then send me a pull request.")
        exit()

    # Calculate it!
    uniProp(t=time, Runit=Runit, cosmo=cosmo)

if __name__ == '__main__':
    main()
