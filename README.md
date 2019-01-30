# Properties of the Universe from Big Bang to now &mdash; Documentation for the cosmological calculator `timeline`
  
## Purpose

Calculate various properties of the Universe at a given time *t*.

The advantage of `timeline` over other cosmological calculators on the internet (such as
[Ned Wright's famous CosmoCalc](http://www.astro.ucla.edu/~wright/CosmoCalc.html),
[it's Python wrapper](http://cxc.harvard.edu/contrib/cosmocalc/),
[CC](https://home.fnal.gov/~gnedin/cc/),
[WolframAlpha](https://www.wolframalpha.com),
[Astropy](http://www.astropy.org),
[CosmoTools](http://www.bo.astro.it/~cappi/cosmotools),
[cosmo](http://faraday.uwyo.edu/~chip/misc/Cosmo2/cosmo.cgi),
[cosmocalc_2013](http://www.einsteins-theory-of-relativity-4engineers.com/cosmocalc_2013.htm),
[cosmo_calc](http://www.oa-roma.inaf.it/zappacosta/cosmo_calc.html), and
[cosmo_calc](http://srini.ph.unimelb.edu.au/cosmo_calc.php))
is that it goes all the way back to inflation, 0.000...[31 zeros]..1 seconds after Big Bang.

Furthermore, where all these calculators calculate the properties for an input *redshift*, `timeline` calculates for an input *age of the Universe*, which is more intuitive for non-astronomers.

At late epochs, `timeline` makes use of Python's [Astropy](http://www.astropy.org) library, but for the early epochs it "calculates backward" from radation-matter equality, assuming a radiation-dominated Universe during that time.


## Background

The code was written in conjunction with the (Danish) popular science article "[Big Bang &mdash; en øjenvidneberetning](https://videnskab.dk/naturvidenskab/big-bang-en-oejenvidneberetning)" (translated to English: "[Big Bang &mdash; an eyewitness account](http://sciencenordic.com/big-bang-–-eyewitness-account)", and awarded best research outreach 2018 in Denmark at [ForskerZonen](https://videnskab.dk/forskerzonen)), in order to calculate various properties of the Universe. It is meant for being run from the command line, since this is more tractable to non-experts, but can also be run from a Python environment is a slightly different way.

## Prerequisites (mostly for non-experts)

### A terminal
First and foremost, you'll need to use a "terminal" which is a program that allows you to give commands to your computer. On most computers, you'll have an app called something like "Terminal".

### A Python installation
The code is written in Python, so you'll need an installation of Python. If you don't have that, you can get it [here](https://www.python.org) (click "Downloads" and choose the one matching your computer).

### The `astropy` package
In addition to the standard Python installation, you will need the `astropy` library.  When you have Python installed, you can install `astropy` by typing, in the terminal

```
$ pip install astropy
```
(the `$` is just a way of showing that here comes a command; it shouldn't be included in the command. If you get an error when trying to install, try writing the word `sudo` in front of the above command, and then type your computer's password.)

<!-- ## Usage from the command line: -->

## Usage

In a terminal, type the following command (in the same directory where you put `timeline.py`):

    $ python timeline.py time unit [-Runit distance_unit] [-cosmo cosmology]

### Arguments (for command line)

In the above command, `python` is the command to make Python run the program, `timeline.py` is the name of the progam. Additionally, there are two mandatory arguments (i.e. words that *must* be written):

    time    Time quantity, i.e. a number

    unit    Unit of time. Allowed values are
            s:      Seconds
            min:    Minutes
            h:      Hours
            day:    Days
            yr:     Years
            kyr:    kilo-years (i.e. 1000 years)
            Myr:    mega-years
            Gyr:    giga-years

### Optional arguments

Optional arguments are words that *may* be written. There are two; one for outputting the result in your preferred distance units, and one for using your preferred set of cosmological parameters. The syntax is

    -Runit my_dist_unit    Units for output distances.
                           Allowed values include:
                           angstrom (or AA), nm, mm, cm, m, km, AU,
                           lightyears (or lyr), parsec (or pc), klyr,
                           kpc, Mlyr, Mpc, Glyr, and Gpc.

    -cosmo my_cosmology    Set of cosmological parameters.
                           Allowed values are:
                           Planck15 (default value), Planck13, WMAP9,
                           WMAP7, WMAP5.

The words `-Runit` and `-cosmo` (don't forget the dash "`-`") are written after the command, followed by your preferred value.

That is, if you want distances to be written in, say, parsec, you write `-Runit pc` after your command, and if you prefer a WMAP 2009 cosmology rather than a Planck 2015 cosmology, you append your command with `-cosmo WMAP9`.



## Examples

(again, the "`$`" in the examples shouldn't be included)

Calculate the properties just after inflation:

    $ python timeline.py 1e-32 s

Calculate the properties today:

    $ python timeline.py 13.79 Gyr

Calculate the properties 500 million years after Big Bang, but use Gpc
(giga-parsec, i.e. billion parsec) for distances:

    $ python timeline.py 500 Myr -Runit Gpc

Calculate properties a microsecond after Big Bang, with distances written in
cm, using a Planck 2013 cosmology:

    $ python timeline.py 1e-6 s -Runit cm -cosmo Planck13

<!--

## Usage from a Python environment

Instead of typing the command directly in the terminal, you can give the
command 

    $ python

to start a "Python interpreter". Now you're "inside" Python. Here you should
first give the commands 

    >>> from astropy import units as u
    >>> import timeline

where now the "`>>>`" shouldn't be included.

The general syntax for running the code is

    >>> timeline.uniProp(time*unit [,Runit=my_dist_unit] [cosmo=my_cosmology])

and the four examples above can then be executed with

`>>> timeline.uniProp(1e-32*u.s)`

`>>> timeline.uniProp(13.79*u.Gyr)`

`>>> timeline.uniProp(500*u.Myr, Runit=u.Gpc)`

and

    >>> from astropy.cosmology import Planck15
    >>> timeline.uniProp(1e-6*u.s, Runit=u.cm, cosmo=Planck13)`

-->


## Output

The following values are written out for the Universe at the chosen time *t*

* Expansion:
  - Scale factor (*a*, size of the Universe relative to today)
  - Redshift (*z*, how much light emitted from a source is "stretched" before it reaches us)
  - Hubble parameter (*H(t)*, expansion rate of the Universe)
* Size:
  - Radius of observable Universe at *t* (*d*<sub>P</sub>; how far away could an
    observer at *t* theoretically see. This is called *the particle horizon*,
    and this calculation involves an integral that takes several seconds for
    late epochs, so if you're impatient, you may want to delete this line from
    the code)
  - Radius of today's obs. Universe at t (*a d*<sub>P,0</sub>; how big was the
    part of the Universe that we can can see *today* at that time)
  - Hubble distance (*c* / *H(t)*; distance at which the expansion makes stuff
    recede faster than the speed of light)
* Gas and radiation:
  - Temperature (*T*; average temperature of stuff in the Universe)
  - Energy (*E* = *k*<sub>B</sub>*T*; the corresponding energy of particles)
  - Energy density (the total, average energy of atoms, radiation, and everything
    else per volume)
  - Ionized fraction (*x<sub>e</sub>*; fraction of hydrogen atoms that are ionized)
  - Photon mean free path (how far can a photon travel before it hits an electron): mfp / *d*<sub>P</sub> (if this ratio is < 1, radation is coupled to matter;
      if it is > 1, photons free-stream through the entire Universe)
  - Photon no. density (Number of photons per cubic centimeter)
  - Baryon no. density (Number of atoms per cubic centimeter)
  - Photon pressure
  - Baryon pressure
