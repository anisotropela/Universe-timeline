# Calculating the properties of the Universe from Big Bang to now
# --- Documentation for the code `timeline.py`

## Purpose

Calculate various properties of the Universe at a given time _t_.



## History

The code was written in conjunction with the (Danish) popular science article
"Big Bang --- en øjenvidneberetning" ("Big Bang --- an eyewitness account"), in
order to calculate various properties of the Universe. It is meant for being
run from the command line, since this is more tractable to non-experts, but can
also be run from a Python environment is a slightly different way.



## Prerequisites (mostly for non-experts)

##### A terminal
First and foremost, you'll need to use a "terminal" which is a program that
allows you to give commands. On most computer, you'll have an app called
something like "Terminal".

##### A Python installation
The code is written in Python, so you'll need an installation of Python. If you
don't have that, you can get it [here](https://www.python.org) (click
"Downloads" and choose the one matching your computer).  

##### The `astropy` package
In addition to standard Python libraries, you will need the `astropy` library.
When you have Python installed, you can install `astropy` by typing 

    > pip install astropy
(the `>` is just a way of showing that here comes a command; it shouldn't be
included in the command. If you get an error when trying to install, try
writing the word `sudo` in front of the above command, and then type your
computer's password.) 



## Usage from the command line:

In a terminal, type the following command (in the same directory where you put
`timeline.py`): 

    > python timeline.py time unit [-Runit distance_unit] [-cosmo cosmology]

##### Arguments (for command line):

In the above command, there are two mandatory arguments (or "keywords", i.e.
words that must be written): 

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

##### Optional arguments

In addition to `time` and `unit`, you may give two more arguments; one for
outputting the result in your preferred distance units, and one for using your
preferred set of cosmological parameters. The syntax is 

    -Runit my_dist_unit    Units for output distances.
                           Allowed values include:
                           angstrom (or AA), nm, mm, cm, m, km, AU,
                           lightyears (or lyr), parsec (or pc), klyr,
                           kpc, Mlyr, Mpc, Glyr, and Gpc.

    -cosmo my_cosmology    Set of cosmological parameters.
                           Allowed values are:
                           Planck15 (default value), Planck13, WMAP9,
                           WMAP7, WMAP5.

The words `-Runit` and `-cosmo` (don't forget the hyphen "`-`") are written
after the command, followed by your preferred value.  

That is, if you want distances to be written in, say, parsec, you write `-Runit
pc` after your command, and if you prefer a WMAP 2009 cosmology rather than a
Planck 2015 cosmology, you append your command with `-cosmo WMAP9`.



## Examples

(again, the "`>`" in the examples shouldn't be included)

Calculate the properties just after inflation:

    > python timeline.py 1e-32 s

Calculate the properties today:

    > python timeline.py 13.79 Gyr

Calculate the properties 500 million years after Big Bang, but use Gpc
(giga-parsec, i.e. billion parsec) for distances: 

    > python timeline.py 500 Myr -Runit Gpc

Calculate properties a microsecond after Big Bang, with distances written in
cm, using a Planck 2013 cosmology: 

    > python timeline.py 1e-6 s -Runit cm -cosmo Planck13



## Usage from a Python environment

Instead of typing the command directly in the terminal, you can give the
command 

    > python

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
