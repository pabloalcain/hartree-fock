# -*- coding: utf-8 -*-
"""
This class will work calculating hartree-fock of single atoms. But we
won't be able to work very cleanly. Unfortunatley --yet again-- the
fotran back-end isn't robust enough to construct the ideal class
set. Even worse than with the previous case, now the errors cannot be
backtraced easily. Eventually this will be rewritten in a more modern
fashion, but until then please report any error you see!
"""

import numpy as np
import johnson
import warnings
from collections import OrderedDict
import sys

this = sys.modules[__name__]
this.orbs = []

_ang = ['s', 'p', 'd', 'f', 'g']


def _nl2orb(n, l):
  """
  Helper function, converts n and l to a readable orbital

  Parameters
  ----------

  n : int
      quantum principal number

  l : int
      angular momentum

  Returns
  -------
  orb : string
      visual representation of the orbital (1s, 2p, etc)
  """
  return str(n)+_ang[l]

def create_atom(z, c=0):
  """
  Create an atom

  Parameters
  ----------

  z : integer
      Charge of the atom

  c : integer, optional
      Index of last core shell. If c = 0, all of the orbitals that are
      created will be considered as core
  """
  johnson.atom(z, c)

def add_orbitals(orbital_list):
  """Add orbitals to the simulation

  Parameters
  ----------

  orbital_list : iterable
      Any iterable that has for entry a tuple (n,
      l, o), with n the prinicipal quantum number, l the angular
      momentum and o the occupancy. if o is 0, the orbital will be
      considered full. For example, if we want to add the orbitals 1s
      and 2s fully occupied, `orbital_list = [(1, 0, 0), (2, 0, 0)]

  Raises
  ------

  Warning
      When the occupancy is higher than the theoretical (2*(2*k + 1))
      or when the angular momentum is higher than the principal
      quantum number. That specific orbital is not added, but the rest
      is.
  """

  for orb in orbital_list:
    if orb[2] > 4*orb[1] + 2:
      msg = ("Occupancy of shell {0} too high for orbital of angular "
             "momentum {1}").format(orb[2], orb[1])
      warnings.warn(msg, UserWarning)
      continue
    if orb[1] > orb[0]:
      msg = ("Angular momentum {0} too high for quantum principal "
             "number {1}").format(orb[1], orb[0])
      warnings.warn(msg, UserWarning)
      continue
    this.orbs.append(orb)
    johnson.add_orbital(*orb)

def set_grid(r0=5e-4, h=0.03):
  """
  Set an exponential interpolation grid, of the form
  r(i) = r0 [exp((i-1)*h) -1]

  Parameters
  ----------
  r0 : float, optional
      first position of the grid

  h : float, optional
      exponential step
  """
  johnson.set_grid(r0, h)

def hartree_fock(damp=0.5):
  """
  Main routine that calculates hartree-fock of the atom

  Parameters
  ----------

  damp : float, optional
      Damping parameter for the calculation of hartree-fock. Ideally,
      any parameter lower than 1 will converge, but usually a sweet
      spot is near 0.5

  Returns
  -------

  etot, r, result: float, array, dict
      etot is the total Hartree-Fock energy, r is the radial grid and
      result is a dictionary of dictionaries, one per orbital. The
      entries of results are, for example, '1s' and '2s'. The
      dictionaries of each orbital, res['1s'] has three keys:
        - 'energy': the energy of the orbital
        - 'p': the wavefunction u
        - 'q': the derivative of u

  Notes
  -----

  `p` and `q` *are not* the radial wavefunctions and derivative, but
  rather

  .. math::
    p(r) = u_{nl}(r) = R_{nl}(r)\,r\\
    q(r) = u'_{nl}(r) = R'_{nl}(r)\,r + R_{nl}

  """
  johnson.hartnr(damp, 0)
  johnson.pickup(0)
  johnson.output(1.0, 0)
  result = {}
  energy = {}
  p = {}
  q = {}
  r = johnson.radial.r[:]
  result = OrderedDict()
  etot = float(johnson.frontend.etot)
  for i, orb in enumerate(this.orbs):
    result[_nl2orb(*orb[:2])] = {}
    d = result[_nl2orb(*orb[:2])]
    d['energy'] = float(johnson.fixdat.wf[i])
    d['p'] = johnson.wavefn.pf[:, i]
    d['q'] = johnson.wavefn.qf[:, i]
  return etot, r, result
