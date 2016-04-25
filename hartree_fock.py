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
  """
  johnson.hartnr(damp, 0)
  johnson.pickup(0)
  johnson.output(1.47749, 0)
