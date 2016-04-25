import hartree_fock as hf
import pylab as pl

orbitals = [(1, 0, 0), (2, 0, 0), (2, 1, 0), (3, 0, 0),
            (3, 1, 0), (3, 2, 0), (4, 0, 0), (4, 1, 0),]

Z = 36
name = 'Kr'
hf.create_atom(Z, 0)
hf.add_orbitals(orbitals)
hf.set_grid()
etot, r, res = hf.hartree_fock()
print 'Finished calculation for Z = {0} ({1})'.format(Z, name)
print 'Hartree-Fock energy: {0}'.format(etot)

fig, ax = pl.subplots()
# This plots every orbital
for name, orb in zip(res.keys(), res.values()):
  lbl = '{0}, E = {1}'.format(name, orb['energy'])
  ax.plot(r, orb['p']/r, label=lbl)
pl.legend()
ax.set_xlim(0, 1)

#If we want to plot only the 1s, 2s and 2p

fig, ax = pl.subplots()
orbs = ['1s', '2s', '2p']
for orb in orbs:
  info = res[orb]
  lbl = '{0}, E = {1}'.format(orb, info['energy'])
  ax.plot(r, info['p']/r, label=lbl)
pl.legend()
ax.set_xlim(0, 1)

#Or a plot of all l=0
fig, ax = pl.subplots()
# This plots every orbital
for name, orb in zip(res.keys(), res.values()):
  if name[1] != 's': continue
  lbl = '{0}, E = {1}'.format(name, orb['energy'])
  ax.plot(r, orb['p']/r, label=lbl)
pl.legend()
ax.set_xlim(0, 1)
pl.show()
