# Using the Hartree-Fock Method to Find the Ground State Energy of Helium

## Background
Quantum-Mechanical systems are naturally complicated. For systems bigger than hydrogen, the particle-particle dependence of the potential energies makes the Schrodinger equation unsolvable. The Hartree-Fock approximation neglects some of the inter-particle dependence, which simplfies the equation. If we apply to the Hartree-Fock approximation to the helium atom, we can solve the eigenvalue problem for the ground state using the self-consistent field method.

## Using the script
The script will output successive approximations for the ground state energy until two successive results are within *epsilon* of each other. *Epsilon* is currently set at 1E-8. 

## Graphics
Next the script will output 2 graphs. The first shows the radial probability distribution of the electrons around the nucleus. The second graphs the radial probability distributions at successive steps of the self-consistent field method. The number of x data points is currently set at 100. This can be increased (e.g. to  ) for a higher resolution graph.

## Sources
This project was done as part of an independent study research class at Lee University. Much help was received from mentor Dr. David Pigg and student collaborator Rebekah Petrosky. Parameterization for wavefunction was taken from Jos Thijssen's book *Computational Physics* ( 2007).
