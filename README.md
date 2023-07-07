# Adiabatic Spin Transport

Code based off the work by Jeff Martin found at https://github.com/jmartin1454/asr, just using a different coordinate system.

Calculates the example of adiabatic spin rotation from

Use of Rotating Coordinates in Magnetic Resonance Problems
I. I. Rabi, N. F. Ramsey, and J. Schwinger
Rev. Mod. Phys. 26, 167 (1954).

https://doi.org/10.1103/RevModPhys.26.167

The calculation is based on Equation (12). For a more extensive list of references see the LaTeX document. The exact analytic solution is derived and plotted, as well as a plot of the solution found using Runge-Kutta integration.

----

Initially, the spin points along the +y direction, and so does the magnetic field.  Then, the magnetic field begins to slowly rotate until pointing along the +z direction.  In the adiabatic limit, where the magnetic field vector changes slowly compared to the Larmor frequency, the spin is dragged along with the magnetic field.

<p align="middle">
<img align="middle" src="./reference photos/frames.jpg" width="700" title="Different frames of reference" />  
</p>

For a full explanation of the theory see [Adiabtic_Transport_theory_work.pdf](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/Adiabatic_Transport_theory_work.pdf) or the working LaTeX document at https://www.overleaf.com/read/rfkkhhwwfjct.

## [RotationTransformations.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/RotationTransformations.ipynb)
This notebook containts the analytic work for solving Bloch's equation from a doubly rotating frame of reference where the Bloch equation is trivial to solve. Then using the intial condition of the polarization vector, the full time dependant solution is transformed back into the lab frame.

## [adiabatic-spin-rotation.pdf](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/adiabatic-spin-rotation.pdf)

Derivation by Jeff Martin for solving in only the singly rotating reference frame.   

*Note this work is for a different coordinate set up, rotating from +x.*

## [asr.py](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/asr.py) & [asr_rk.py](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/asr_rk.py)

These files compares the actual spin dynamics with the adiabatic limit.

## [AdiabaticParameterRequirements.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/AdiabaticParameterRequirements.ipynb)
This notebook containts numerical work looking at the minimum requirements for the adiabaticity condition being fulfilled.
