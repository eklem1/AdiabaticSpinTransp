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

For a full explanation of the theory see Chapter 3 of my thesis: [ubc_2023_november_klemets_emma.pdf](https://ucn.triumf.ca/ucn-group-publications/student-theses/ubc_2023_november_klemets_emma.pdf/view).

Previous work (now outdated) can be found at [Adiabtic_Transport_theory_work.pdf](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/Adiabatic_Transport_theory_work.pdf) or the working LaTeX document at https://www.overleaf.com/read/rfkkhhwwfjct.

---

# Files

## [adiabatic-spin-rotation.pdf](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/adiabatic-spin-rotation.pdf)
Derivation by Jeff Martin for solving in only the singly rotating reference frame.   
*Note this work is for a different coordinate set up, rotating from +x.*

## [Adiabatic_Transport_theory_work.pdf](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/Adiabatic_Transport_theory_work.pdf)
The now outdated theory work document, with Derek's edits in [DerekEdit_Adiabatic_Transport_theory_work_df[1394].pdf](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/DerekEdit_Adiabatic_Transport_theory_work_df[1394].pdf)

## [Slichter_Textbook_PrinciplesofMagneticResonance.pdf](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/Slichter_Textbook_PrinciplesofMagneticResonance.pdf)
Useful NMR textbook: Principles of Magnetic Resonance by Springer

## [VV_Vladimirskii_JETP12.pdf](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/VV_Vladimirskii_JETP12.pdf)
The elusive 1960s paper by V.V. Vladimirskii (scan of hardcopy found at the UBC library). This translated version of the russian paper does not seem to be available online, with each of the scanned pages found in [/Vladimirski_paper](https://github.com/eklem1/AdiabaticSpinTransp/tree/master/Vladimirski_paper).

---

# Python Scripts

## [AdiabaticEquations.py](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/AdiabaticEquations.py)
This python notebook contains all the version of the equations used to calculate k, and their (analytic) errors. This file is used extensively in the mapping work done for the [Magnetic Test Enviroment](https://njord.triumf.ca:3000/eklemets/SimModels_work/src/branch/main/MiniEnviro), and explained a fair amount in chapters 4 & 7 in my [thesis](https://ucn.triumf.ca/ucn-group-publications/student-theses/ubc_2023_november_klemets_emma.pdf/view).


## [asr.py](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/asr.py) & [asr_rk.py](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/asr_rk.py)
These files compares the actual spin dynamics with the adiabatic limit, creating a video that plots this comparison out with vectors. Can be used to make a [video](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/asr.mp4) or take a [single image](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/fig72.pdf) from it.

---

# Jupyter notebooks
This contain most of the actual work done here. All figures are set to be saved in either [/Photos](https://github.com/eklem1/AdiabaticSpinTransp/tree/master/Photos) or [/thesisPlots](https://github.com/eklem1/AdiabaticSpinTransp/tree/master/thesisPlots)

## [AdiabaticParameterRequirements.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/AdiabaticParameterRequirements.ipynb)
This notebook containts numerical work looking at the minimum requirements for the adiabaticity condition being fulfilled.

## [AnewMess.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/AnewMess.ipynb)
Just some code testing for plotting, and figuring out constants.

## [CalculatingK_fromMaps.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/CalculatingK_fromMaps.ipynb)
Looking at how $k$ can be calculated from existing magnetic field map data.

## [CalculatingK.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/CalculatingK.ipynb)
Looking at adiabatic transport from the side of the physical parameters, mainly the magnetic field. Makes 2D color maps to show the balance between the magnitude of the fields and their gradients.

## [CollectingFields.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/CollectingFields.ipynb)
My attempt to collect existing data that I can use to get a general idea of the kind of field in the region of the guides, that can then be used to try to figure out how adiabatic the transport would be without any guiding fields in place. Some of this is used in [CalculatingK_fromMaps.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/CalculatingK_fromMaps.ipynb).

## [RotationTransformations.ipynb](https://github.com/eklem1/AdiabaticSpinTransp/blob/master/RotationTransformations.ipynb)
This notebook containts the analytic work for solving Bloch's equation from a doubly rotating frame of reference where the Bloch equation is trivial to solve. Then using the intial condition of the polarization vector, the full time dependant solution is transformed back into the lab frame.
