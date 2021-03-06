# Toronto2015
___

## Project
The Kepler satellite has found several compact solar systems containing earth mass planets. Yet, these systems do not appear to conatain any giant planets. This purpose of this project is to investigate whether or not this lack of giant planets is due to the bias in Kepler to find planets closer to the star or if we don't find them because their presence would destabilize the system. 

Destabilization due to an MMR has been ruled out because for a giant planet to be in a 2:1 resonance with the outermost planet would require it to have a semi-major axis of ~0.7 AU at which point it would have been detected by Kepler. Therefore, the most likely cause for destabilization would be a secular perturbation due to two giant planets further out, ~5 AU, with an eigenfrequency resonant with an eigenfrequency of one of the inner planets.

The majority of the code uses [Rebound](https://github.com/hannorein/rebound) which is an open-source multi-purpose N-body code developed by:

* Hanno Rein, University of Toronto, <hanno@hanno-rein.de>
* Shangfei Liu, Kavli Institute for Astronomy and Astrophysics at Peking University (KIAA-PKU), Beijing, <liushangfei@pku.edu.cn>
* David S. Spiegel, Institute for Advanced Study (IAS), Princeton, <dave@ias.edu>
* Akihiko Fujii, National Astronomical Observatory of Japan/University of Tokyo, Tokyo, <akihiko.fujii@nao.ac.jp>
* Dan Tamayo, University of Toronto, <dtamayo@cita.utoronto.ca>
