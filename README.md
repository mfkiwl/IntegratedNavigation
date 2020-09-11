# Integrated Navigation
integrated navigation using INS/GNSS

Due to I know little about computer copyright issue, this repository is not publicly available.

# INS/GNSS integrated navigation uses Openblas(http://www.openblas.net/), I really appreciate.

# How to use?
please see main.cc in src, main.cc is a example how to read data and how to set some arguments.

Body coordinate: front-right-down

navigation coordinate: north-east-down


ins.h/ins.cc have implemented Allan Variance and coarse alignment. Allan Variance is not used in integrated navigation.

# reference
Shin, Eun-Hwan, and Naser El-Sheimy. "Accuracy improvement of low cost INS/GPS for land applications." Proceedings of the 2002 national technical meeting of the institute of navigation. 2002.

Shin, Eun-Hwan. "Estimation techniques for low-cost inertial navigation." UCGE report 20219 (2005).

Quinchia, Alex G., et al. "A comparison between different error modeling of MEMS applied to GPS/INS integrated systems." Sensors 13.8 (2013): 9549-9588.

