
# kinematic utilities
#
# davide.gerbaudo@gmail.com
# Jan 2014

import math

def phi_mpi_pi(phi) :
    pi = math.pi
    while phi < -pi : phi += pi
    while phi > +pi : phi -= pi
    return phi
