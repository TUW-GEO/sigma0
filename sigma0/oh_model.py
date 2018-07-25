# Copyright (c) 2018, TU Wien, Department of Geodesy and Geoinformation (GEO).
# All rights reserved.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import cmath
import numpy as np


def gammah(eps, theta):
    """
    Incoherent reflectivity for H-pol.

    Parameters
    ----------
    eps : complex number
        Complex dielectric constant.
    theta : float
        Angle in degrees.

    Returns
    -------
    gammah : float
        Incoherent reflectivity for H-pol.
    """

    theta_rad = theta / 180.0 * cmath.pi
    costheta = np.cos(theta_rad)
    neliojuuri = np.sqrt(eps - np.sin(theta_rad) ** 2)
    gammah = np.abs((costheta - neliojuuri) / (costheta + neliojuuri)) ** 2

    return gammah


def gammav(eps, theta):
    """
    Incoherent reflectivity for V-pol.

    Parameters
    ----------
    eps : complex number
        Complex dielectric constant.
    theta : float
        Angle in degrees.

    Returns
    -------
    gammav : float
        Incoherent reflectivity for V-pol.
    """
    theta_rad = theta / 180.0 * cmath.pi
    costheta = np.cos(theta_rad)
    neliojuuri = np.sqrt(eps - np.sin(theta_rad) ** 2)
    gammav = np.abs((eps * costheta - neliojuuri) /
                    (eps * costheta + neliojuuri)) ** 2

    return gammav


def sigma0_bare(theta, eps_low, f_rms, f, eps_top=1):
    """
    Oh et.al. (1992) surface backscatter calculations

    This functions calculations surface backscatter using the
    Oh et al. (1992) surface model.

    References
    Oh et al., 1992, An empirical model and an inversion technique for
    rader scattering from bare soil surfaces. IEEE Trans. Geos. Rem.,
    30, pp. 370-380

    Parameters
    ----------
    theta : float
        incidence angle in degrees
    eps_low : complex number
        complex permittivity of lower medium
    f : float
        frequency in hertz
    f_rms: float
        fractional rms height. rms_height in meters divided by wavelength
    eps_top : complex number, optional
        complex permittivity of upper(incoming) medium

    Returns
    -------
    sigma0 : float
        Surface backscatter.
    """
    k_rms = (2. * cmath.pi) * f_rms

    eps_eff = eps_low / eps_top

    gam0 = gammah(eps_eff, 0)
    gamh = gammah(eps_eff, theta)
    gamv = gammav(eps_eff, theta)
    theta = theta / 180.0 * cmath.pi

    # precalulcate cosines of angles
    ct = np.cos(theta)

    # Oh model equations
    g = 0.7 * (1. - np.exp(-0.65 * k_rms ** 1.8))
    root_p = 1. - ((2. * theta / cmath.pi) **
                   (1. / (3. * gam0))) * np.exp(-k_rms)

    # p = (1 - (2 * theta / cmath.pi) ** (1 / (3 * gam0)) * np.exp(-1 * k_rms)) ** 2

    q = 0.23 * np.sqrt(gam0) * (1. - np.exp(-k_rms))

    sigma0 = {}

    sigma0['vv'] = g * (ct * (ct * ct)) * (gamv + gamh) / root_p

    # sig0vv = ((g*(cos(theta))^3)/sqrt(p))*(gamv+gamh);

    sigma0['vh'] = q * sigma0['vv']
    sigma0['hh'] = root_p * root_p * sigma0['vv']

    return sigma0
