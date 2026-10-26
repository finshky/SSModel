r'''GW_models_2f.py extends GW_models.py from cosmoGW 
to the 2-fluid case (arXiv xx.yy). It computes
- the double time-integral Delta_mn_2f (Eq. 46)
- the GW spectrum Omega_GW_2f (Eq. 45)
- the single-fluid spectrum Omega_GW, corresponding to
    Omega_GW_2f for gamma=r=0.
    This function implements the same equation as 
    effective_ET_correlator from cosmoGW, with
    some (documented) modifications.
The two parameters introduced in the 2-fluid scenario
are
- gamma = \Gamma_\text{int}
Similarly In addition, 

Author
------
- **Marco Finetti**
  (`marco.finetti@hotmail.com <mailto:marco.finetti@hotmail.com>`_)

Dates
-----
- Created: **26/10/2026**

Contributors
------------
- **Daniil Krichevskii**, **Helena Kolesova**

References
----------
arXiv xx.yy
'''
from cosmoGW import GW_models

def Omega_GW(
    k, EK, tfin,
    k_value=None,
    HR=1.,
    verbose=True,
    Np=Np_ref, Nk=Nkconv_ref, Nz=500, double_Ns=False,
    plot=False, expansion=True,
    extend=False, largek=3, smallk=-3,
    tini=tini_ref, cs2=cs2_ref,
    terms="all", inds_m=None, inds_n=None,
    ):
    r'''Modification of effective_ET_correlator_stat: 
    - saving memory by iterating over k (instead of meshgrid)
    - including prefactor 3 Gamma_ad^2 K^3
    - scaling moment in Delta(t,k,p,p) by HR instead of kstar
    - printing result at each iteration
    Compute normalized GW spectrum :math:`\zeta_{\rm GW}(k)`
    from velocity field spectrum for compressional anisotropic stresses,
    assuming stationary UETC
    (e.g., sound waves under the sound-shell model).

    Parameters
    ----------
    k : array_like
        Wave numbers.
    EK : array_like
        Kinetic spectrum.
    tfin : float
        Final time.
    Np : int, optional
        Number of p discretizations.
    Nk : int, optional
        Number of k discretizations.
    plot : bool, optional
        Plot interpolated spectrum.
    expansion : bool, optional
        Include expansion of Universe.
    kstar : float, optional
        Peak position of the kinetic spectrum.
    extend : bool, optional
        Extend wave number array.
    largek, smallk : int, optional
        Range for extended k (in log10).
    tini : float, optional
        Initial time of GW production.
    cs2 : float, optional
        Speed of sound squared (default 1/3).
    terms : str, optional
        Terms to compute ('all').
    inds_m, inds_n : list, optional
        Indices for Delta_mn.
    corr_Delta_0 : bool, optional
        Correct Delta_0.

    Returns
    -------
    kp : ndarray
        Wave numbers.
    Omm : ndarray
        Normalized GW spectrum, :math:`\zeta_{\rm GW}(k)`.

    References
    ----------
    RoperPol:2023dzg, eq. 93
    '''
    Gamma_ad = 1+cs2            # adiabatic index
    pref     = 3*Gamma_ad**2    # prefactor of Omega_GW
    
    k = np.asarray(k)
    EK = np.asarray(EK)

    p = np.logspace(np.log10(k[0]), np.log10(k[-1]), Np)
    kp = np.logspace(np.log10(k[0]), np.log10(k[-1]), Nk)

    if extend:
        kp = np.logspace(smallk, largek, Nk)
        # k1, k2, k3, k4 = np.sort([smallk,np.log10(k[0]),np.log10(k[-1]),largek])
        # Nk_ext = int(Nk / 6)
        # kp = np.concatenate([
        #     np.logspace(k1, k2, Nk_ext),
        #     np.logspace(k2, k3, 4 * Nk_ext),
        #     np.logspace(k3, k4, Nk_ext),
        #     ])

    if k_value:
        kp = [k_value]

    z = np.linspace(-1, 1, Nz)

    EK_p = np.interp(p, k, EK)

    if plot:
        plt.plot(p, EK_p)
        plt.xscale("log")
        plt.yscale("log")

    if terms == "all":
        pairs = [(m, n) for m in (-1, 1) for n in (-1, 1)]
    else:
        inds_m = [] if inds_m is None else inds_m
        inds_n = [] if inds_n is None else inds_n
        pairs = [(m, n) for m in inds_m for n in inds_n]

    l   = len(pairs)
    Omm = np.zeros((l + 1, len(kp)))

    # Reusable 2D broadcasting arrays
    p2d = p[:, None]
    z2d = z[None, :]
    angular = (1 - z2d**2)**2
    pweight = p**2 * EK_p
    
    doubled = False
    if verbose: print(f"{'K':>10}{'Omega_GW':>10}")
    for ik, kk in enumerate(kp):
        if double_Ns and kk>10. and not doubled:
            z    = np.linspace(-1, 1, 2*Nz)
            p    = np.logspace(np.log10(k[0]), np.log10(k[-1]), 2*Np)
            EK_p = np.interp(p, k, EK)
            # Reusable 2D broadcasting arrays
            p2d = p[:, None]
            z2d = z[None, :]
            angular = (1 - z2d**2)**2
            pweight = p**2 * EK_p
            doubled = True
            
        # Arrays have shape (Np, Nz), not (Nk, Np, Nz)
        ptilde = np.sqrt(
            p2d**2 + kk**2 - 2 * kk * p2d * z2d
        )

        EK_ptilde = np.interp(ptilde, k, EK)
        np.maximum(ptilde, 1e-50, out=ptilde)

        # Shared integrand, independent of m and n
        weight = EK_ptilde / ptilde**4
        weight *= angular

        for j, (m, n) in enumerate(pairs):

            Delta = compute_Delta_mn(
                tfin,
                kk/HR,
                p2d/HR,
                ptilde/HR,
                cs2=cs2,
                m=m,
                n=n,
                tini=tini,
                expansion=expansion,
            )

            # Same normalization as the original implementation
            Delta /= l

            # Integrate over z
            Pi_1 = safe_trapezoid(
                weight * Delta, z, axis=1
            )

            # Integrate over p
            Omm[j, ik] = pref*kk**3*safe_trapezoid(
                Pi_1 * pweight, p
            )
        Omm[-1,ik] = np.sum(Omm[:-1,ik])
        if verbose: print(f"{kk:>10.2e}{Omm[-1, ik]:>10.2e}")

    return kp, Omm


def Omega_GW_2f(
    k, EK, tfin,
    gamma,r,
    HR=1.,
    k_value=None,
    verbose=True,
    Np=Np_ref, Nk=Nkconv_ref, Nz=500, double_Ns=False,
    plot=False, expansion=True,
    extend=False, largek=3, smallk=-3,
    tini=tini_ref, cs2=cs2_ref,
    terms="all", inds_m=None, inds_n=None,
    ):
    r'''gamma,r-depenendent gravitational wave spectrum Omega_GW for 2 fluids.

    Compute normalized GW spectrum :math:`\zeta_{\rm GW}(k)`
    from velocity field spectrum for compressional anisotropic stresses,
    assuming stationary UETC
    (e.g., sound waves under the sound-shell model).

    Parameters
    ----------
    k : array_like
        Wave numbers.
    EK : array_like
        Kinetic spectrum.
    tfin : float
        Final time.
    Np : int, optional
        Number of p discretizations.
    Nk : int, optional
        Number of k discretizations.
    plot : bool, optional
        Plot interpolated spectrum.
    expansion : bool, optional
        Include expansion of Universe.
    HR : float, optional
        Hubble parameter*bubble radius at t_*.
    extend : bool, optional
        Extend wave number array.
    largek, smallk : int, optional
        Range for extended k (in log10).
    tini : float, optional
        Initial time of GW production.
    cs2 : float, optional
        Speed of sound squared (default 1/3).
    terms : str, optional
        Terms to compute ('all').
    inds_m, inds_n : list, optional
        Indices for Delta_mn.
    corr_Delta_0 : bool, optional
        Correct Delta_0.

    Returns
    -------
    kp : ndarray
        Wave numbers.
    Omm : ndarray
        Normalized GW spectrum, :math:`\zeta_{\rm GW}(k)`.

    References
    ----------
    RoperPol:2023dzg, eq. 93
    '''
    Gamma_ad = 1+cs2            # adiabatic index
    pref     = 3*Gamma_ad**2    # prefactor of Omega_GW
    
    k = np.asarray(k)
    EK = np.asarray(EK)

    p = np.logspace(np.log10(k[0]), np.log10(k[-1]), Np)
    kp = np.logspace(np.log10(k[0]), np.log10(k[-1]), Nk)

    if extend:
        kp = np.logspace(smallk, largek, Nk)
        # k1, k2, k3, k4 = np.sort([smallk,np.log10(k[0]),np.log10(k[-1]),largek])
        # Nk_ext = int(Nk / 6)
        # kp = np.concatenate([
        #     np.logspace(k1, k2, Nk_ext),
        #     np.logspace(k2, k3, 4 * Nk_ext),
        #     np.logspace(k3, k4, Nk_ext),
        #     ])

    if k_value:
        kp = [k_value]
        
    z = np.linspace(-1, 1, Nz)

    EK_p = np.interp(p, k, EK)

    if plot:
        plt.plot(p, EK_p)
        plt.xscale("log")
        plt.yscale("log")

    if terms == "all":
        pairs = [(m, n) for m in (-1, 1) for n in (-1, 1)]
    else:
        inds_m = [] if inds_m is None else inds_m
        inds_n = [] if inds_n is None else inds_n
        pairs = [(m, n) for m in inds_m for n in inds_n]

    l   = len(pairs)
    Omm = np.zeros((l + 1, len(kp)))

    # Reusable 2D broadcasting arrays
    p2d = p[:, None]
    z2d = z[None, :]
    angular = (1 - z2d**2)**2
    pweight = p**2 * EK_p

    doubled = False
    if verbose: print(f"{'K':>10}{'Omega_GW':>10}")
    for ik, kk in enumerate(kp):
        if double_Ns and kk>10. and not doubled:
            z    = np.linspace(-1, 1, 2*Nz)
            p    = np.logspace(np.log10(k[0]), np.log10(k[-1]), 2*Np)
            EK_p = np.interp(p, k, EK)
            # Reusable 2D broadcasting arrays
            p2d = p[:, None]
            z2d = z[None, :]
            angular = (1 - z2d**2)**2
            pweight = p**2 * EK_p
            doubled = True

        # Arrays have shape (Np, Nz), not (Nk, Np, Nz)
        ptilde = np.sqrt(
            p2d**2 + kk**2 - 2 * kk * p2d * z2d
        )

        EK_ptilde = np.interp(ptilde, k, EK)
        np.maximum(ptilde, 1e-50, out=ptilde)

        # Shared integrand, independent of m and n
        weight = EK_ptilde / ptilde**4
        weight *= angular

        for j, (m, n) in enumerate(pairs):

            Delta = compute_Delta_mn_2f(
                tfin,
                kk/HR,
                p2d/HR,
                ptilde/HR,
                gamma,r,   # additional gamma, r dependence
                cs2=cs2,
                m=m,
                n=n,
                tini=tini,
                expansion=expansion,
            )

            # Same normalization as the original implementation
            Delta /= l

            # Integrate over z
            Pi_1 = safe_trapezoid(
                weight * Delta, z, axis=1
            )

            # Integrate over p
            Omm[j, ik] = pref*kk**3*safe_trapezoid(
                Pi_1 * pweight, p
                )
        Omm[-1,ik] = np.sum(Omm[:-1,ik])
        if verbose: print(f"{kk:>10.2e}{Omm[-1, ik]:>10.2e}")
        if k_value: return kp[0], Omm[:,0]

    return kp, Omm

# 2-fluid extension of compute_Delta_mn
def compute_Delta_mn_2f(
    t, k, p, ptilde, gamma, r,
    cs2=1/3, m=1, n=1,
    tini=1.0, expansion=True
):
    cs = np.sqrt(cs2)

    # Original Python convention for the frequencies
    q = n*k + cs*(m*ptilde + p)
    q = np.asarray(q, dtype=float)

    if not expansion:
        raise NotImplementedError(
            "The Mathematica comparison assumes expansion=True."
        )

    # Undamped integral: Ci and Si differences
    qsafe = np.where(q == 0, 1.0, q)

    Si_t, Ci_t = special.sici(qsafe*t)
    Si_ini, Ci_ini = special.sici(qsafe*tini)

    DCi = np.where(q == 0, np.log(t/tini),
                   Ci_t - Ci_ini)
    DSi = np.where(q == 0, 0.0,
                   Si_t - Si_ini)

    I0 = DCi + 1j*DSi

    # Damped integral: use E1 rather than Ei
    a = 2*gamma - 1j*q

    if gamma == 0:
        Igamma = I0
    else:
        Igamma = np.exp(2*gamma) * (
            special.exp1(a*tini)
            - special.exp1(a*t)
        )

    return np.abs(I0 + r*Igamma)**2 / (1+r)**2