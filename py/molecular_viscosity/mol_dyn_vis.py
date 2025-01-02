def mol_dyn_vis(p,T,S):
    """
    mol_dyn_vis calculates the molecular dynamic viscosity according to
    Richards 1998: The effect of temperature, pressure, and salinity
    on sound attenuation in turbid seawater. J. Acoust. Soc. Am. 103 (1),
    originally published by  Matthaeus, W. (1972): Die Viskositaet des
    Meerwassers. Beitraege zur Meereskunde, Heft 29 (in German).

    returns: molecular dynamic viscosity [kg/(m*s)]

    p:  pressure [dbar]
    T:  temperature [deg C]
    S:  salinity [psu]


    Valid range: 
    0 degC < T < 30 degC, 0 < S < 36, and 1 dbar < p < 1000 dbar.
    """
    # Unit conversion: g / (cm*s) -> kg / (m*s)
    mu = 0.1   *(1.79e-2                                                                           \
                 - 6.1299e-4*T + 1.4467e-5*T**2                                                    \
                 - 1.6826e-7*T**3                                                                  \
                 - 1.8266e-7*p  + 9.8972e-12*p**2                                                  \
                 + 2.4727e-5*S                                                                     \
                 + S*(4.8429e-7*T - 4.7172e-8*T**2                                                 \
                      + 7.5986e-10*T**3)                                                           \
                 + p*(1.3817e-8*T - 2.6363e-10*T**2)                                               \
                 - p**2 * (6.3255e-13*T - 1.2116e-14*T**2))

    # Setting a minimum viscosity value 
    #mu = max(0.00043,mu)
    # this refers to ca. 10 dbar, 45 degC and 0 psu which would be: 0.0004277485639235001 
    # with the formula above
    return mu
