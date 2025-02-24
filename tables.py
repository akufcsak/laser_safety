import numpy as np


def table3(lambda_, t):
    """
    Accessible emission limits for Class 1 and Class 1M laser products and C6 = 1
    :param lambda_: wavelength [nm]
    :param t: emission duration [s]
    :return: AEL_single [J], AEL_single [W]
    """
    # get the correction factors
    c4 = get_c4(lambda_)
    c7 = get_c7(lambda_)
    # set the outputs to infinity in case we need to
    # evaluate AEL at a border value of t
    ael_j = np.inf
    ael_w = np.inf
    # set minimum emission duration
    if t < 1e-13:
        t = 1e-13   # as suggested by table note b
    if 3e4 < t:
        raise ValueError(["No AEL is defined "
                          "for the specified emission duration ({t} s)".format(t=t)])
    # implemented only for 700-1050 nm for now...
    if 700 <= lambda_ <= 1050:
        if 1e-13 <= t <= 1e-11:
            ael_j = np.min([3.8e-8, ael_j])
            ael_w = np.min([ael_j / t, ael_w])
        if 1e-11 <= t <= 5e-6:
            ael_j = np.min([7.7e-8 * c4, ael_j])
            ael_w = np.min([ael_j / t, ael_w])
        if 5e-6 <= t <= 10:
            ael_j = np.min([7e-4 * c4 * t**0.75, ael_j])
            ael_w = np.min([ael_j / t, ael_w])
        if 10 <= t <= 3e4:
            ael_w = np.min([3.9e-4 * c4 * c7, ael_w])
            ael_j = np.min([ael_w * t, ael_j])

    return ael_j, ael_w


def table4(lambda_, t, alpha):
    """
    Accessible emission limits for Class 1 and Class 1M laser products
    in the wavelength range from 400 nm to 1 400 nm (retinal hazard region):
    extended sources
    :param lambda_: wavelength [nm]
    :param t: emission duration [s]
    :param alpha: angular subtense of the apparent source [rad]
    :return: AEL_single [J], AEL_single [W]
    """
    # get the correction factors
    c4 = get_c4(lambda_)
    c6 = get_c6(lambda_, alpha, t)
    t2 = get_t2(lambda_, alpha)
    # set the outputs to infinity in case we need to
    # evaluate AEL at a border value of t
    ael_j = np.inf
    ael_w = np.inf
    # set minimum emission duration
    if t < 1e-13:
        t = 1e-13   # as suggested by table note b
    if 3e4 < t:
        raise ValueError(["No AEL is defined "
                          "for the specified emission duration ({t} s)".format(t=t)])
    # implemented only for 700-1050 nm for now...
    if 700 <= lambda_ <= 1050:
        if 1e-13 <= t <= 1e-11:
            ael_j = np.min([3.8e-8 * c6, ael_j])
            ael_w = np.min([ael_j / t, ael_w])
        if 1e-11 <= t <= 5e-6:
            ael_j = np.min([7.7e-8 * c4 * c6, ael_j])
            ael_w = np.min([ael_j / t, ael_w])
        if 5e-6 <= t <= 10:
            ael_j = np.min([7e-4 * c4 * c6 * t**0.75, ael_j])
            ael_w = np.min([ael_j / t, ael_w])
        if 10 <= t <= 3e4:
            if t <= t2:
                ael_j = np.min([7e-4 * c4 * c6 * t**0.75, ael_j])
                ael_w = np.min([ael_j / t, ael_w])
            else:  # t > t2
                ael_w = np.min([7e-4 * c4 * c6 * t2**-0.25, ael_w])
                ael_j = np.min([ael_w * t, ael_j])

    return ael_j, ael_w


def table_a1(lambda_, t):
    """
    Maximum permissible exposure (MPE) for C6 = 1 at the cornea expressed as irradiance or radiant exposure
    :param lambda_: wavelength [nm]
    :param t: emission duration [s]
    :return: MPE [J/m2], MPE [W/m2]
    """
    # get the correction factors
    c4 = get_c4(lambda_)
    c7 = get_c7(lambda_)
    # set minimum emission duration
    if t < 1e-13:
        t = 1e-13   # as suggested by table note b
    if 3e4 < t:
        raise ValueError(["No MPE is defined "
                          "for the specified emission duration ({t} s)".format(t=t)])
    # set the outputs to infinity in case we need to
    # evaluate AEL at a border value of t
    mpe_re = np.inf     # MPE as radiant exposure
    mpe_i = np.inf      # MPE as irradiance
    # implemented only for 700-1050 nm for now...
    if 700 <= lambda_ <= 1050:
        if 1e-13 <= t <= 1e-11:
            mpe_re = np.min([1e-3, mpe_re])
            mpe_i = np.min([mpe_re / t, mpe_i])
        if 1e-11 <= t <= 5e-6:
            mpe_re = np.min([2e-3 * c4, mpe_re])
            mpe_i = np.min([mpe_re / t, mpe_i])
        if 5e-6 <= t <= 10:
            mpe_re = np.min([18 * c4 * t ** 0.75, mpe_re])
            mpe_i = np.min([mpe_re / t, mpe_i])
        if 10 <= t <= 3e4:
            mpe_i = np.min([10 * c4 * c7, mpe_i])
            mpe_re = np.min([mpe_i * t, mpe_re])

    return mpe_re, mpe_i


def table_a2(lambda_, t, alpha):
    """
    Maximum permissible exposure (MPE) at the cornea for extended sources
    in the wavelength range from 400 nm to 1 400 nm (retinal hazard region)
    expressed as irradiance or radiant exposure
    :param lambda_: wavelength [nm]
    :param t: emission duration [s]
    :param alpha: angular subtense of the apparent source [rad]
    :return: MPE [J/m2], MPE [W/m2]
    """
    # get the correction factors
    c4 = get_c4(lambda_)
    c6 = get_c6(lambda_, alpha, t)
    t2 = get_t2(lambda_, alpha)
    # set minimum emission duration
    if t < 1e-13:
        t = 1e-13  # as suggested by table note b
    if 3e4 < t:
        raise ValueError(["No MPE is defined "
                          "for the specified emission duration ({t} s)".format(t=t)])
    # set the outputs to infinity in case we need to
    # evaluate AEL at a border value of t
    mpe_re = np.inf  # MPE as radiant exposure
    mpe_i = np.inf  # MPE as irradiance
    # implemented only for 700-1050 nm for now...
    if 700 <= lambda_ <= 1050:
        if 1e-13 <= t <= 1e-11:
            mpe_re = np.min([1e-3 * c6, mpe_re])
            mpe_i = np.min([mpe_re / t, mpe_i])
        if 1e-11 <= t <= 5e-6:
            mpe_re = np.min([2e-3 * c4 * c6, mpe_re])
            mpe_i = np.min([mpe_re / t, mpe_i])
        if 5e-6 <= t <= 10:
            mpe_re = np.min([18 * c4 * c6 * t ** 0.75, mpe_re])
            mpe_i = np.min([mpe_re / t, mpe_i])
        if 10 <= t <= 3e4:
            if t <= t2:
                mpe_re = np.min([18 * c4 * c6 * t ** 0.75, mpe_re])
                mpe_i = np.min([mpe_re / t, mpe_i])
            else:  # t > t2
                mpe_i = np.min([18 * c4 * c6 * t2 ** -0.25, mpe_i])
                mpe_re = np.min([mpe_i * t, mpe_re])

    return mpe_re, mpe_i


def table_a5(lambda_, t):
    """
    Maximum permissible exposure (MPE) of the skin to laser radiation
    :param lambda_: wavelength [nm]
    :param t: emission duation [s]
    :return: MPE [J/m2], MPE [W/m2]
    """
    # get the correction factors
    c4 = get_c4(lambda_)
    # set the outputs to infinity in case we need to
    # evaluate AEL at a border value of t
    mpe_re = np.inf  # MPE as radiant exposure
    mpe_i = np.inf  # MPE as irradiance
    if 3e4 < t:
        raise ValueError(["No MPE is defined "
                          "for the specified emission duration ({t} s)".format(t=t)])
    # implemented only for 700-1400 nm for now...
    if t <= 1e-9:
        mpe_i = np.min([2e11 * c4, mpe_i])
        mpe_re = np.min([mpe_i * t, mpe_re])
    if 1e-9 <= t <= 1e-7:
        mpe_re = np.min([200 * c4, mpe_re])
        mpe_i = np.min([mpe_re / t, mpe_i])
    if 1e-7 <= t <= 10:
        mpe_re = np.min([1.1e4 * c4 * t ** 0.25, mpe_re])
        mpe_i = np.min([mpe_re / t, mpe_i])
    if 10 <= t <= 3e4:
        mpe_i = np.min([2000 * c4, mpe_i])
        mpe_re = np.min([mpe_i * t, mpe_re])

    return mpe_re, mpe_i


def get_alpha_max(t):
    """
    :param t: emission duration considered [s] (see Table 9 in standard)
    :return: The value of angular subtense of the apparent source
                above which the MPEs and AELs are independent of
                the source size (alpha_max)
    """
    # alpha_max depends on the emission duration
    if t < 625e-6:
        alpha_max = 5e-3
    elif 625e-6 <= t <= 0.25:
        alpha_max = 200 * t ** 0.5 * 1e-3
    else:  # 0.25 < t
        alpha_max = 100e-3
    return alpha_max


def get_t1(lambda_):
    """
    :param lambda_: wavelength [nm]
    :return: breakpoint T1 [s]
    """
    # t1
    if 302.5 <= lambda_ <= 315:
        return 10**(0.8 * (lambda_-295)) * 1e-15
    else:
        raise ValueError(["No breakpoint is defined"
                          "for the specified wavelength"])


def get_t2(lambda_, alpha):
    """
    :param lambda_: wavelength [nm]
    :param alpha: angular subtense of the apparent source
    :return: breakpoint T2 [s]
    """
    alpha_min = 1.5e-3
    if 400 <= lambda_ <= 1400:
        if alpha <= alpha_min:
            t2 = 10
        elif alpha_min <= alpha <= 100e-3:
            t2 = 10 * 10**((alpha-alpha_min)/98.5e-3)
        else:  # 100e-3 <= alpha
            t2 = 100
    else:
        raise ValueError(["No breakpoint is defined for"
                          "the specified wavelength"])
    return t2


def get_c1(lambda_, t):
    """
    :param lambda_: wavelength [nm]
    :param t: emission duration [s]
    :return: correction factor C1
    """
    # c1
    if 180 <= lambda_ <= 400:
        return 5.6e-3 * t**0.25
    else:
        raise ValueError(["No correction factor is defined"
                          "for the specified wavelength"])


def get_c2(lambda_):
    """
    :param lambda_: wavelength [nm]
    :return: correction factor C2
    """
    if 180 <= lambda_ < 302.5:
        c2 = 30
    elif 302.5 <= lambda_ <= 315:
        c2 = 10**(0.2 * (lambda_-295))
    else:
        raise ValueError(["No correction factor is defined"
                          "for the specified wavelength"])
    return c2


def get_c3(lambda_):
    """
    :param lambda_: wavelength [nm]
    :return: correction factor C3
    """
    if 400 <= lambda_ < 450:
        c3 = 1.0
    elif 450 <= lambda_ < 600:
        c3 = 10**(0.02*(lambda_-450))
    else:
        raise ValueError(["No correction factor is defined"
                          "for the specified wavelength"])
    return c3


def get_c4(lambda_):
    """
    :param lambda_: wavelength [nm]
    :return: correction factor C4
    """
    if 700 <= lambda_ < 1050:
        c4 = 10**(0.002*(lambda_-700))
    elif 1050 <= lambda_ < 1400:
        c4 = 5
    else:
        raise ValueError(["No correction factor is defined"
                          "for the specified wavelength"])
    return c4


def get_c5(lambda_, t_pulse, t_period, t_base, alpha):
    """
    Finds the value of correction factor C5 for reduced pulse analysis
    :param lambda_: wavelength [nm]
    :param t_pulse: pulse duration
    :param t_period: length of the laser pulse period (1/f)
    :param t_base: time base
    :param alpha: angular subtense of the apparent source
    :return: C5, effective pulse period,
            number of real pulses in effective pulse
    """
    # time below which pulse groups are summed
    #   (i.e., min effective pulse length)
    t_i = get_t_i(lambda_)
    # number of real pulses in min effective pulse duration
    n = max(t_i / t_period, 1)
    # if there are multiple pulses in t_i (see 4.3.f)
    # (i.e., real pulses are too short),
    # we need to form effective pulses of length t_i
    t_pulse_eff = max(t_i, t_period)

    # breakpoint for calculating c5 (from table 9)
    t2 = get_t2(lambda_, alpha)
    if t2 is None:
        print('T2 not provided, will use T (time base)')
        t_max = t_base
    else:
        t_max = min(t2, t_base)

    # Find the number of effective pulses in the considered
    # emission duration.
    n_eff = t_max / t_pulse_eff

    if t_pulse > 0.25:
        return 1, t_pulse_eff, n    # c5 is not applicable

    if t_pulse <= t_i:
        if t_max <= 0.25:
            c5 = 1
        else:  # t_max > 0.25
            if n_eff <= 600:
                c5 = 1
            else:  # n_eff > 600
                c5 = max(0.4, 5 * n_eff**-0.25)
    else:  # t_pulse > t_i
        # the value of angular subtense of the apparent source
        #     above which the MPEs and AELs are independent of the source size
        alpha_max = get_alpha_max(t_pulse)
        if alpha <= 5e-3:
            c5 = 1
        elif 5e-3 < alpha <= alpha_max:
            if n_eff <= 40:
                c5 = n_eff**-0.25
            else:  # n_eff > 40
                c5 = 0.4
        else:  # alpha > alpha_max
            if n_eff <= 625:
                c5 = n_eff**-0.25
            else:  # n_eff > 625
                c5 = 0.2
        if alpha > 100e-3:
            c5 = 1

    return c5, t_pulse_eff, n


def get_c6(lambda_, alpha, t):
    """
    :param lambda_: wavelength [nm]
    :param alpha: angular subtense of the apparent source
    :param t: emission duration
    :return: correction factor C6
    """
    alpha_min = 1.5e-3
    alpha_max = get_alpha_max(t)
    # alpha = max(alpha, alpha_min)
    # Is this only if lambda in [302.5, 4000], due
    # to limiting the angle of acceptance
    # (for thermal hazard limits)?
    alpha = min(alpha, alpha_max)
    if 180 <= lambda_ < 400 or 1400 <= lambda_ <= 1e-6:
        c6 = 1
    elif 400 <= lambda_ < 1400:
        if alpha <= alpha_min:
            c6 = 1
        elif alpha_min <= alpha <= alpha_max:
            c6 = alpha / alpha_min
        else:  # alpha_max <= alpha
            c6 = alpha_max / alpha_min
    else:
        raise ValueError(["No correction factor is defined"
                          "for the specified wavelength"])
    return c6


def get_c7(lambda_):
    """
    :param lambda_: wavelength [nm]
    :return: correction factor C7
    """
    if 700 <= lambda_ < 1150:
        c7 = 1
    elif 1150 <= lambda_ < 1200:
        c7 = 10**(0.018*(lambda_-1150))
    elif 1200 <= lambda_ <= 1400:
        c7 = 8 + 10 ** (0.04 * (lambda_ - 1250))
    else:
        raise ValueError(["No correction factor is defined"
                          "for the specified wavelength"])
    return c7


def get_t_i(lambda_):
    """
    Times below which pulse groups are summed
    :param lambda_: wavelength [nm]
    :return: t_i [s]
    """
    t_i = None
    if 400 <= lambda_ < 1050:
        t_i = 5e-6
    if 1050 <= lambda_ < 1400:
        t_i = 13e-6
    if 1400 <= lambda_ < 1500:
        t_i = 1e-3
    if 1500 <= lambda_ < 1800:
        t_i = 10
    if 1800 <= lambda_ < 2600:
        t_i = 1e-3
    if 2600 <= lambda_ < 1e6:
        t_i = 1e-7

    return t_i
