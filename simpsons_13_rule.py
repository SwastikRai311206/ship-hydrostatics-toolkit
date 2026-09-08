"""
simpsons_13_rule.py
----------------
Simpson's 1/3 Rule for numerical integration, applied to a table
of half-ordinates (half-breadths) at equally spaced stations.

Commonly used in naval architecture to calculate hull properties
- such as waterplane area, displacement, centers of flotation (LCF)
and buoyancy (LCB), and moments of inertia -- from a table
of offsets. Unlike closed-form geometric formulas(for eg -> Wall sided formula), Simpson's Rule
requires no assumption about hull shape, making it applicable to
any hull form.

THEORY
------
simpson's 1/3rd rule calculates the area under the cuve , while approximating
the curve to be parabola 

Conditions:
1. Theres must be an odd number of ordinates
2. The ordinates must be equally spaced, with spacing h.

Area = (h / 3) * [y0 + 4*y1 + 2*y2 + 4*y3 + ... + 4*y(n-2) + y(n-1)]

The repeating multiplier pattern  1, 4, 2, 4, 2, ..., 4, 1  is called
the "Simpson's multipliers" (SM).

Parameters Analyzed using Simpson's 1/3rd rule
----------------------------------------------
Waterplane Area, Moment of Inertia , Longitudinal Centre of Flotation,


For a set of HALF-ordinates (breadth measured from the centerline
out to one side of the hull), the results must be DOUBLED to get
the full waterplane property (since the hull is symmetric port/
starboard):
 
    Waterplane Area   A_w  = 2 * (h/3) * Sum(SM * y)
    Moment of Inertia I_CL = (2h/9) * Sum(SM * y^3)
 
    (I_CL comes from integrating the strip inertia (2y)^3/12 along
    the length -- see derivation in the docstring of
    moment_of_inertia_simpsons below.)

"""

import numpy as np

def simpsons_multipliers (n):
    """
    Generate pattern of Simpson's multiplier ( 1,4,2,4,...,4,1)
    for n number of ordinates (n must be odd).
    """

    if (n<3):
        raise ValueError ( f"Only {n} ordinates were entered. "

            "At least three ordinates are required for Simson's 1/3 rule")
    if (n%2 == 0):
        raise ValueError ("The number of ordinates must be odd to apply Simpson's 1/3 rule")
    sm = np.ones(n)
    sm[1:-1:2] = 4              # Assign multiplier 4 to odd-indexed interior ordinates
    sm[2:-1:2] = 2              # Assign multiplier 2 to even-indexed interior ordinates
    return sm

def simpsons_integral(y,h):
    '''
    Integrates a set of ordinates y, equally spaced by a distance h, 
    using Simpson's 1/3rd Rule.

    Parameters
    ----------

    y : array-like
        ordinate values
    h : float
        Equal spacing between ordinates

    Returns
    ------
    float
         Approximate area under the curve 

    '''
    y = np.asarray(y, dtype=float)
    sm = simpsons_multipliers(len(y))

    return h/3 * np.sum(sm*y)

def waterplane_area_simpsons(half_ordinates, h):
    '''
    Waterplane area from a table of HALF-ordinates (half-breadths).
 
    A_w = 2 * (h/3) * Sum(SM * y)
 
    The factor of 2 accounts for both port and starboard sides,
    since half_ordinates only measure from the centerline outward
    '''
    return 2.0 * simpsons_integral(half_ordinates,h)

def moment_of_inertia_simpsons(half_ordinates, h):
    '''
    Tranverse moment of inertia of the waterpalne about the
    centreline, from a table of Half-ordinates.

    Derivation : a thin strip of the waterplane at a station, of
    length dx along the ship and FULL breadth (2y), has its own
    moment of inertia about the centerline of (2y)^3 / 12 * dx
    (standard rectangle inertia formula, b^3*dx/12).
 
        I_CL = Integral over length of  (2y)^3 / 12  dx
             = (2/3) * Integral of y^3 dx
             = (2/3) * (h/3) * Sum(SM * y^3)      [Simpson's Rule]
             = (2h/9) * Sum(SM * y^3)
 
    This I_CL is exactly the quantity needed for BM = I_CL / V.
    '''
    half_ordinates = np.asarray(half_ordinates,dtype=float)
    sm = simpsons_multipliers(len(half_ordinates))
    return (2.0*h/9.0) * np.sum(sm*half_ordinates**3)

def lcf_simpsons(half_ordinates, h, reference="midship"):
    """
    Longitudinal Center of Flotation (LCF): the longitudinal centroid of the
    waterplane area, calculated using Simpson's Rule moments.

    It is the longitudinal point about which the ship trims for small changes
    in trim.
    
    Theory
    ------
    Area is proportional to  Sum(SM * y).
    First moment of area about a reference is proportional to
    Sum(SM * y * lever), where lever is each station's distance from that
    reference point.
 
    LCF (distance from reference) = Sum(SM * y * lever) / Sum(SM * y)
 
    Note: the "doubling for port/starboard" and "h/3" factors used
    for Area cancel out of this ratio, so LCF can be computed
    directly from the raw Simpson-weighted sums.
 
    Parameters
    ----------
    half_ordinates : array-like
        Half-breadth (m) at each station.
    h : float
        Common spacing between stations (m).
    reference : "midship" or "AP"
        "AP"      -> LCF reported as distance forward of the aft
                     perpendicular (station 0).
        "midship" -> LCF reported as +/- distance from midships
                     (positive = forward of midship, negative = aft),
                     the convention normally used on stability reports.
 
    Returns
    -------
    float : LCF position (m), per the chosen reference/sign convention.
    """
    y = np.asarray(half_ordinates,dtype=float)
    n = len(y)
    sm = simpsons_multipliers(n)

    levers = np.arange(n) * h     #  Lever arms measured from AP (Aft Perpendicular)

    moment_function = np.sum(y * sm * levers)
    Area_function = np.sum(sm*y)

    lcf_from_AP = moment_function/Area_function

    if (reference == "AP"):
        return lcf_from_AP
    elif (reference == "midship"):
        midship_from_AP = (n - 1) / 2.0 * h
        return lcf_from_AP - midship_from_AP   # + fwd, - aft of midship
    else:
        raise ValueError("reference must be 'AP' or 'midship'")
 

if (__name__ == "__main__"):
    # -----------------------------------------------------------
    # Worked example: half-ordinates (half-breadths, in metres)
    # measured at 7 equally spaced stations along a barge's
    # waterline (0 = one end, 6 = the other end). 7 ordinates
    # -> 6 intervals -> even, so Simpson's 1/3rd Rule applies.
    # -----------------------------------------------------------
    # Asymmetric hull (fuller aft, finer forward) so LCF is not
    # trivially at midship -- more representative of a real vessel.

    half_ordinates = [0.0, 4.5, 5.8, 6.0, 5.9, 4.8, 0.0]   #meters
    L = 40.0                    # total waterline length, m
    n_intervals = len(half_ordinates) - 1
    h = L/n_intervals        # common spacing between ordinates


    print("=" * 50)
    print("SIMPSON'S 1/3rd RULE - WATERPLANE FROM HALF-ORDINATES")
    print("=" * 50)
    print(f"Stations (half-ordinates, m): {half_ordinates}")
    print(f"Number of ordinates          : {len(half_ordinates)}")
    print(f"Common spacing h              : {h:.3f} m")

    sm = simpsons_multipliers(len(half_ordinates))
    print(f"Simpson's multipliers : {sm.astype(int).tolist()}")

    A_w = waterplane_area_simpsons(half_ordinates,h)
    I_CL = moment_of_inertia_simpsons(half_ordinates,h)
    LCF_ap = lcf_simpsons(half_ordinates, h, reference="AP")
    LCF_midship = lcf_simpsons(half_ordinates, h, reference="midship")
    Load_TPC = A_w/97.56
    

    print("-"*50)
    print(f"waterplane area (A_w)         : {A_w:.3f} m^2")
    print(f"Moment of inertia (I_CL)      : {I_CL:.3f} m^4")
    print(f"LCF from AP (station 0)       : {LCF_ap:.3f} m")
    print(f"LCF from midship              : {LCF_midship:+.3f} m"
           f"({'fwd of' if LCF_midship >= 0 else 'aft of'} midship)")
    print(f"Tonnes per centimeter Immersion (TPC) : {Load_TPC:.3f} t/cm")

    # ------------------------------------------------------------
    # Compare against the equivalent box-hull formula, using the
    # max breadth, as a sanity check (won't match exactly, since
    # this hull tapers at the ends -- that's the whole point of
    # using Simpson's Rule instead of a box formula).
    # ------------------------------------------------------------
    B_max = 2 * max(half_ordinates)
    A_box = L * B_max
    I_box = L * B_max**3 / 12
    print("-" * 50)
    print("For comparison, an equivalent full RECTANGULAR box")
    print(f"(L={L} m x B={B_max} m) would give:")
    print(f"  Box area                    : {A_box:.3f} m^2")
    print(f"  Box moment of inertia        : {I_box:.3f} m^4")
    print("(Simpson's result is smaller because this hull tapers")
    print(" at bow/stern, unlike a full rectangular barge.)")

 

    






    
