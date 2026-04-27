"""
get_parameters.py
=================
Derives ITU-R channel parameters for a satellite link given
latitude, longitude, and frequency.

What this script computes:
    - k, alpha coefficients          (ITU-R P.838)
    - Rainfall rate R_001            (ITU-R P.837)
    - 0 deg C isotherm height h_0    (ITU-R P.839)
    - [Optional] Gaseous atten.      (ITU-R P.676)
    - [Optional] Cloud atten.        (ITU-R P.840)
    - [Optional] Scintillation       (ITU-R P.618)

What is NOT computed here (done in channel_tropical.m):
    - Specific rain attenuation gamma_R = k * R^alpha
    - Rain attenuation A_p via CDF scaling
    - Channel coefficient h_k

Assumed fixed:
    - Elevation angle = 45 degrees
    - Polarization    = V (vertical)

Usage:
    python get_parameters.py --lat 5.17 --lon 100.4 --freq 20
    python get_parameters.py --lat 5.17 --lon 100.4 --freq 20 --optional
    python get_parameters.py --lat 5.17 --lon 100.4 --freq 20 --output penang.json

Install:
    pip install itur numpy
"""

import argparse
import json
import sys
import numpy as np

try:
    import itur
    import itur.models.itu837 as itu837
    import itur.models.itu838 as itu838
    import itur.models.itu839 as itu839
    import itur.models.itu840 as itu840
    import itur.models.itu618 as itu618
    import itur.models.itu676 as itu676
except ImportError:
    print("ERROR: itur library not found.")
    print("Install with:  pip install itur")
    sys.exit(1)

# Fixed assumptions
ELEVATION_DEG = 45.0
POLARIZATION  = 'V'    # vertical -> tau = 0 
TAU           = 0      # tilt angle for P.838 (0=vertical, 90=horizontal)


def get_parameters(lat, lon, freq_GHz, include_optional=False,
                   altitude_m=0.0, D_antenna_m=1.0, eta=0.6):
    """
    Compute ITU-R parameters for a given location and frequency.

    Parameters
    ----------
    lat             : float  — latitude (deg N, negative = South)
    lon             : float  — longitude (deg E, negative = West)
    freq_GHz        : float  — frequency in GHz
    include_optional: bool   — if True, also compute A_gas, A_cloud, A_scint
    altitude_m      : float  — station altitude above sea level (m)
    D_antenna_m     : float  — antenna diameter for scintillation (m)
    eta             : float  — antenna efficiency for scintillation

    Returns
    -------
    dict — all computed parameters
    """

    alt_km = altitude_m / 1000.0
    f      = freq_GHz

    print(f"\nComputing ITU-R parameters for:")
    print(f"  Location  : ({lat}N, {lon}E)  Alt: {altitude_m} m")
    print(f"  Frequency : {f} GHz")
    print(f"  Elevation : {ELEVATION_DEG} deg (fixed)")
    print(f"  Pol.      : {POLARIZATION} (fixed)\n")

    # ── 1. k and alpha  (ITU-R P.838) ────────────────────────────────────
    # Returns (k, alpha) tuple for given frequency, elevation, tilt angle
    k_v, alpha_v = itu838.rain_specific_attenuation_coefficients(
        f, el=ELEVATION_DEG, tau=0)      # vertical pol (tau=0)

    k_h, alpha_h = itu838.rain_specific_attenuation_coefficients(
        f, el=ELEVATION_DEG, tau=90)     # horizontal pol (tau=90)

    print(f"[P.838] k_V     = {k_v:.8f}")
    print(f"[P.838] alpha_V = {alpha_v:.8f}")
    print(f"[P.838] k_H     = {k_h:.8f}")
    print(f"[P.838] alpha_H = {alpha_h:.8f}")

    # ── 2. Rainfall rate R_001  (ITU-R P.837) ────────────────────────────
    R_001 = itu837.rainfall_rate(lat, lon, p=0.01).value
    print(f"\n[P.837] R_001   = {R_001:.4f} mm/h")

    # ── 3. 0 deg C isotherm height h_0  (ITU-R P.839) ────────────────────
    # P.839 returns the mean annual 0 deg C isotherm height above sea level
    h_0 = itu839.isoterm_0(lat, lon).value
    print(f"[P.839] h_0     = {h_0:.4f} km  (0 deg C isotherm height)")

    # Also get rain height h_R = h_0 + 0.36  (ITU-R P.839 Eq. 2)
    h_R = itu839.rain_height(lat, lon).value
    print(f"[P.839] h_R     = {h_R:.4f} km  (rain height = h_0 + 0.36)")

    # ── 4. Optional: gaseous, cloud, scintillation ────────────────────────
    A_gas   = None
    A_cloud = None
    A_scint = None

    if include_optional:
        print("\n[Optional parameters]")

        # Gaseous attenuation (ITU-R P.676)
        A_gas = itu676.gaseous_attenuation_slant_path(
            f   = f,
            el  = ELEVATION_DEG,
            rho = 7.5,       # surface water vapour density g/m3 (standard)
            P   = 1013.25,   # surface pressure hPa
            T   = 15,        # surface temperature deg C
            h   = alt_km,
        ).value
        print(f"[P.676] A_gas   = {A_gas:.4f} dB")

        # Cloud attenuation (ITU-R P.840) at p=0.01%
        A_cloud = itu840.cloud_attenuation(
            lat, lon,
            el = ELEVATION_DEG,
            f  = f,
            p  = 0.01,
        ).value
        print(f"[P.840] A_cloud = {A_cloud:.4f} dB")

        # Scintillation attenuation (ITU-R P.618) at p=0.01%
        A_scint = itu618.scintillation_attenuation(
            lat, lon,
            f   = f,
            el  = ELEVATION_DEG,
            p   = 0.01,
            D   = D_antenna_m,
            eta = eta,
        ).value
        print(f"[P.618] A_scint = {A_scint:.4f} dB")

    # ── Pack results ──────────────────────────────────────────────────────
    def get_val(item):
        if hasattr(item, 'value'):
            return item.value
        return item
    params = {
        "metadata": {
            "latitude_deg":   lat,
            "longitude_deg":  lon,
            "altitude_m":     altitude_m,
            "frequency_GHz":  f,
            "elevation_deg":  ELEVATION_DEG,
            "polarization":   POLARIZATION,
        },
        "itu_p838": {
            "k_V":      round(float(k_v),     8),
            "alpha_V":  round(float(alpha_v), 8),
            "k_H":      round(float(k_h),     8),
            "alpha_H":  round(float(alpha_h), 8),
        },
        "itu_p837": {
            "R_001_mmh": round(float(R_001), 6),
        },
        "itu_p839": {
            "h_0_km": round(float(get_val(h_0)), 6),
            "h_R_km": round(float(get_val(h_R)), 6),
        },
        "itu_optional": {
            "A_gas_dB":   round(float(A_gas),   4) if A_gas   is not None else None,
            "A_cloud_dB": round(float(A_cloud), 4) if A_cloud is not None else None,
            "A_scint_dB": round(float(A_scint), 4) if A_scint is not None else None,
            "computed":   include_optional,
        },
        "matlab_ready": {
            # Flat struct — load directly in MATLAB via load_itu_params.m
            "k_V":        round(float(k_v),     8),
            "alpha_V":    round(float(alpha_v), 8),
            "k_H":        round(float(k_h),     8),
            "alpha_H":    round(float(alpha_h), 8),
            "R_001_mmh":  round(float(R_001),   6),
            "h_0_km":     round(float(h_0),     6),
            "h_R_km":     round(float(h_R),     6),
            "A_gas_dB":   round(float(A_gas),   4) if A_gas   is not None else 0.0,
            "A_cloud_dB": round(float(A_cloud), 4) if A_cloud is not None else 0.0,
            "A_scint_dB": round(float(A_scint), 4) if A_scint is not None else 0.0,
        }
    }
    return params


def print_summary(params):
    mr   = params["matlab_ready"]
    meta = params["metadata"]
    opt  = params["itu_optional"]

    print("\n" + "="*55)
    print("  ITU-R Parameters Summary")
    print(f"  ({meta['latitude_deg']}N, {meta['longitude_deg']}E) | "
          f"{meta['frequency_GHz']} GHz | El: {meta['elevation_deg']} deg")
    print("="*55)
    rows = [
        ("k_V   (P.838)",   mr["k_V"]),
        ("alpha_V (P.838)", mr["alpha_V"]),
        ("k_H   (P.838)",   mr["k_H"]),
        ("alpha_H (P.838)", mr["alpha_H"]),
        ("R_001  (P.837, mm/h)", mr["R_001_mmh"]),
        ("h_0   (P.839, km)",    mr["h_0_km"]),
        ("h_R   (P.839, km)",    mr["h_R_km"]),
    ]
    for name, val in rows:
        print(f"  {name:<28} {val:>14.8f}")

    if opt["computed"]:
        print(f"\n  {'A_gas   (P.676, dB)':<28} {mr['A_gas_dB']:>14.4f}")
        print(f"  {'A_cloud (P.840, dB)':<28} {mr['A_cloud_dB']:>14.4f}")
        print(f"  {'A_scint (P.618, dB)':<28} {mr['A_scint_dB']:>14.4f}")
    else:
        print("\n  Optional params (A_gas, A_cloud, A_scint): not computed")
        print("  Run with --optional to include them.")

    print("="*55)
    


def main():
    parser = argparse.ArgumentParser(
        description="Derive ITU-R parameters (k, alpha, R_001, h_0) from lat/lon/freq.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--lat",      type=float, required=True,
                        help="Latitude in degrees N (negative = South)")
    parser.add_argument("--lon",      type=float, required=True,
                        help="Longitude in degrees E (negative = West)")
    parser.add_argument("--freq",     type=float, required=True,
                        help="Frequency in GHz")
    parser.add_argument("--alt",      type=float, default=0.0,
                        help="Station altitude above sea level (m)")
    parser.add_argument("--optional", action="store_true",
                        help="Also compute A_gas, A_cloud, A_scint")
    parser.add_argument("--diameter", type=float, default=1.0,
                        help="Antenna diameter in m (for scintillation)")
    parser.add_argument("--eta",      type=float, default=0.6,
                        help="Antenna efficiency (for scintillation)")
    parser.add_argument("--output",   type=str,
                        default="channel_params.json",
                        help="Output JSON filename")

    args = parser.parse_args()

    params = get_parameters(
        lat              = args.lat,
        lon              = args.lon,
        freq_GHz         = args.freq,
        include_optional = args.optional,
        altitude_m       = args.alt,
        D_antenna_m      = args.diameter,
        eta              = args.eta,
    )

    print_summary(params)

    with open(args.output, "w") as f:
        json.dump(params, f, indent=2)

    print(f"Saved to: {args.output}")
    print(f"\nLoad in MATLAB:")
    print(f"  itu = load_itu_params('{args.output}');")


if __name__ == "__main__":
    main()
