"""
build_d3bj_and_solvation_inputs.py

Generates the single-point input files for the D3BJ-dispersion-only and
D3BJ+SMD(water)-solvation correction levels used in this project. These
are simpler than the BSSE files: no ghost atoms, just the existing
converged geometry (isolated monomer or full complex) re-run as a single
point with extra keywords on the "!" line.

D3BJ only:
    ! B3LYP D3BJ 6-31G* SP

D3BJ + implicit water (SMD model):
    ! B3LYP D3BJ SMD(water) 6-31G* SP

Note: solvation here is a single point on the gas-phase-optimized
geometry (not re-optimized in solvent), combined with D3BJ. It was
deliberately NOT combined with the BSSE ghost-atom counterpoise method
in this project -- CP correction and continuum solvation correction are
reported as two separate, non-additive refinements on top of the
raw B3LYP binding energy, not stacked together.

Usage:
    python build_d3bj_and_solvation_inputs.py curcumin_on_go_restart.xyz curcumin_on_go --out inputs/d3bj
    python build_d3bj_and_solvation_inputs.py curcumin_on_go_restart.xyz curcumin_on_go --out inputs/solvation --solvation
"""
import argparse
import os

D3BJ_HEADER = "! B3LYP D3BJ 6-31G* SP\n%maxcore 800\n\n* xyz 0 1\n"
SOLV_HEADER = "! B3LYP D3BJ SMD(water) 6-31G* SP\n%maxcore 800\n\n* xyz 0 1\n"


def read_xyz_block(path):
    with open(path) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    return "".join(lines[2:2 + natoms])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz_file", help="Converged geometry (monomer or complex)")
    ap.add_argument("label", help="Base name for the output .inp, e.g. curcumin_on_go")
    ap.add_argument("--out", required=True)
    ap.add_argument("--solvation", action="store_true",
                     help="Add SMD(water) on top of D3BJ (default: D3BJ only)")
    args = ap.parse_args()

    header = SOLV_HEADER if args.solvation else D3BJ_HEADER
    suffix = "_solv" if args.solvation else "_D3"

    os.makedirs(args.out, exist_ok=True)
    coords = read_xyz_block(args.xyz_file)

    out_path = os.path.join(args.out, f"{args.label}{suffix}.inp")
    with open(out_path, "w") as f:
        f.write(header)
        f.write(coords)
        f.write("*\n")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
