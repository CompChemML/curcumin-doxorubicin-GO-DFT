"""
extract_energy.py

Utility functions used throughout this project to pull data out of raw
ORCA .out files.

1. final_single_point_energy(path)
   Returns the last "FINAL SINGLE POINT ENERGY" value (in Hartree) found
   in an ORCA output file. Used for every drug / GO / complex energy in
   this project.

2. extract_last_cartesian(path, natoms)
   Returns the last block of `natoms` Cartesian coordinate lines
   (element, x, y, z) from an ORCA output file's
   "CARTESIAN COORDINATES (ANGSTROEM)" sections.

   This was needed specifically for curcumin: the original isolated
   geometry optimization job did not have a companion .xyz file saved
   alongside curcumin.out, so the converged geometry had to be pulled
   directly out of the .out file's last coordinate block instead.

Usage:
    python extract_energy.py path/to/file.out
    python extract_energy.py path/to/curcumin.out --xyz 47 > curcumin.xyz
"""
import argparse
import re


def final_single_point_energy(path: str) -> float:
    energy = None
    with open(path, "r", errors="ignore") as f:
        for line in f:
            if "FINAL SINGLE POINT ENERGY" in line:
                m = re.search(r"(-?\d+\.\d+)", line)
                if m:
                    energy = float(m.group(1))
    if energy is None:
        raise ValueError(f"No FINAL SINGLE POINT ENERGY found in {path}")
    return energy


def extract_last_cartesian(path: str, natoms: int):
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()

    last_start = None
    for i, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            last_start = i

    if last_start is None:
        raise ValueError(f"No CARTESIAN COORDINATES block found in {path}")

    # ORCA format: header line, dashed underline, then natoms coordinate lines
    coord_start = last_start + 2
    coords = []
    for line in lines[coord_start:coord_start + natoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Malformed coordinate line in {path}: {line!r}")
        sym, x, y, z = parts[0], parts[1], parts[2], parts[3]
        coords.append((sym, float(x), float(y), float(z)))

    if len(coords) != natoms:
        raise ValueError(
            f"Expected {natoms} atoms, got {len(coords)} in {path}"
        )
    return coords


def write_xyz(coords, comment=""):
    lines = [str(len(coords)), comment]
    for sym, x, y, z in coords:
        lines.append(f"{sym:2s}  {x:14.8f}  {y:14.8f}  {z:14.8f}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("out_file")
    parser.add_argument("--xyz", type=int, metavar="NATOMS",
                         help="Instead of printing the energy, extract the "
                              "last geometry block with NATOMS atoms as XYZ")
    args = parser.parse_args()

    if args.xyz:
        coords = extract_last_cartesian(args.out_file, args.xyz)
        print(write_xyz(coords, comment=f"extracted from {args.out_file}"))
    else:
        print(f"{final_single_point_energy(args.out_file):.12f}")
