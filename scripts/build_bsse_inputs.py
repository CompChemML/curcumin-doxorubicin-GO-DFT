"""
build_bsse_inputs.py

Generates the 8 counterpoise (Boys-Bernardi) BSSE input files used in this
project (4 per drug: monobasis-drug, monobasis-GO, dimerbasis-drug,
dimerbasis-GO), plus the corresponding BSSE+D3BJ variants.

Method
------
Ghost atoms (basis functions only, no nucleus/electrons) are marked in
ORCA by appending ":" directly after the element symbol in the coordinate
block, e.g.:

    O :   1.234000   -0.567000   2.345000

For each drug-on-GO complex, at the FROZEN complex geometry:
  - monobasis-drug   = real drug atoms only, own basis only
  - monobasis-GO     = real GO atoms only, own basis only
  - dimerbasis-drug  = real drug atoms + ghost GO atoms (drug computed in
                       the full dimer basis set)
  - dimerbasis-GO    = real GO atoms + ghost drug atoms

BSSE = (E_monobasis_drug - E_dimerbasis_drug) + (E_monobasis_GO - E_dimerbasis_GO)
CP-corrected E_bind = E_complex - (E_dimerbasis_drug + E_dimerbasis_GO)

The GO cluster is always the first GO_ATOM_COUNT=28 atoms in every combined
xyz file in this project; the drug is always listed second.

Usage:
    python build_bsse_inputs.py curcumin_on_go_restart.xyz curcumin --out inputs/bsse
    python build_bsse_inputs.py doxorubicin_on_go_restart2.xyz doxorubicin --out inputs/bsse --d3bj
"""
import argparse
import os

GO_ATOM_COUNT = 28

RAW_HEADER = "! B3LYP 6-31G* SP\n%maxcore 800\n\n* xyz 0 1\n"
D3BJ_HEADER = "! B3LYP D3BJ 6-31G* SP\n%maxcore 800\n\n* xyz 0 1\n"


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    atoms = []
    for line in lines[2:2 + natoms]:
        parts = line.split()
        sym, x, y, z = parts[0], parts[1], parts[2], parts[3]
        atoms.append((sym, x, y, z))
    return atoms


def fmt_atom(sym, x, y, z, ghost=False):
    marker = " :" if ghost else "  "
    return f"  {sym}{marker}  {x}   {y}   {z}"


def build_block(go_atoms, drug_atoms, go_real, drug_real):
    lines = []
    for sym, x, y, z in go_atoms:
        lines.append(fmt_atom(sym, x, y, z, ghost=not go_real))
    for sym, x, y, z in drug_atoms:
        lines.append(fmt_atom(sym, x, y, z, ghost=not drug_real))
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("complex_xyz", help="Combined GO+drug xyz at the converged complex geometry")
    ap.add_argument("drug_name")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--d3bj", action="store_true", help="Also add D3BJ dispersion")
    args = ap.parse_args()

    atoms = read_xyz(args.complex_xyz)
    go_atoms = atoms[:GO_ATOM_COUNT]
    drug_atoms = atoms[GO_ATOM_COUNT:]

    header = D3BJ_HEADER if args.d3bj else RAW_HEADER
    suffix = "_D3" if args.d3bj else ""

    os.makedirs(args.out, exist_ok=True)

    jobs = {
        f"{args.drug_name}_monobasis_{args.drug_name}{suffix}.inp":
            build_block([], drug_atoms, go_real=False, drug_real=True),
        f"{args.drug_name}_monobasis_GO{suffix}.inp":
            build_block(go_atoms, [], go_real=True, drug_real=False),
        f"{args.drug_name}_dimerbasis_{args.drug_name}{suffix}.inp":
            build_block(go_atoms, drug_atoms, go_real=False, drug_real=True),
        f"{args.drug_name}_dimerbasis_GO{suffix}.inp":
            build_block(go_atoms, drug_atoms, go_real=True, drug_real=False),
    }

    for fname, block in jobs.items():
        path = os.path.join(args.out, fname)
        with open(path, "w") as f:
            f.write(header)
            f.write(block)
            f.write("\n*\n")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
