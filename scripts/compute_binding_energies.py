"""
compute_binding_energies.py

Recomputes every binding-energy row in data/results.csv directly from the
raw ORCA .out files in outputs/, so the numbers in the README/dashboard
are always reproducible from source rather than hand-copied.

E_binding = E_complex - (E_drug + E_GO)     [simple / D3BJ-only / solvation levels]
E_binding_CP = E_complex - (E_dimerbasis_drug + E_dimerbasis_GO)   [BSSE-corrected levels]

1 Hartree = 627.5094740631 kcal/mol

Usage (run from the repo root):
    python scripts/compute_binding_energies.py
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_energy import final_single_point_energy as E

HARTREE_TO_KCAL = 627.5094740631

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "outputs")


def p(*parts):
    return os.path.join(OUT, *parts)


def kcal(hartree):
    return hartree * HARTREE_TO_KCAL


def report(label, e_complex, e_a, e_b):
    binding = e_complex - (e_a + e_b)
    print(f"{label:32s} {kcal(binding):8.2f} kcal/mol")
    return kcal(binding)


def main():
    print("=== Curcumin - GO ===")
    report("raw B3LYP",
           E(p("optimization", "curcumin_on_go_FINAL.out")),
           E(p("optimization", "curcumin.out")),
           E(p("optimization", "go_cluster.out")))
    report("BSSE-corrected (CP)",
           E(p("optimization", "curcumin_on_go_FINAL.out")),
           E(p("bsse", "curcumin_dimerbasis_curcumin.out")),
           E(p("bsse", "curcumin_dimerbasis_GO.out")))
    report("D3BJ only",
           E(p("d3bj", "curcumin_on_go_D3.out")),
           E(p("d3bj", "curcumin_D3.out")),
           E(p("d3bj", "go_cluster_D3.out")))
    report("BSSE + D3BJ",
           E(p("d3bj", "curcumin_on_go_D3.out")),
           E(p("bsse", "curcumin_dimerbasis_curcumin_D3.out")),
           E(p("bsse", "curcumin_dimerbasis_GO_D3.out")))
    report("D3BJ + SMD(water)",
           E(p("solvation", "curcumin_on_go_solv.out")),
           E(p("solvation", "curcumin_solv.out")),
           E(p("solvation", "go_cluster_solv.out")))

    print("\n=== Doxorubicin - GO ===")
    report("raw B3LYP",
           E(p("optimization", "doxorubicin_on_go_FINAL.out")),
           E(p("optimization", "doxorubicin.out")),
           E(p("optimization", "go_cluster.out")))
    report("BSSE-corrected (CP)",
           E(p("optimization", "doxorubicin_on_go_FINAL.out")),
           E(p("bsse", "doxorubicin_dimerbasis_doxorubicin.out")),
           E(p("bsse", "doxorubicin_dimerbasis_GO.out")))
    report("D3BJ only",
           E(p("d3bj", "doxorubicin_on_go_D3.out")),
           E(p("d3bj", "doxorubicin_D3.out")),
           E(p("d3bj", "go_cluster_D3.out")))
    report("BSSE + D3BJ",
           E(p("d3bj", "doxorubicin_on_go_D3.out")),
           E(p("bsse", "doxorubicin_dimerbasis_doxorubicin_D3.out")),
           E(p("bsse", "doxorubicin_dimerbasis_GO_D3.out")))
    report("D3BJ + SMD(water)",
           E(p("solvation", "doxorubicin_on_go_solv.out")),
           E(p("solvation", "doxorubicin_solv.out")),
           E(p("solvation", "go_cluster_solv.out")))


if __name__ == "__main__":
    main()
