"""
generate_figures.py

Generates every image in figures/:
  - 2D skeletal structures of curcumin and doxorubicin (RDKit, from SMILES)
  - 3D renders of the converged DFT geometries: bare GO cluster,
    curcumin-on-GO complex, doxorubicin-on-GO complex (ASE, matplotlib backend)

Run from the repo root:
    python scripts/generate_figures.py
"""
import os
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import numpy as np
from ase.io import read
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.neighborlist import neighbor_list
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
OUT = os.path.join(ROOT, "outputs", "optimization")
os.makedirs(FIG, exist_ok=True)

SMILES = {
    "curcumin": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
    "doxorubicin": "COc1cccc2c1C(=O)c1c(O)c3c(c(O)c1C2=O)C[C@](O)(C(C)=O)C[C@@H]3O[C@H]1C[C@H](N)[C@H](O)[C@H](C)O1",
}


def draw_2d():
    for name, smi in SMILES.items():
        mol = Chem.MolFromSmiles(smi)
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(900, 700))
        path = os.path.join(FIG, f"{name}_2D_structure.png")
        img.save(path)
        print(f"wrote {path}")


def draw_3d(xyz_path, title, out_name, rotation="20x,10y,0z"):
    """Simple self-contained ball-and-stick renderer (painter's algorithm):
    rotate atoms, project to x/y, draw bonds then atoms back-to-front by z."""
    atoms = read(xyz_path)
    atoms.center()

    for token in rotation.split(","):
        token = token.strip()
        if not token:
            continue
        angle = float(token[:-1])
        axis = token[-1]
        if angle != 0:
            atoms.rotate(angle, axis)

    pos = atoms.get_positions()
    numbers = atoms.numbers
    radii = np.array([covalent_radii[n] for n in numbers])
    colors = [jmol_colors[n] for n in numbers]

    cutoffs = radii * 1.25
    i_list, j_list = neighbor_list("ij", atoms, cutoffs)
    bonds = [(i, j) for i, j in zip(i_list, j_list) if i < j]

    # z-order draw list: bonds tagged with their midpoint depth, atoms with
    # their own depth; sort everything back-to-front (lowest z first).
    draw_items = []
    for i, j in bonds:
        z_mid = (pos[i, 2] + pos[j, 2]) / 2
        draw_items.append((z_mid, "bond", (i, j)))
    for idx in range(len(atoms)):
        draw_items.append((pos[idx, 2], "atom", idx))
    draw_items.sort(key=lambda t: t[0])

    span = max(np.ptp(pos[:, 0]), np.ptp(pos[:, 1]), 1.0)
    fig, ax = plt.subplots(figsize=(8, 8))

    for _, kind, payload in draw_items:
        if kind == "bond":
            i, j = payload
            ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]],
                    color="#888888", linewidth=2.2, solid_capstyle="round", zorder=1)
        else:
            idx = payload
            r = max(0.28 * radii[idx], 0.14)
            circ = Circle((pos[idx, 0], pos[idx, 1]), radius=r,
                          facecolor=colors[idx], edgecolor="black",
                          linewidth=0.6, zorder=2)
            ax.add_patch(circ)

    pad = 0.08 * span
    ax.set_xlim(pos[:, 0].min() - pad, pos[:, 0].max() + pad)
    ax.set_ylim(pos[:, 1].min() - pad, pos[:, 1].max() + pad)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(title, fontsize=14)
    path = os.path.join(FIG, out_name)
    fig.savefig(path, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"wrote {path}")


def draw_all_3d():
    draw_3d(os.path.join(OUT, "go_cluster.xyz"),
            "Graphene Oxide cluster (28 atoms, DFT-optimized)",
            "go_cluster_3D.png", rotation="20x,10y,0z")
    draw_3d(os.path.join(OUT, "curcumin_on_go_FINAL.xyz"),
            "Curcumin adsorbed on GO (converged complex)",
            "curcumin_on_GO_3D.png", rotation="20x,10y,0z")
    draw_3d(os.path.join(OUT, "doxorubicin_on_go_FINAL.xyz"),
            "Doxorubicin adsorbed on GO (converged complex)",
            "doxorubicin_on_GO_3D.png", rotation="20x,10y,0z")


if __name__ == "__main__":
    draw_2d()
    draw_all_3d()
