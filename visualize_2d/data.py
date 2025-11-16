"""
Data preparation utilities for visualize_2d.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List
import re
import xml.etree.ElementTree as ET

from rdkit import Chem
from rdkit.Chem import AllChem

INTERACTION_KIND_MAP = {
    "hb": ("Hydrogen bond", "HBond"),
    "hp": ("Hydrophobic contact", "Hydrophobic"),
    "pc": ("Pi–cation", "PiCation"),
    "ps": ("Pi–stacking", "PiStack"),
    "wb": ("Water bridge", "WaterBridge"),
    "ha": ("Halogen bond", "Halogen"),
    "sb": ("Salt bridge", "SaltBridge"),
    "me": ("Metal coordination", "Metal"),
}


@dataclass
class Interaction:
    kind: str
    ligand_atoms: list[str]
    protein_atom: str
    protein_label: str
    energy: float | None = None

    @property
    def label(self) -> str:
        human, short = INTERACTION_KIND_MAP.get(self.kind, (self.kind, self.kind))
        suffix = f" · {self.protein_label}" if self.protein_label else ""
        if self.energy is not None:
            return f"{short}{suffix} ({self.energy:.1f})"
        return f"{short}{suffix}"


@dataclass
class ComplexAssets:
    complex_dir: Path
    receptor_id: str
    ligand_id: str
    ligand: Chem.Mol
    interactions: List[Interaction]
    snapshot: Path | None
    ligand_atom_map: dict[str, int]
    smiles_map: dict[str, str]


def _normalize_color_hex(hex_str: str) -> str:
    hex_str = hex_str.strip().lower()
    if hex_str.startswith("0x"):
        return f"#{hex_str[2:]}"
    if hex_str.startswith("#"):
        return hex_str
    if re.fullmatch(r"[0-9a-f]{6}", hex_str):
        return f"#{hex_str}"
    raise ValueError(f"Unsupported color value {hex_str}")


def _find_first_png(directory: Path) -> Path | None:
    for candidate in sorted(directory.glob("*.png")):
        return candidate
    return None


def _iter_ligand_candidates(complex_dir: Path):
    patterns = [
        "ligand.*",
        "*_ligand.*",
        "ligand/ligand.*",
        "ligand/*",
    ]
    seen: set[Path] = set()
    for pattern in patterns:
        for candidate in complex_dir.glob(pattern):
            if candidate in seen or not candidate.is_file():
                continue
            seen.add(candidate)
            yield candidate


def _load_ligand(complex_dir: Path) -> Chem.Mol:
    candidates = list(_iter_ligand_candidates(complex_dir))
    preferred = [c for c in candidates if c.suffix.lower() in {".sdf", ".mol"}]
    if preferred:
        mol = Chem.MolFromMolFile(str(preferred[0]), removeHs=False)
        if mol:
            AllChem.Compute2DCoords(mol)
            return mol
    pdb_candidates = [c for c in candidates if c.suffix.lower() in {".pdb", ".pdbqt"}]
    for pdb in pdb_candidates:
        mol = Chem.MolFromPDBFile(str(pdb), removeHs=False)
        if mol:
            AllChem.Compute2DCoords(mol)
            return mol
    raise FileNotFoundError(f"Cannot locate ligand structure in {complex_dir}")


def _residue_label(node: ET.Element) -> str:
    restype = node.findtext("restype") or node.findtext("restype_lig") or ""
    resnr = node.findtext("resnr") or node.findtext("resnr_prot") or node.findtext("resnr_lig") or ""
    chain = node.findtext("reschain") or node.findtext("chain") or node.findtext("reschain_lig") or ""
    pieces = [restype.strip(), resnr.strip(), chain.strip()]
    return "".join(filter(None, pieces))


def _parse_plip(plip_xml: Path) -> tuple[list[Interaction], dict[str, str]]:
    tree = ET.parse(plip_xml)
    root = tree.getroot()
    out: list[Interaction] = []

    def add(kind: str, lig: Iterable[str], prot: Iterable[str], label_node: ET.Element | None = None):
        label = _residue_label(label_node) if label_node is not None else ""
        for l in lig:
            for p in prot:
                out.append(Interaction(kind=kind, ligand_atoms=[l], protein_atom=p, protein_label=label))

    for node in root.findall(".//hydrogen_bonds/hydrogen_bond"):
        prot_is_donor = (node.findtext("protisdon") or "").lower() == "true"
        lig_idx = node.findtext("acceptoridx") if prot_is_donor else node.findtext("donoridx")
        prot_idx = node.findtext("donoridx") if prot_is_donor else node.findtext("acceptoridx")
        add("hb", [lig_idx], [prot_idx], node)
    for node in root.findall(".//hydrophobic_interactions/hydrophobic_interaction"):
        add("hp", [node.findtext("ligcarbonidx")], [node.findtext("protcarbonidx")], node)
    for node in root.findall(".//pi_cation_interactions/pi_cation_interaction"):
        ligand_atoms = [idx.text for idx in node.findall("./lig_idx_list/idx") if idx.text]
        prot = [idx.text for idx in node.findall("./prot_idx_list/idx") if idx.text]
        label = _residue_label(node)
        for p in prot:
            out.append(Interaction(kind="pc", ligand_atoms=ligand_atoms, protein_atom=p, protein_label=label))
    for node in root.findall(".//pi_stacks/pi_stack"):
        add("ps", [node.findtext("ligcentroididx")], [node.findtext("protcentroididx")], node)
    for node in root.findall(".//water_bridges/water_bridge"):
        add("wb", [node.findtext("a_idx")], [node.findtext("d_idx")], node)
    for node in root.findall(".//halogen_bonds/halogen_bond"):
        add("ha", [node.findtext("halogenidx")], [node.findtext("acceptoridx")], node)
    for node in root.findall(".//salt_bridges/salt_bridge"):
        lig = [idx.text for idx in node.findall("./lig_idx_list/idx")]
        prot = [idx.text for idx in node.findall("./prot_idx_list/idx")]
        add("sb", lig, prot, node)
    for node in root.findall(".//metal_complexes/metal_complex"):
        add("me", [node.findtext("ligidx")], [node.findtext("metalidx")], node)

    mapping = {}
    mapping_node = root.find(".//mappings/smiles_to_pdb")
    if mapping_node is not None and mapping_node.text:
        pairs = mapping_node.text.strip().split(",")
        for pair in pairs:
            if ":" not in pair:
                continue
            smiles_idx, pdb_idx = pair.split(":", 1)
            mapping[smiles_idx.strip()] = pdb_idx.strip()

    return [i for i in out if i.ligand_atoms and i.protein_atom], mapping


def build_complex_assets(complex_dir: Path) -> ComplexAssets:
    complex_dir = complex_dir.resolve()
    pdb_files = list(complex_dir.glob("complex__*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No complex__*.pdb in {complex_dir}")
    base = pdb_files[0].stem.replace("complex__", "")
    try:
        receptor_id, ligand_id, _ = base.split("__")
    except ValueError:
        receptor_id = base
        ligand_id = "ligand"
    plip_xml = complex_dir / "plip" / "report.xml"
    if plip_xml.exists():
        interactions, smiles_map = _parse_plip(plip_xml)
    else:
        interactions, smiles_map = [], {}

    ligand = _load_ligand(complex_dir)
    atom_map: dict[str, int] = {}
    for atom in ligand.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info and info.GetSerialNumber():
            atom_map[str(info.GetSerialNumber())] = atom.GetIdx()
        atom_map.setdefault(str(atom.GetIdx()), atom.GetIdx())
    snapshot = _find_first_png(complex_dir)

    return ComplexAssets(
        complex_dir=complex_dir,
        receptor_id=receptor_id,
        ligand_id=ligand_id,
        ligand=ligand,
        interactions=interactions,
        snapshot=snapshot,
        ligand_atom_map=atom_map,
        smiles_map=smiles_map,
    )

