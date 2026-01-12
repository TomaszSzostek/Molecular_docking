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


def _load_ligand(complex_dir: Path, plip_xml: Path | None = None) -> Chem.Mol:
    """
    Load ligand molecule, preferring SMILES from PLIP report for accurate 2D structure.
    
    If PLIP XML is provided and contains SMILES, use it to create the molecule.
    Otherwise, fall back to loading from PDB/SDF files.
    """
    # First, try to load from SMILES in PLIP report (most accurate for 2D rendering)
    if plip_xml and plip_xml.exists():
        try:
            tree = ET.parse(plip_xml)
            root = tree.getroot()
            smiles_node = root.find(".//smiles")
            if smiles_node is not None and smiles_node.text:
                smiles = smiles_node.text.strip()
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Remove hydrogens for cleaner 2D visualization (consistent with PDB loading)
                    mol = Chem.RemoveHs(mol)
                    AllChem.Compute2DCoords(mol)
                    return mol
        except Exception:
            # If SMILES loading fails, fall through to file-based loading
            pass
    
    # Fallback: load from structure files
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


def _load_smiles_from_csv(ligand_id: str, config_path: Path | None = None) -> str | None:
    """
    Load original SMILES from ligands.csv file.
    
    This is more reliable than PLIP-extracted SMILES, which may have incorrect
    bond orders (e.g., S1(O)[O] instead of S1(=O)=O for thiazolidine derivatives).
    """
    import yaml
    import pandas as pd
    
    # Try to find config.yaml to get ligands folder path
    if config_path is None:
        # Try current directory and parent directories
        for parent in [Path.cwd(), Path(__file__).parent.parent / "run_pipeline"]:
            config_candidate = parent / "config.yaml"
            if config_candidate.exists():
                config_path = config_candidate
                break
    
    if config_path is None or not config_path.exists():
        return None
    
    try:
        with config_path.open() as f:
            cfg = yaml.safe_load(f)
        ligands_csv = Path(cfg["paths"]["ligands_folder"]) / "ligands.csv"
        if not ligands_csv.exists():
            return None
        
        # Read CSV with proper encoding and delimiter detection
        df = pd.read_csv(ligands_csv, sep=";", encoding="utf-8-sig")
        df.columns = [c.strip().lstrip('\ufeff') for c in df.columns]
        
        # Find ID and SMILES columns (case-insensitive)
        id_col = next((c for c in df.columns if c.lower() == "id"), None)
        smiles_col = next((c for c in df.columns if c.lower() == "smiles"), None)
        
        if not id_col or not smiles_col:
            return None
        
        # Find matching ligand
        match = df[df[id_col].astype(str).str.strip() == str(ligand_id).strip()]
        if not match.empty:
            smiles = match.iloc[0][smiles_col]
            if pd.notna(smiles) and str(smiles).strip():
                return str(smiles).strip()
    except Exception:
        # If anything fails, return None (fallback to PLIP SMILES)
        pass
    
    return None


def build_complex_assets(complex_dir: Path, config_path: Path | None = None) -> ComplexAssets:
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

    # Try to load original SMILES from CSV first (most accurate)
    # If not available, fall back to PLIP SMILES, then PDB file
    original_smiles = _load_smiles_from_csv(ligand_id, config_path)
    if original_smiles:
        # Use original SMILES from CSV (most reliable)
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            mol = Chem.MolFromSmiles(original_smiles)
            if mol:
                # Remove hydrogens for cleaner 2D visualization (like PDB loading does)
                mol = Chem.RemoveHs(mol)
                AllChem.Compute2DCoords(mol)
                ligand = mol
            else:
                # Fallback to PLIP/PDB loading
                ligand = _load_ligand(complex_dir, plip_xml=plip_xml if plip_xml.exists() else None)
        except Exception:
            # If SMILES parsing fails, fallback
            ligand = _load_ligand(complex_dir, plip_xml=plip_xml if plip_xml.exists() else None)
    else:
        # Load ligand using SMILES from PLIP if available (more accurate for 2D rendering)
        ligand = _load_ligand(complex_dir, plip_xml=plip_xml if plip_xml.exists() else None)
    
    # Build atom map: prefer PDB serial numbers if available, otherwise use SMILES mapping
    atom_map: dict[str, int] = {}
    
    # Check if ligand was loaded from SMILES (no PDB residue info) or from PDB (has residue info)
    has_pdb_info = any(atom.GetPDBResidueInfo() and atom.GetPDBResidueInfo().GetSerialNumber() 
                       for atom in ligand.GetAtoms())
    
    if has_pdb_info:
        # Ligand loaded from PDB - use PDB serial numbers
        for atom in ligand.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info and info.GetSerialNumber():
                atom_map[str(info.GetSerialNumber())] = atom.GetIdx()
            # Also map by index as fallback
            atom_map.setdefault(str(atom.GetIdx()), atom.GetIdx())
    else:
        # Ligand loaded from SMILES - use SMILES mapping from PLIP
        # smiles_map maps SMILES atom indices (from PLIP) to PDB atom indices (from docked structure)
        # We need to map these to molecule atom indices
        if smiles_map:
            # Create mapping: SMILES index -> molecule index
            # Since molecule from SMILES has atoms in SMILES order, we can use direct mapping
            # But we need to account for the fact that PLIP uses 1-based indices
            for smiles_idx, pdb_idx in smiles_map.items():
                try:
                    # Try to use SMILES index directly (convert to 0-based if needed)
                    smiles_int = int(smiles_idx.strip())
                    # PLIP SMILES indices are typically 1-based, molecule indices are 0-based
                    mol_idx = smiles_int - 1 if smiles_int > 0 else smiles_int
                    if 0 <= mol_idx < ligand.GetNumAtoms():
                        atom_map[smiles_idx] = mol_idx
                        atom_map[pdb_idx] = mol_idx  # Also map PDB index
                except (ValueError, IndexError):
                    pass
        
        # Always add index-based mapping as fallback
        for idx in range(ligand.GetNumAtoms()):
            atom_map.setdefault(str(idx), idx)
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

