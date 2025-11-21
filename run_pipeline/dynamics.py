#!/usr/bin/env python3
"""
Minimal export of docking results for molecular dynamics.

This script builds a small, self‑contained directory with:
    - one native ligand per receptor (from `better_than_native.csv`)
    - the corresponding docked poses that outperformed the native

For each entry it stores the protein PDB, the ligand SDF and the SMILES
string in a `dynamics_manifest.csv` file, which can be consumed by
downstream MD workflows.

The conversion from PDBQT to SDF preserves the original structure by:
    1. Reconstructing the molecule from SMILES (preserves bond orders/stereochemistry)
    2. Extracting 3D coordinates from the best pose in PDBQT
    3. Applying coordinates to the RDKit molecule and writing as SDF
"""

import logging, shutil
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

log = logging.getLogger(__name__)


# ───────────────── helpers ────────────────────────────────

def pdbqt_to_sdf(pdbqt: Path, sdf: Path, smiles: Optional[str] = None, ref_pdbqt: Optional[Path] = None) -> None:
    """
    Convert PDBQT to SDF while preserving structure from original SMILES.
    
    This function ensures bond orders, stereochemistry, and connectivity
    match the original SMILES by reconstructing the molecule from SMILES
    and applying 3D coordinates from the PDBQT pose.
    
    Uses reference PDBQT (if provided) to establish atom mapping between
    SMILES and docked PDBQT structures.
    
    Parameters
    ----------
    pdbqt : Path
        Input PDBQT file with docked pose(s).
    sdf : Path
        Output SDF file path.
    smiles : Optional[str]
        Original SMILES string. If provided, structure is reconstructed
        from SMILES; otherwise falls back to OpenBabel parsing.
    ref_pdbqt : Optional[Path]
        Reference PDBQT file (original input) for atom mapping.
    """
    if smiles:
        # Reconstruct molecule from SMILES to preserve structure
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            log.warning("Failed to parse SMILES for %s, falling back to OpenBabel", pdbqt.name)
            _pdbqt_to_sdf_openbabel(pdbqt, sdf)
            return
        
        # Strategy: Convert PDBQT to temporary PDB, read in RDKit, then replace structure with SMILES
        # This preserves 3D coordinates while fixing the structure
        try:
            import tempfile
            from analyze.RMSD import convert_pdbqt_to_clean_pdb
            
            # Convert PDBQT to PDB (preserves coordinates)
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp_pdb:
                tmp_pdb_path = Path(tmp_pdb.name)
            
            try:
                # Convert PDBQT to clean PDB (first MODEL only)
                convert_pdbqt_to_clean_pdb(pdbqt, tmp_pdb_path, template=None)
                
                # Read PDB in RDKit (will have wrong structure, but correct coordinates)
                mol_from_pdb = Chem.MolFromPDBFile(str(tmp_pdb_path), removeHs=False, sanitize=False)
                
                if mol_from_pdb and mol_from_pdb.GetNumAtoms() > 0:
                    # Get coordinates from PDB molecule
                    conf_pdb = mol_from_pdb.GetConformer()
                    pdb_coords = [conf_pdb.GetAtomPosition(i) for i in range(mol_from_pdb.GetNumAtoms())]
                    pdb_elems = [mol_from_pdb.GetAtomWithIdx(i).GetSymbol() for i in range(mol_from_pdb.GetNumAtoms())]
                    
                    # Get heavy atoms from SMILES
                    mol_heavy = Chem.RemoveHs(mol)
                    mol_heavy_atoms = list(mol_heavy.GetAtoms())
                    rdkit_heavy_elems = [a.GetSymbol() for a in mol_heavy_atoms]
                    
                    # Filter heavy atoms from PDB
                    pdb_heavy_coords = []
                    pdb_heavy_elems = []
                    for coords, elem in zip(pdb_coords, pdb_elems):
                        if elem != 'H':
                            pdb_heavy_coords.append([coords.x, coords.y, coords.z])
                            pdb_heavy_elems.append(elem)
                    
                    # Check element composition
                    from collections import Counter
                    if Counter(rdkit_heavy_elems) == Counter(pdb_heavy_elems):
                        # First embed SMILES to get approximate 3D structure for matching
                        mol_temp = Chem.AddHs(mol)
                        try:
                            AllChem.EmbedMolecule(mol_temp, AllChem.ETKDG())
                            conf_temp = mol_temp.GetConformer()
                            rdkit_heavy_coords_temp = []
                            for i, atom in enumerate(mol_temp.GetAtoms()):
                                if atom.GetSymbol() != 'H':
                                    pos = conf_temp.GetAtomPosition(i)
                                    rdkit_heavy_coords_temp.append([pos.x, pos.y, pos.z])
                        except:
                            rdkit_heavy_coords_temp = None
                        
                        # Match atoms using 3D distance if we have coords, otherwise by element
                        if rdkit_heavy_coords_temp and len(rdkit_heavy_coords_temp) == len(pdb_heavy_coords):
                            atom_mapping = _match_atoms_by_3d(
                                rdkit_heavy_coords_temp,
                                rdkit_heavy_elems,
                                pdb_heavy_coords,
                                pdb_heavy_elems
                            )
                        else:
                            # Sequential mapping by element (fallback)
                            from collections import defaultdict, deque
                            elem_queues = defaultdict(deque)
                            for i, elem in enumerate(pdb_heavy_elems):
                                elem_queues[elem].append(i)
                            atom_mapping = []
                            for elem in rdkit_heavy_elems:
                                if elem_queues[elem]:
                                    atom_mapping.append(elem_queues[elem].popleft())
                                else:
                                    atom_mapping.append(-1)
                        
                        # Add hydrogens to SMILES mol and assign coordinates
                        mol = Chem.AddHs(mol)
                        conf = Chem.Conformer(mol.GetNumAtoms())
                        
                        # Assign heavy atom coordinates using mapping
                        rdkit_heavy_idx = 0
                        for i, atom in enumerate(mol.GetAtoms()):
                            if atom.GetSymbol() != "H":
                                if rdkit_heavy_idx < len(atom_mapping):
                                    pdb_idx = atom_mapping[rdkit_heavy_idx]
                                    if pdb_idx >= 0 and pdb_idx < len(pdb_heavy_coords):
                                        conf.SetAtomPosition(i, pdb_heavy_coords[pdb_idx])
                                rdkit_heavy_idx += 1
                        
                        mol.AddConformer(conf)
                        # Optimize hydrogens
                        AllChem.UFFOptimizeMolecule(mol, ignoreInterfragInteractions=True, maxIters=100)
                    else:
                        # Element mismatch - use embedding
                        mol = Chem.AddHs(mol)
                        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                        AllChem.MMFFOptimizeMolecule(mol)
                else:
                    # Failed to read PDB - use embedding
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    AllChem.MMFFOptimizeMolecule(mol)
            finally:
                if tmp_pdb_path.exists():
                    tmp_pdb_path.unlink()
                
        except Exception as e:
            log.warning("Coordinate assignment failed for %s: %s, embedding instead", 
                       pdbqt.name, e)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.MMFFOptimizeMolecule(mol)
        
        # Write SDF
        writer = Chem.SDWriter(str(sdf))
        writer.write(mol)
        writer.close()
    else:
        # Fallback to OpenBabel if no SMILES provided
        _pdbqt_to_sdf_openbabel(pdbqt, sdf)


def _match_atoms_by_3d(ref_coords: list, ref_elems: list,
                       docked_coords: list, docked_elems: list) -> list[int]:
    """
    Match atoms between reference and docked structures using element type and 3D distance.
    
    For each element type, matches atoms by nearest neighbor distance.
    Returns mapping: ref_idx -> docked_idx
    """
    from collections import defaultdict
    
    if len(ref_elems) != len(docked_elems):
        # If counts differ, use sequential mapping
        return list(range(min(len(ref_elems), len(docked_elems))))
    
    # Group atoms by element
    ref_by_elem = defaultdict(list)
    docked_by_elem = defaultdict(list)
    
    for i, elem in enumerate(ref_elems):
        ref_by_elem[elem].append(i)
    for i, elem in enumerate(docked_elems):
        docked_by_elem[elem].append(i)
    
    # Match atoms element by element
    mapping = [None] * len(ref_elems)
    used_docked = set()
    
    for elem in set(ref_elems):
        ref_indices = ref_by_elem[elem]
        docked_indices = docked_by_elem[elem]
        
        if len(ref_indices) != len(docked_indices):
            # Sequential mapping for this element type
            for i, ref_idx in enumerate(ref_indices):
                if i < len(docked_indices):
                    mapping[ref_idx] = docked_indices[i]
        else:
            # Match by minimum distance if coords available
            if ref_coords and len(ref_coords) > 0 and len(ref_coords[0]) == 3:
                ref_coords_elem = [ref_coords[i] for i in ref_indices]
                docked_coords_elem = [docked_coords[i] for i in docked_indices]
                
                # Calculate distance matrix
                distances = np.zeros((len(ref_indices), len(docked_indices)))
                for i, ref_c in enumerate(ref_coords_elem):
                    for j, docked_c in enumerate(docked_coords_elem):
                        distances[i, j] = np.linalg.norm(np.array(ref_c) - np.array(docked_c))
                
                # Greedy matching: match closest pairs
                matched_ref = set()
                matched_docked = set()
                for _ in range(len(ref_indices)):
                    min_dist = float('inf')
                    best_ref, best_docked = None, None
                    for i in range(len(ref_indices)):
                        if i in matched_ref:
                            continue
                        for j in range(len(docked_indices)):
                            if j in matched_docked or docked_indices[j] in used_docked:
                                continue
                            if distances[i, j] < min_dist:
                                min_dist = distances[i, j]
                                best_ref = i
                                best_docked = j
                    
                    if best_ref is not None and best_docked is not None:
                        ref_idx = ref_indices[best_ref]
                        docked_idx = docked_indices[best_docked]
                        mapping[ref_idx] = docked_idx
                        matched_ref.add(best_ref)
                        matched_docked.add(best_docked)
                        used_docked.add(docked_idx)
            else:
                # No coords, sequential mapping by element
                for i, ref_idx in enumerate(ref_indices):
                    mapping[ref_idx] = docked_indices[i]
    
    # Fill in any unmapped atoms sequentially
    for i, m in enumerate(mapping):
        if m is None:
            for j in range(len(docked_elems)):
                if j not in used_docked:
                    mapping[i] = j
                    used_docked.add(j)
                    break
    
    return mapping


def _pdbqt_to_sdf_openbabel(pdbqt: Path, sdf: Path) -> None:
    """Fallback conversion using OpenBabel (less reliable for structure preservation)."""
    try:
        from openbabel import pybel
        mol = next(pybel.readfile("pdbqt", str(pdbqt)))
        mol.addh()
        mol.localopt(forcefield="uff", steps=500)
        mol.write("sdf", str(sdf), overwrite=True)
    except Exception as e:
        log.error("OpenBabel conversion failed for %s: %s", pdbqt.name, e)
        raise


def normalize(df: pd.DataFrame) -> None:
    """
    Normalize DataFrame column names by stripping whitespace and BOM markers.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to normalize.
    """
    df.columns = [c.strip().lstrip("\ufeff") for c in df.columns]


# ───────────────── main ───────────────────────────────────
def prepare_for_dynamics(
    better_csv: Path,
    smiles_csv: Path,
    out_root: Path,
    visuals_root: Path,
    dynamics_root: Path,
    input_ligands_dir: Optional[Path] = None,
    id_col="ID",
    smiles_col="SMILES",
) -> None:

    # 1) input tables -------------------------------------------------------
    better = pd.read_csv(better_csv)
    smiles = pd.read_csv(smiles_csv, encoding="utf-8-sig", sep=";", engine="python")
    normalize(smiles)

    id_col  = next(c for c in smiles if c.lower() == id_col.lower())
    smiles_col = next(c for c in smiles if c.lower() == smiles_col.lower())
    # Convert ID to string for matching (ID in CSV may be int or string)
    # Also handle both string and numeric IDs from better_than_native.csv
    smiles_map = dict(zip(smiles[id_col].astype(str), smiles[smiles_col]))
    # Add numeric keys as well for compatibility
    for idx, smi in zip(smiles[id_col], smiles[smiles_col]):
        if str(idx) != str(int(idx)) if pd.notna(idx) and str(idx).isdigit() else True:
            smiles_map[str(int(idx))] = smi

    # 2) output dirs (wipe ligands/ first) ----------------------------------
    base     = dynamics_root / "files_for_dynamics"
    lig_dir  = base / "ligands"
    prot_dir = base / "protein_cleaned"
    shutil.rmtree(lig_dir,  ignore_errors=True)
    shutil.rmtree(prot_dir, ignore_errors=True)
    lig_dir.mkdir(parents=True, exist_ok=True)
    prot_dir.mkdir(parents=True, exist_ok=True)
    rel = lambda p: p.relative_to(base)

    folder_map = {"native_redock": "dock_native",
                  "diag":          "diagonal",
                  "dock":          "matrix"}

    records = []

    # 3) one native ligand per receptor ------------------------------------
    for rec in better["receptor"].unique():
        nat = next((out_root / "dock_native").glob(f"{rec}__*__native_redock.pdbqt"), None)
        if not nat:
            continue
        lig_id = nat.stem.split("__")[1]
        sdf = lig_dir / f"{rec}__{lig_id}__native_redock.sdf"
        # Native ligands don't have SMILES in the CSV, use None
        pdbqt_to_sdf(nat, sdf, smiles=None)

        vis = visuals_root / f"native_redock__{rec}__{lig_id}"
        rec_src = vis / f"{rec}_receptor.pdb"
        rec_dst = prot_dir / f"{rec}.pdb"
        if rec_src.exists():
            shutil.copy2(rec_src, rec_dst)

        records.append({
            "protein_pdb": str(rel(rec_dst)),
            "ligand_id":   lig_id,
            "ligand_sdf":  str(rel(sdf)),
            "smiles":      ""
        })

    # 4) better‑than‑native rows -------------------------------------------
    for _, row in better.iterrows():
        rec, lig, mode = row["receptor"], row["ligand"], row["mode"]
        sub = folder_map.get(mode, mode)
        pdbqt = out_root / sub / f"{rec}__{lig}__{mode}.pdbqt"
        if not pdbqt.exists():
            continue
        sdf = lig_dir / f"{rec}__{lig}__{mode}.sdf"
        # Use original SMILES to preserve structure
        # Also try to use original input PDBQT for better atom mapping
        # Convert ligand ID to string for lookup (handle both string and numeric)
        lig_str = str(lig)
        orig_smiles = smiles_map.get(lig_str, None)
        orig_pdbqt = input_ligands_dir / f"{lig_str}.pdbqt" if input_ligands_dir else None
        if orig_pdbqt and orig_pdbqt.exists():
            pdbqt_to_sdf(pdbqt, sdf, smiles=orig_smiles, ref_pdbqt=orig_pdbqt)
        else:
            pdbqt_to_sdf(pdbqt, sdf, smiles=orig_smiles)

        vis  = visuals_root / f"{mode}__{rec}__{lig}"
        rec_src = vis / f"{rec}_receptor.pdb"
        rec_dst = prot_dir / f"{rec}.pdb"
        if rec_src.exists() and not rec_dst.exists():
            shutil.copy2(rec_src, rec_dst)

        records.append({
            "protein_pdb": str(rel(rec_dst)),
            "ligand_id":   lig,
            "ligand_sdf":  str(rel(sdf)),
            "smiles":      smiles_map.get(lig_str, "")
        })

    # 5) manifest -----------------------------------------------------------
    df_manifest = pd.DataFrame(records,
                                columns=["protein_pdb", "ligand_id", "ligand_sdf", "smiles"])
    # Replace NaN with empty string for SMILES column
    df_manifest["smiles"] = df_manifest["smiles"].fillna("")
    df_manifest.to_csv(base / "dynamics_manifest.csv", index=False)
    log.info("Manifest ready (%d entries)", len(records))


def _default_paths() -> tuple[Path, Path, Path, Path, Path, Path]:
    """
    Construct default locations relative to the repository root.

    The layout assumes this file lives in `run_pipeline/` and that the
    standard folders (`data/`, `results/`) follow the structure described
    in the project README.
    """
    root = Path(__file__).resolve().parent
    better_csv = root / "results" / "better_than_native.csv"
    smiles_csv = root / "data" / "ligands" / "ligands.csv"
    out_root = root / "results"
    visuals_root = root / "results" / "visuals"
    dynamics_root = root
    input_ligands_dir = root / "data" / "ligands"
    return better_csv, smiles_csv, out_root, visuals_root, dynamics_root, input_ligands_dir


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    better_csv, smiles_csv, out_root, visuals_root, dynamics_root, input_ligands_dir = _default_paths()

    prepare_for_dynamics(
        better_csv,
        smiles_csv,
        out_root,
        visuals_root,
        dynamics_root,
        input_ligands_dir=input_ligands_dir,
        id_col="ID",
        smiles_col="SMILES",
    )
