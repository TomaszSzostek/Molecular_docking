import gzip, shutil, requests
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select


# ── helpers ───────────────────────────────────────────
def _guess_biggest_het_ligand(pdb_file: Path) -> str:
    """Return 3-letter resname of the largest HETATM ligand."""
    sizes = {}
    with open(pdb_file) as fh:
        for line in fh:
            if line.startswith("HETATM"):
                res = line[17:20].strip()
                sizes[res] = sizes.get(res, 0) + 1
    return max(sizes, key=sizes.get)


def fetch(pdb_id: str, out_dir: Path) -> Path:
    """Download and unzip PDB file."""
    out_dir.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb.gz"
    gz_path = out_dir / f"{pdb_id}.pdb.gz"
    pdb_path = out_dir / f"{pdb_id}.pdb"

    r = requests.get(url, timeout=30)
    r.raise_for_status()
    gz_path.write_bytes(r.content)
    with gzip.open(gz_path, "rb") as gzf, open(pdb_path, "wb") as out:
        shutil.copyfileobj(gzf, out)
    gz_path.unlink()
    return pdb_path

# ── selectors ────────────────────────────────────────

class LigandSelect(Select):
    def __init__(self, het_id):
        self.het_id = het_id.strip().upper()
        self.found = 0

    def accept_residue(self, residue):
        resname = residue.resname.strip().upper()
        is_lig = resname == self.het_id
        if is_lig:
            self.found += 1
        return is_lig


class ProteinSelectOneChain(Select):
    def __init__(self, ligand_resname, chain_id="A"):
        self.ligand_resname = ligand_resname.strip().upper()
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        hetfield = residue.id[0]
        resname = residue.resname.strip().upper()
        return (
            hetfield == " "
            and resname != self.ligand_resname
            and resname not in {"HOH", "WAT"}
        )

# ── main splitter ─────────────────────────────────────

def split_receptor_ligand(pdb_path: Path, ligand_resname: str, receptors_dir: Path, ligands_dir: Path):
    """Split a PDB file into receptor (chain A only) and ligand (by resname)."""
    struct = PDBParser(QUIET=True).get_structure('X', str(pdb_path))

    pdb_id = pdb_path.stem
    rec_path = receptors_dir / f"{pdb_id}.pdb"
    lig_path = ligands_dir / f"{pdb_id}_{ligand_resname}.pdb"

    io = PDBIO()

    # Receptor: only chain A, no ligand/water
    io.set_structure(struct)
    io.save(str(rec_path), ProteinSelectOneChain(ligand_resname))

    # Ligand
    io.set_structure(struct)
    lig_sel = LigandSelect(ligand_resname)
    io.save(str(lig_path), lig_sel)

    if lig_sel.found == 0:
        raise RuntimeError(f"❌ Ligand '{ligand_resname}' not found in {pdb_path.name} (check HETATM lines)")

    for extra in receptors_dir.glob(f"{pdb_id}_*.pdb"):
        extra.unlink(missing_ok=True)


# ── batch driven by config ───────────────────────────

def fetch_and_split_batch(cfg, log):
    pdb_ids = cfg["fetch_crystals"]["pdb_ids"]
    res_map = cfg["fetch_crystals"].get("ligand_resnames", {})
    raw_dir = Path(cfg["paths"]["crystals_folder"]).expanduser()

    receptors_dir = Path(cfg["paths"]["receptors_folder"]).expanduser()
    ligands_dir = Path(cfg["paths"]["native_ligands_folder"]).expanduser()
    receptors_dir.mkdir(parents=True, exist_ok=True)
    ligands_dir.mkdir(parents=True, exist_ok=True)

    for pdb_id in pdb_ids:
        log.info(f"Fetch {pdb_id} …")
        pdb_path = fetch(pdb_id, raw_dir)

        ligand_resname = res_map.get(pdb_id)
        if not ligand_resname:
            ligand_resname = _guess_biggest_het_ligand(pdb_path)
            log.warning(f"No ligand_resname for {pdb_id}, guessed: {ligand_resname}")

        log.info(f"Split {pdb_id} → receptor + {ligand_resname}")
        split_receptor_ligand(
            pdb_path,
            ligand_resname,
            receptors_dir,
            ligands_dir
        )
