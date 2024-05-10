from rdkit import Chem
from rdkit.Chem import AllChem
import logging as log
from vina import Vina
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def optimize_conformation(mol):
    # Votre fonction d'optimisation de conformation ici
    pass

def load(pdb_file):
    atom_coordinates = []

    # Ouvrir le fichier PDB en mode lecture
    with open(pdb_file, 'r') as file:
        # Parcourir chaque ligne du fichier
        for line in file:
            # Vérifier si la ligne commence par "ATOM" (les lignes contenant les coordonnées atomiques)
            if line.startswith('ATOM'):
                # Extraire les coordonnées x, y, z des atomes
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                # Ajouter les coordonnées à la liste des coordonnées des atomes
                atom_coordinates.append((x, y, z))

    return atom_coordinates

def get_extent(atom_coordinates):
    # Initialiser les coordonnées minimales et maximales avec les premières coordonnées atomiques
    min_x, min_y, min_z = atom_coordinates[0]
    max_x, max_y, max_z = atom_coordinates[0]

    # Parcourir les coordonnées de tous les atomes pour trouver les coordonnées minimales et maximales dans chaque dimension
    for x, y, z in atom_coordinates:
        min_x = min(min_x, x)
        min_y = min(min_y, y)
        min_z = min(min_z, z)
        max_x = max(max_x, x)
        max_y = max(max_y, y)
        max_z = max(max_z, z)

    # Retourner les coordonnées minimales et maximales sous forme de tuples
    min_coordinates = (min_x, min_y, min_z)
    max_coordinates = (max_x, max_y, max_z)
    return min_coordinates, max_coordinates

def get_box(pdb_file, extending=5.0):
    atom_coordinates = load(pdb_file)

    min_x, min_y, min_z = atom_coordinates[0]
    max_x, max_y, max_z = atom_coordinates[0]
    for x, y, z in atom_coordinates[1:]:
        min_x = min(min_x, x)
        min_y = min(min_y, y)
        min_z = min(min_z, z)
        max_x = max(max_x, x)
        max_y = max(max_y, y)
        max_z = max(max_z, z)

    minX = min_x - float(extending)
    minY = min_y - float(extending)
    minZ = min_z - float(extending)
    maxX = max_x + float(extending)
    maxY = max_y + float(extending)
    maxZ = max_z + float(extending)

    SizeX = maxX - minX
    SizeY = maxY - minY
    SizeZ = maxZ - minZ
    CenterX = (maxX + minX) / 2
    CenterY = (maxY + minY) / 2
    CenterZ = (maxZ + minZ) / 2

    return {'center_x': CenterX, 'center_y': CenterY, 'center_z': CenterZ}, \
           {'size_x': SizeX, 'size_y': SizeY, 'size_z': SizeZ}
    
def fix_protein(filename='',addHs_pH=7.4,output='',try_renumberResidues=False):

    fix = PDBFixer(filename=filename)
    fix.findMissingResidues()
    fix.findNonstandardResidues()
    fix.replaceNonstandardResidues()
    fix.removeHeterogens(True)
    fix.findMissingAtoms()
    fix.addMissingAtoms()
    fix.addMissingHydrogens(addHs_pH)
    PDBFile.writeFile(fix.topology, fix.positions, open(output, 'w'))

    if try_renumberResidues == True:
        try:
            original=mda.Universe(filename)
            from_fix=mda.Universe(output)

            resNum=[res.resid for res in original.residues]
            for idx,res in enumerate(from_fix.residues):
                res.resid = resNum[idx]

            save=PDB.PDBWriter(filename=output)
            save.write(from_fix)
            save.close()
        except Exception:
            print('Not possible to renumber residues, check excepton for extra details')
            

def dock_molecule(target_smile: str, protein_pdb: str, res_file: str, center_x: float, center_y: float, center_z: float,
                  size_x: float = 30, size_y: float = 30, size_z: float = 30, exhaustiveness: int = 8) -> float or None:
    mol = Chem.MolFromSmiles(target_smile)
    mol = optimize_conformation(mol)
    if mol is None:
        return None
    
    # Initialiser l'objet Vina
    v = Vina(sf_name='vina')
    v.set_receptor(protein_pdb)
    v.set_ligand_from_mol(mol)

    # Calculer les cartes Vina
    v.compute_vina_maps(center=[center_x, center_y, center_z], box_size=[size_x, size_y, size_z])

    # Effectuer le docking
    v.dock(exhaustiveness=exhaustiveness)

    # Obtenir le score de docking
    result = v.result
    docking_score = result.get_lowest_binding_energy()

    return docking_score
