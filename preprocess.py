import argparse
import copy
import json
import os

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def get_parser():
    parser = argparse.ArgumentParser(
        description="",
        usage=f"python {os.path.basename(__file__)}"
    )
    parser.add_argument(
        "-f", "--json_file", type=str, required=True,
        help="path to a ASCKOS json file"
    )
    parser.add_argument(
        "-d", "--debug", action='store_true',
        help="debug mode"
    )
    return parser.parse_args()


def create_molfile(mol: Chem.Mol, node_id: str, output_dir: str)->None:
    _mol = copy.deepcopy(mol)
    AllChem.EmbedMolecule(_mol)
    output_path = f"{output_dir}/{node_id}.txt"  # because .mol is not supported in TextAsset
    Chem.MolToMolFile(_mol, output_path)
    

def draw_2dmol(mol: Chem.Mol, node_id: str, output_dir: str) -> None:
    _mol = copy.deepcopy(mol)
    AllChem.Compute2DCoords(_mol)

    draw_opts = rdMolDraw2D.MolDrawOptions()
    draw_opts.clearBackground = False
    draw_opts.fixedFontSize = 60
    draw_opts.bondLineWidth = 3

    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
    drawer.SetDrawOptions(draw_opts)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    output_path = f"{output_dir}/{node_id}.png"
    with open(output_path, 'wb') as f:
        f.write(drawer.GetDrawingText())


def main():
    print("Preprocess ASKCOS JSON file...")
    args = get_parser()
    if not args.debug:
        RDLogger.DisableLog("rdApp.*")
    path_json = args.json_file
    with open(path_json) as f:
        data = json.load(f)

    output_parentdir = os.path.dirname(path_json)
    output_molfile_dir = os.path.join(output_parentdir, "molfiles")
    output_img_dir = os.path.join(output_parentdir, "images")

    for node in data['dispGraph']['nodes']:
        if node['type'] == 'reaction':
            continue
        mol = Chem.MolFromSmiles(node['smiles'])
        create_molfile(mol, node['id'], output_molfile_dir)
        draw_2dmol(mol, node['id'], output_img_dir)
    
    print(f"Mol files were saved in {output_molfile_dir}")
    print(f"PNG files were saved in {output_img_dir}")
    print("DONE!")


if __name__ == "__main__":
    main()
