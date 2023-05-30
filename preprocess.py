import argparse
import copy
import glob
import json
import os

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw
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
    parser.add_argument(
        "--delete", action='store_true',
        help="delete png and txt(mol) files related to the json file"
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
    draw_opts.clearBackground = True
    draw_opts.fixedFontSize = 55
    draw_opts.bondLineWidth = 3

    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
    drawer.SetDrawOptions(draw_opts)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    output_path = f"{output_dir}/{node_id}.png"
    with open(output_path, 'wb') as f:
        f.write(drawer.GetDrawingText())


def draw_rxn_image(rxn_smiles: str, rxn_id: str, output_dir: str) -> None:
    rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
    
    draw_opts = rdMolDraw2D.MolDrawOptions()
    draw_opts.clearBackground = True
    draw_opts.fixedFontSize = 30
    draw_opts.bondLineWidth = 2.5
    
    drawer = rdMolDraw2D.MolDraw2DCairo(600, 300)
    drawer.SetDrawOptions(draw_opts)
    drawer.DrawReaction(rxn)
    drawer.FinishDrawing()

    output_path = f"{output_dir}/{rxn_id}.png"
    with open(output_path, 'wb') as f:
        f.write(drawer.GetDrawingText())


def delete_all_files(directory: str) -> None:
    files = glob.glob(f"{directory}/*")
    for f in files:
        if os.path.isfile(f):
            os.remove(f)


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

    if args.delete:
        print(f"Mol files were deleted: {output_molfile_dir}")
        delete_all_files(output_molfile_dir)
        print(f"PNG files were deleted: {output_img_dir}")
        delete_all_files(output_img_dir)
        return

    for node in data['dispGraph']['nodes']:
        if node['type'] == 'reaction':
            draw_rxn_image(node['smiles'], node['id'], output_img_dir)
        else:  # 'chemical'
            mol = Chem.MolFromSmiles(node['smiles'])
            create_molfile(mol, node['id'], output_molfile_dir)
            draw_2dmol(mol, node['id'], output_img_dir)
    
    print(f"Mol files were saved in {output_molfile_dir}")
    print(f"PNG files were saved in {output_img_dir}")
    print("DONE!")


if __name__ == "__main__":
    main()
