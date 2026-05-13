#!/usr/bin/env python3
"""
Route data extractor from pickle files.

This module deserializes Syntheseus route objects from pickle files and converts them
into JSON-serializable dictionaries containing node information (molecules, reactions, metadata).
"""

import json
import pickle
from typing import Any, Dict, Union


def extract_data(file_path: str) -> str:
    """
    Extract Mols, Reactions, and relevant information from a pickled route.

    Deserializes a pickled Syntheseus route object and extracts node information,
    including molecule SMILES, reaction SMILES, purchasability, and node depth.

    Args:
        file_path: Path to the pickle file containing a Syntheseus route object.

    Returns:
        JSON string representation of extracted node data.
        Keys are "node_1", "node_2", etc. (indexed by depth order).
        Each node contains:
            - depth: Node depth in the route graph
            - is_mol: Boolean indicating if node is a molecule
            - mol_smiles: SMILES string if molecule, None otherwise
            - is_rxn: Boolean indicating if node is a reaction
            - rxn_smiles: SMILES string if reaction, None otherwise
            - rxn_class: Reserved for reaction classification (currently None)
            - rxn_name: Reserved for reaction name (currently None)
            - is_purchasable: Purchasability metadata if molecule, None otherwise

    Raises:
        pickle.UnpicklingError: If the file cannot be unpickled.
        Exception: If node extraction fails.
    """
    # Load the pickle file
    with open(file_path, "rb") as f:
        route = pickle.load(f)  # set of syntheseus nodes

    syntheseus_route_data: Dict[str, Any] = {}
    # Sort nodes by depth
    sorted_nodes = sorted(enumerate(route), key=lambda x: x[1].depth)
    for idx, (_, node) in enumerate(sorted_nodes):
        syntheseus_route_data[f"node_{idx+1}"] = {
            "depth": node.depth,
            "is_mol": hasattr(node, "mol"),
            "mol_smiles": getattr(node, "mol", None).smiles if hasattr(node, "mol") else None,
            "is_rxn": hasattr(node, "reaction"),
            "rxn_smiles": str(getattr(node, "reaction", None)) if hasattr(node, "reaction") else None,
            "rxn_class": None,  # Dummy value
            "rxn_name": None,  # Dummy value
            "is_purchasable": node.mol.metadata["is_purchasable"] if hasattr(node, "mol") else None
        }

    # Convert the data to JSON format
    return json.dumps(syntheseus_route_data)


if __name__ == "__main__":
    import sys
    pickle_file = sys.argv[1]  # The path to the pickle file is passed as an argument
    data = extract_data(pickle_file)
    print(data)  # Capture result
