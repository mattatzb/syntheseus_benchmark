"""The generator must import in the pinned Python environment without network access."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from types import ModuleType
from unittest.mock import Mock, patch

import pandas as pd


class GeneratorCompatibilityTests(unittest.TestCase):
    def test_generator_import_and_unresolved_identifiers_need_no_network(self):
        # Substitute remote-service clients and the HTML parser so this import
        # check can also run in the existing analysis-only environment.
        cirpy = ModuleType("cirpy")
        cirpy.resolve = Mock(side_effect=AssertionError("Unexpected CIR lookup"))
        pubchempy = ModuleType("pubchempy")
        pubchempy.get_compounds = Mock(side_effect=AssertionError("Unexpected PubChem lookup"))
        bs4 = ModuleType("bs4")
        bs4.BeautifulSoup = object
        source = Path(__file__).resolve().parents[1] / "src/data_prep/generate_pesticide.py"
        spec = importlib.util.spec_from_file_location("generator_import_check", source)
        generator = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {"cirpy": cirpy, "pubchempy": pubchempy, "bs4": bs4}):
            spec.loader.exec_module(generator)
            self.assertEqual(generator.get_smiles_from_identifiers(None, None), "Failed")
        cirpy.resolve.assert_not_called()
        pubchempy.get_compounds.assert_not_called()

    def test_live_generator_preserves_cir_then_pubchem_fallback(self):
        # Exercise the refactored generator without scraping or resolving online.
        source = Path(__file__).resolve().parents[1] / "src/data_prep/generate_pesticide.py"
        spec = importlib.util.spec_from_file_location("generator_resolution_check", source)
        generator = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(generator)
        identifiers = pd.DataFrame({"name": ["cir-hit", "fallback"],
                                    "cas": ["cas-a", "cas-b"],
                                    "inchi_key": ["key-a", "key-b"]})
        session = Mock()
        session.__enter__ = Mock(return_value=session)
        session.__exit__ = Mock(return_value=False)
        with patch.object(generator, "_create_session", return_value=session), \
                patch.object(generator, "fetch_or_load_raw", return_value=identifiers), \
                patch.object(generator.cirpy, "resolve", side_effect=["CCO", None]) as cir, \
                patch.object(generator.pcp, "get_compounds",
                             return_value=[Mock(isomeric_smiles="CCN")]) as pubchem, \
                patch.object(pd.DataFrame, "to_csv", side_effect=AssertionError("Unexpected file write")):
            result = generator.generate_raw_dataset()
        self.assertEqual(result.smiles.tolist(), ["CCO", "CCN"])
        self.assertEqual([call.args for call in cir.call_args_list],
                         [("cas-a", "smiles"), ("cas-b", "smiles")])
        pubchem.assert_called_once_with("key-b", "inchikey")
        session.__exit__.assert_called_once()


if __name__ == "__main__":
    unittest.main()
