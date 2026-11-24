"""Check path-only parameter preparation and safe refusal of existing output directories."""
import contextlib
import hashlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from prepare_paper_parameters import prepare, DATA


class ParameterCopies(unittest.TestCase):
    def test_only_paths_change(self):
        source_files = list(DATA.rglob('*'))
        before = {p: hashlib.sha256(p.read_bytes()).digest() for p in source_files if p.is_file()}
        with tempfile.TemporaryDirectory(prefix='hosc_parameters_') as folder:
            output = Path(folder) / 'runs'
            with contextlib.redirect_stdout(io.StringIO()):
                prepare(output)
            copies = list((output / 'Parameters').rglob('*.par'))
            self.assertEqual(len(copies), 138)
            for copy in copies:
                source = DATA / copy.relative_to(output)
                original = ET.parse(source).getroot()
                prepared = ET.parse(copy).getroot()
                self.assertEqual(len(original), len(prepared))
                for old, new in zip(original, prepared):
                    self.assertEqual(old.tag, new.tag)
                    if old.attrib != new.attrib:
                        self.assertIn(old.tag, ('audit_output_directory', 'origin_offsets_file'))
                        self.assertEqual(set(old.attrib), set(new.attrib))
                    if new.tag == 'audit_output_directory':
                        self.assertTrue(Path(new.get('value')).is_dir())
                    if new.tag == 'origin_offsets_file' and new.get('value'):
                        self.assertTrue(Path(new.get('value')).is_file())
            for source, target in [
                (DATA / 'Parameters/origin_offsets_500.csv', output / 'protocol/origin_offsets_500.csv'),
                (DATA / 'Parameters/Supp_weighted_attributes/Supp_samples.dat', output / 'Parameters/Supp_weighted_attributes/Supp_samples.dat')
            ]:
                self.assertEqual(source.read_bytes(), target.read_bytes())
        for path, digest in before.items():
            self.assertEqual(hashlib.sha256(path.read_bytes()).digest(), digest)

    def test_existing_directory_is_untouched(self):
        with tempfile.TemporaryDirectory(prefix='hosc_existing_') as folder:
            existing = Path(folder)
            with self.assertRaisesRegex(ValueError, 'existing files'):
                prepare(existing)
            self.assertEqual(list(existing.iterdir()), [])

    def test_non_ascii_path_is_rejected_before_creation(self):
        with tempfile.TemporaryDirectory(prefix='hosc_ascii_') as folder:
            destination = Path(folder) / '\u6d4b\u8bd5'
            with self.assertRaisesRegex(ValueError, 'ASCII'):
                prepare(destination)
            self.assertFalse(destination.exists())


if __name__ == '__main__':
    unittest.main()
