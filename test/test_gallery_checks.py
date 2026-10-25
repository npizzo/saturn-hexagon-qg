#!/usr/bin/env python3
"""Small local checks; no Git operations, credentials, network, or big fixtures."""

import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


spec = importlib.util.spec_from_file_location('gallery_checks', Path(__file__).resolve().parents[1] / 'scripts/check_gallery.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


class GalleryChecks(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name) / 'gallery'
        self.root.mkdir()
        (self.root / 'index.html').write_text('<!doctype html><title>Gallery</title>')

    def test_static_assets_and_nested_additions_are_accepted(self):
        nested = self.root / 'experiments-01'
        nested.mkdir()
        (nested / 'index.html').write_text('<title>Results</title>')
        (nested / 'metadata.json').write_text('{"status": "illustrative"}')
        (self.root / '.nojekyll').touch()
        count, size = checker.check_gallery(self.root)
        self.assertEqual(count, 4)
        self.assertGreater(size, 0)

    def test_raw_outputs_and_hidden_credentials_are_rejected(self):
        for name in ('checkpoint.jld2', 'snapshot.h5', 'archive.tar.gz', '.env', 'credentials.json'):
            with self.subTest(name=name):
                path = self.root / name
                path.write_bytes(b'not public')
                with self.assertRaises(ValueError):
                    checker.check_gallery(self.root)
                path.unlink()

    def test_symbolic_directory_links_and_hard_links_are_rejected(self):
        link = self.root / 'linked.html'
        link.symlink_to(self.root / 'index.html')
        with self.assertRaisesRegex(ValueError, 'symbolic link'):
            checker.check_gallery(self.root)
        link.unlink()
        link.symlink_to(self.root.parent, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'symbolic link'):
            checker.check_gallery(self.root)
        link.unlink()
        os.link(str(self.root / 'index.html'), str(link))
        with self.assertRaisesRegex(ValueError, 'hard links'):
            checker.check_gallery(self.root)

    def test_common_secret_markers_are_rejected_without_printing_values(self):
        path = self.root / 'metadata.json'
        for secret in (b'-----BEGIN OPENSSH PRIVATE KEY-----', b'github_pat_' + b'a' * 40,
                       b'ghp_' + b'b' * 36, b'AKIA' + b'A' * 16):
            path.write_bytes(secret)
            with self.assertRaisesRegex(ValueError, 'possible credential') as caught:
                checker.check_gallery(self.root)
            self.assertNotIn(secret.decode(), str(caught.exception))

    def test_size_limits_without_large_files(self):
        with patch.object(checker, 'FILE_LIMIT', 4):
            with self.assertRaisesRegex(ValueError, 'exceeds'):
                checker.check_gallery(self.root)
        with patch.object(checker, 'TOTAL_LIMIT', 4):
            with self.assertRaisesRegex(ValueError, 'smaller than'):
                checker.check_gallery(self.root)

    def test_missing_index_and_linked_root_are_rejected(self):
        (self.root / 'index.html').unlink()
        with self.assertRaisesRegex(ValueError, 'index.html'):
            checker.check_gallery(self.root)
        alias = self.root.parent / 'alias'
        alias.symlink_to(self.root, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'ordinary directory'):
            checker.check_gallery(alias)

    def test_renderer_toml_and_exact_completion_marker(self):
        nested = self.root / 'maintenance-01'
        nested.mkdir()
        (nested / 'config.toml').write_text('[grid]\nnx = 1024\n')
        marker = nested / 'COMPLETED'
        marker.write_bytes(b'Static gallery build completed.\n')
        self.assertEqual(checker.check_gallery(self.root)[0], 3)
        for contents in (b'1000000\n', b'Complete\n', b'Static gallery build completed.',
                         b'Static gallery build completed.\nextra data\n'):
            with self.subTest(contents=contents):
                marker.write_bytes(contents)
                with self.assertRaisesRegex(ValueError, 'exact static-gallery marker'):
                    checker.check_gallery(self.root)
        marker.unlink()
        (nested / 'OTHER_MARKER').write_bytes(b'Static gallery build completed.\n')
        with self.assertRaisesRegex(ValueError, 'unsupported static asset'):
            checker.check_gallery(self.root)


if __name__ == '__main__':
    unittest.main(verbosity=2)
