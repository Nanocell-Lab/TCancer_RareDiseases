from __future__ import annotations

import unittest
from pathlib import Path

import rare_oncology


class IsolationTests(unittest.TestCase):
    def test_import_resolves_inside_21(self):
        package_file = Path(rare_oncology.__file__).resolve()
        self.assertIn("2.1", package_file.parts)
        self.assertNotIn("2.0", package_file.parts)

    def test_project_version_is_21(self):
        root = Path(__file__).resolve().parents[1]
        pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
        self.assertIn('version = "2.1.0"', pyproject)


if __name__ == "__main__":
    unittest.main()
