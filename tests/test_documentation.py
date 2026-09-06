import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"


class QuartoDocumentationTests(unittest.TestCase):
    def test_book_configuration_lists_existing_chapters(self):
        config = (DOCS / "_quarto.yml").read_text(encoding="utf-8")
        self.assertIn("type: book", config)
        self.assertIn("number-sections: false", config)
        self.assertIn('part: "User manual"', config)
        self.assertIn('part: "Physical theory"', config)
        self.assertIn('part: "Technical notes"', config)

        chapter_paths = re.findall(r"^\s*-\s+([^\s]+\.qmd)\s*$", config, re.MULTILINE)
        self.assertGreater(len(chapter_paths), 20)
        for chapter_path in chapter_paths:
            with self.subTest(chapter=chapter_path):
                self.assertTrue((DOCS / chapter_path).is_file())

    def test_local_links_and_explicit_anchors_resolve(self):
        link_pattern = re.compile(r"\[[^]]+\]\(([^)]+)\)")
        for source in DOCS.rglob("*.qmd"):
            text = source.read_text(encoding="utf-8")
            for target in link_pattern.findall(text):
                if target.startswith(("http://", "https://", "mailto:")):
                    continue
                path_text, separator, anchor = target.partition("#")
                target_path = source if not path_text else (source.parent / path_text)
                with self.subTest(source=source.name, target=target):
                    self.assertTrue(target_path.resolve().is_file())
                    if separator:
                        target_text = target_path.read_text(encoding="utf-8")
                        self.assertIn(f"{{#{anchor}}}", target_text)

    def test_bibliography_contains_every_citation_key(self):
        bibliography = (DOCS / "references.bib").read_text(encoding="utf-8")
        known_keys = set(re.findall(r"^@\w+\{([^,]+),", bibliography, re.MULTILINE))
        cited_keys = set()
        for source in DOCS.rglob("*.qmd"):
            cited_keys.update(
                re.findall(r"(?<!\w)@([A-Za-z][\w:-]+)", source.read_text(encoding="utf-8"))
            )
        self.assertTrue(cited_keys)
        self.assertEqual(cited_keys - known_keys, set())


if __name__ == "__main__":
    unittest.main()
