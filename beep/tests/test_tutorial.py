import os
import ast
import unittest
"""
Automatically test the code shown in the tutorial. Is skipped if there is no
markdown source doc directory.
"""

this_dir = os.path.dirname(os.path.abspath(__file__))
docs_src_dir = os.path.join(
    this_dir,
    os.path.join(
        os.pardir,
        os.path.join(os.pardir, "docs_src")
    )
)
HAS_DOCS_SRC = os.path.isdir(docs_src_dir)


class DocumentationTutorialTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tutorial_src_path = os.path.join(docs_src_dir, "tutorial.md")

    @unittest.skipIf(not HAS_DOCS_SRC, "Docs directory not found, cannot test tutorial")
    def test_tutorial_valid(self):
        self.assertTrue(os.path.exists(self.tutorial_src_path))
        blocks = read_code_blocks_from_md(self.tutorial_src_path)

        is_valid = True
        try:
            for b in blocks:
                ast.parse(b)
        except SyntaxError:
            is_valid = False
        self.assertTrue(is_valid)

    @unittest.skipIf(not HAS_DOCS_SRC, "Docs directory not found, cannot test tutorial")
    def test_tutorial_code(self):
        blocks = read_code_blocks_from_md(self.tutorial_src_path)
        for b in blocks:
            if "plt.show()" in b:
                block_safe = b.replace("plt.show()", "plt.savefig('fname.png')")
            else:
                block_safe = b
            exec(block_safe)


def read_code_blocks_from_md(md_path):
    """
    Read ```python annotated code blocks from a markdown file.

    Args:
        md_path (str): Path to the markdown fle

    Returns:
        py_blocks ([str]): The blocks of python code.

    """
    with open(md_path, "r") as f:
        full_md = f.read()

    md_py_splits = full_md.split("```python")[1:]
    py_blocks = [split.split("```")[0] for split in md_py_splits]
    return py_blocks


if __name__ == "__main__":
    unittest.main()