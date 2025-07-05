"""Tests for main.py module."""

from pathlib import Path


class TestMain:
    """Test main.py functionality."""

    def test_main_file_exists(self):
        """Test that main.py file exists."""
        main_path = Path("main.py")
        assert main_path.exists()

    def test_main_file_has_correct_structure(self):
        """Test that main.py has the expected structure."""
        main_path = Path("main.py")
        content = main_path.read_text()

        # Check for expected content
        assert "from src.ui.cli.app import main" in content
        assert 'if __name__ == "__main__"' in content
        assert "main()" in content

    def test_main_is_executable_script(self):
        """Test that main.py is properly structured as executable script."""
        main_path = Path("main.py")
        content = main_path.read_text()

        # Should have shebang
        lines = content.strip().split("\n")
        assert lines[0].startswith("#!")

        # Should have docstring
        assert '"""' in content
