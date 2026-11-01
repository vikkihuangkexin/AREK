"""Reproduce continuous PCA and fixed-seed PAM stability without forcing K."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from ac.cli import entry
if __name__ == "__main__":
    entry("19")
