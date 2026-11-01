"""Reconstruct supported specimen taxon-to-trait contributions and summaries."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from ac.cli import entry
if __name__ == "__main__":
    entry("18")
