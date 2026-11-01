"""Run the included original fastp/Bowtie2/Kraken2/Bracken adapter."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from ac.cli import entry
if __name__ == "__main__":
    entry("05")
