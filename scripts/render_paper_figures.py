import hashlib
import json
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs/corrections/scientific_reports_rpa"
OUT = ROOT / "docs/figures"
FIGURES = {
    "paper-sequences": "Fig1_sequence_construction",
    "paper-kappa-controls": "Fig2_kappa_controls",
    "paper-sequence-response": "Fig3_pi_kappa_response",
    "paper-rpa-diagnostic": "Fig4_sequence_response_theory",
    "paper-cross-attraction": "Fig5_cross_attraction",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    checksums = json.loads((SOURCE / "SHA256SUMS.json").read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    records = {}
    for name, figure in FIGURES.items():
        source = SOURCE / "manuscript_source/figures" / f"{figure}.pdf"
        relative = source.relative_to(SOURCE).as_posix()
        digest = sha(source)
        if digest != checksums[relative]:
            raise ValueError(f"Source checksum mismatch: {relative}")
        with tempfile.TemporaryDirectory(dir=OUT) as temporary:
            target = Path(temporary) / name
            svg = target.with_suffix(".svg")
            png = target.with_suffix(".png")
            subprocess.run(["pdftocairo", "-svg", str(source), str(svg)],
                           check=True, capture_output=True)
            subprocess.run(["pdftocairo", "-png", "-singlefile", "-r", "300",
                            str(source), str(target)], check=True, capture_output=True)
            ET.parse(svg)
            if png.read_bytes()[-12:] != b"\0\0\0\0IEND\xaeB`\x82":
                raise ValueError(f"Incomplete PNG export: {name}")
            svg.replace(OUT / svg.name)
            png.replace(OUT / png.name)
        records[name] = {
            "source": source.relative_to(ROOT).as_posix(),
            "source_sha256": digest,
            "svg_sha256": sha(OUT / f"{name}.svg"),
            "png_sha256": sha(OUT / f"{name}.png"),
        }
    provenance = {
        "method": "Unmodified manuscript PDF exports using pdftocairo; PNG at 300 dpi",
        "checksum_source": (SOURCE / "SHA256SUMS.json").relative_to(ROOT).as_posix(),
        "captions": (SOURCE / "manuscript_source/main.tex").relative_to(ROOT).as_posix(),
        "figures": records,
    }
    (OUT / "paper-source-verification.json").write_text(
        json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
