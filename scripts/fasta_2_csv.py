import csv
import sys
from pathlib import Path

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna", ".faa", ".ffn", ".frn"}


def parse_fasta(fasta_path: Path) -> list[dict]:
    records = []
    accession_id = None
    sequence_parts = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if accession_id is not None:
                    records.append({
                        "Accession ID": accession_id,
                        "Sequence": "".join(sequence_parts),
                    })

                # EXtract Accession ID 
                header = line[1:]
                accession_id = header.split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)

        if accession_id is not None:
            records.append({
                "Accession ID": accession_id,
                "Sequence": "".join(sequence_parts),
            })

    return records


def convert_file(fasta_path: Path, output_dir: Path | None = None) -> Path:
    out_dir = output_dir if output_dir else fasta_path.parent
    output_path = out_dir / fasta_path.with_suffix(".csv").name

    records = parse_fasta(fasta_path)

    with open(output_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["Accession ID", "Sequence"])
        writer.writeheader()
        writer.writerows(records)

    print(f"  ✅ {fasta_path.name} → {output_path.name}  ({len(records)} sequences)")
    return output_path


def fasta_to_csv(input_path: str, output_dir: str | None = None) -> list[Path]:
    path = Path(input_path)
    out_dir = Path(output_dir) if output_dir else None

    if not path.exists():
        raise FileNotFoundError(f"Directory not found: {input_path}")

    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    # ── Single File ──────────────────────────────────────────────
    if path.is_file():
        if path.suffix.lower() not in FASTA_EXTENSIONS:
            raise ValueError(
                f"Extensão '{path.suffix}' não reconhecida como FASTA.\n"
                f"Extensões aceitas: {', '.join(sorted(FASTA_EXTENSIONS))}"
            )
        print(f"📄 Convertendo arquivo: {path}")
        return [convert_file(path, out_dir)]

    # ── Directory ──────────────────────────────────────────────────
    if path.is_dir():
        fasta_files = sorted(
            f for f in path.iterdir()
            if f.is_file() and f.suffix.lower() in FASTA_EXTENSIONS
        )

        if not fasta_files:
            print(f"FASTA file not found in: {path}")
            return []

        print(f"Directory: {path}  ({len(fasta_files)} files(s) found(s))\n")
        converted = []
        for fasta_file in fasta_files:
            converted.append(convert_file(fasta_file, out_dir))

        print(f"\n{len(converted)} File(s) converted(s).")
        return converted

    raise ValueError(f"'{input_path}' not valid.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso:")
        print("  python fasta_to_csv.py <arquivo.fasta>            #  single file")
        print("  python fasta_to_csv.py <diretorio/>               # directory")
        print("  python fasta_to_csv.py <entrada> <saida_dir/>     # define dir de saída")
        sys.exit(1)

    input_path  = sys.argv[1]
    output_dir  = sys.argv[2] if len(sys.argv) > 2 else None

    fasta_to_csv(input_path, output_dir)