#!/usr/bin/env python3
"""
Dataset cleaning pipeline for GREAC FASTA datasets, built on top of SeqKit2.

Given a directory containing FASTA files (searched recursively), this script
creates a sibling directory `<input_dir>_clean` with the very same tree
structure and the very same file names, so any downstream GREAC/GRAMEP
workflow can be pointed at the clean copy without any other change.

Cleaning rules (all of them performed by SeqKit2 commands):

  (i)   seqkit rmdup   -> remove duplicated sequences inside each file;
  (ii)  seqkit common  -> find the sequences shared between the files of a
                          same group (the files living in the same directory),
                          compared pairwise;
  (iii) seqkit seq     -> extract those shared sequences into a reference file;
  (iv)  seqkit grep -v -> write the filtered files, without the shared
                          sequences, to avoid over-fitting.

The files are compared group by group: by default a group is the set of FASTA
files living in the same directory; with `--group all` the whole tree is a
single group, which is what the GREAC layout needs when every class has its
own subdirectory. Non-FASTA files (references, metadata, models) are copied
untouched, unless `--skip-extra` is given.

The reference files, the id lists and the report are kept inside the hidden
folder `.clean_meta/` of the output directory, so that they never show up in
`*.fasta` globs of the cleaned dataset.

Usage:
    python clean_dataset.py <fasta_dir> [options]

Examples:
    python clean_dataset.py ~/Desktop/datasets/virus_dataset
    python clean_dataset.py ./virus_dataset -o ./virus_clean --by seq -j 8
    python clean_dataset.py ./virus_dataset --group all --by seq -j 8
"""

from __future__ import annotations

import argparse
import itertools
import json
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna", ".faa", ".ffn", ".frn"}
META_DIRNAME = ".clean_meta"

# `seqkit common` / `seqkit rmdup` comparison criteria
BY_FLAGS = {
    "id": [],            # sequence ID (default, e.g. the accession)
    "name": ["-n"],      # full header name
    "seq": ["-s"],       # sequence content
}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def log(msg: str) -> None:
    print(msg, flush=True)


def is_fasta(path: Path) -> bool:
    suffixes = [s.lower() for s in path.suffixes]
    if not suffixes:
        return False
    if suffixes[-1] in {".gz", ".xz", ".bz2", ".zst"}:
        suffixes = suffixes[:-1]
    return bool(suffixes) and suffixes[-1] in FASTA_EXTENSIONS


def resolve_seqkit(explicit: str | None) -> str:
    candidates = [explicit] if explicit else ["seqkit", "seqkit2"]
    for cand in candidates:
        if cand and shutil.which(cand):
            return cand
    sys.exit(
        "❌ SeqKit2 not found. Install it (https://bioinf.shenwei.me/seqkit/) "
        "or point to the binary with --seqkit /path/to/seqkit"
    )


class SeqKit:
    """Thin wrapper around the seqkit CLI."""

    def __init__(self, binary: str, threads: int, quiet: bool = True):
        self.binary = binary
        self.threads = threads
        self.quiet = quiet

    def run(self, args: list[str], stdout: Path | None = None) -> str:
        cmd = [self.binary, *args, "--threads", str(self.threads)]
        if self.quiet:
            cmd.append("--quiet")
        if stdout is not None:
            with open(stdout, "w") as fh:
                proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
            out = ""
        else:
            proc = subprocess.run(cmd, capture_output=True, text=True)
            out = proc.stdout
        if proc.returncode != 0:
            raise RuntimeError(
                f"seqkit failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr.strip()}"
            )
        return out

    def version(self) -> str:
        proc = subprocess.run([self.binary, "version"], capture_output=True, text=True)
        return proc.stdout.strip() or "unknown"

    def count(self, fasta: Path) -> int:
        out = subprocess.run(
            [self.binary, "stats", "-T", str(fasta)], capture_output=True, text=True
        ).stdout.splitlines()
        if len(out) < 2:
            return 0
        try:
            return int(out[1].split("\t")[3])
        except (IndexError, ValueError):
            return 0


# --------------------------------------------------------------------------- #
# pipeline steps
# --------------------------------------------------------------------------- #
def step_rmdup(sk: SeqKit, src: Path, dst: Path, by: str) -> None:
    """(i) remove duplicated sequences inside a single file."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    sk.run(["rmdup", *BY_FLAGS[by], str(src)], stdout=dst)


def file_keys(sk: SeqKit, fasta: Path, by: str) -> list[tuple[str, str]]:
    """(id, comparison key) of every record of a file, listed by `seqkit seq`.

    The key is what `seqkit common` compares (the id, the full name or the
    sequence itself), while the id is what `seqkit grep` uses to filter, so a
    match found by sequence is still removed from the right records of every
    file, even when their ids differ.
    """
    if by == "seq":
        # sequences are only materialised when the comparison really needs them
        table = sk.run(["fx2tab", str(fasta)])
        pairs = []
        for line in table.splitlines():
            if not line.strip():
                continue
            name, _, seq = line.partition("\t")
            pairs.append((name.split()[0], seq.strip().upper()))
        return pairs

    flags = ["-n"] if by == "name" else ["-n", "-i"]
    out = sk.run(["seq", *flags, str(fasta)])
    return [(l.split()[0], l.strip()) for l in out.splitlines() if l.strip()]


def step_common_keys(sk: SeqKit, files: list[Path], by: str, workdir: Path) -> set[str]:
    """(ii) keys of the sequences shared between the files of a group.

    `seqkit common` returns the sequences present in *all* the given files, so
    the files are compared pairwise and the results are joined: any sequence
    seen in more than one file of the group is considered common.
    """
    common: set[str] = set()
    for idx, (a, b) in enumerate(itertools.combinations(files, 2)):
        out = workdir / f"pair_{idx}.fasta"
        sk.run(["common", *BY_FLAGS[by], str(a), str(b)], stdout=out)
        if out.stat().st_size == 0:
            continue
        # (iii) the shared records are read back with `seqkit seq`
        common.update(key for _, key in file_keys(sk, out, by))
        out.unlink(missing_ok=True)
    return common


def write_ids(ids: set[str], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(sorted(ids)) + "\n" if ids else "")
    return path


def step_extract_reference(sk: SeqKit, ids_files: list[tuple[Path, Path]], ref: Path) -> None:
    """(iii) extract the common sequences into a single reference file."""
    ref.parent.mkdir(parents=True, exist_ok=True)
    tmp = ref.with_suffix(".tmp.fasta")
    with open(tmp, "w") as fh:
        for fasta, ids_file in ids_files:
            proc = subprocess.run(
                [sk.binary, "grep", "-f", str(ids_file), str(fasta), "--quiet"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            )
            if proc.returncode != 0:
                raise RuntimeError(f"seqkit grep failed on {fasta}:\n{proc.stderr.strip()}")
            fh.write(proc.stdout)
    # deduplicate the reference itself, keeping one copy of each shared sequence
    sk.run(["rmdup", "-s", str(tmp)], stdout=ref)
    tmp.unlink(missing_ok=True)


def step_filter(sk: SeqKit, src: Path, dst: Path, ids_file: Path) -> None:
    """(iv) write the file without the common sequences."""
    sk.run(["grep", "-v", "-f", str(ids_file), str(src)], stdout=dst)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Clean a FASTA dataset with SeqKit2, mirroring the directory tree.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("input_dir", type=Path, help="directory containing the FASTA files")
    p.add_argument(
        "-o", "--output-dir", type=Path, default=None,
        help="output directory (default: <input_dir>_clean)",
    )
    p.add_argument(
        "--by", choices=sorted(BY_FLAGS), default="id",
        help="criterion used by rmdup/common: id (default), name or seq",
    )
    p.add_argument(
        "--rmdup-by", choices=sorted(BY_FLAGS), default="seq",
        help="criterion used by rmdup only (default: seq, i.e. identical sequences)",
    )
    p.add_argument(
        "--group", choices=("dir", "all"), default="dir",
        help="how the files are grouped for `seqkit common`: dir (default, the "
             "files of a same directory) or all (every file of the tree, useful "
             "when each class lives in its own subdirectory)",
    )
    p.add_argument(
        "--skip-extra", action="store_true",
        help="do not copy the non-FASTA files (they are copied as-is by default)",
    )
    p.add_argument(
        "--keep-common", action="store_true",
        help="do not remove the common sequences (only rmdup is applied)",
    )
    p.add_argument(
        "-j", "--threads", type=int, default=4, help="threads given to seqkit (default: 4)"
    )
    p.add_argument("--seqkit", default=None, help="path to the seqkit binary")
    p.add_argument(
        "-f", "--force", action="store_true", help="overwrite the output directory if it exists"
    )
    p.add_argument("-v", "--verbose", action="store_true", help="do not silence seqkit")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    src_root: Path = args.input_dir.expanduser().resolve()
    if not src_root.is_dir():
        sys.exit(f"❌ Not a directory: {src_root}")

    out_root: Path = (
        args.output_dir.expanduser().resolve()
        if args.output_dir
        else src_root.with_name(f"{src_root.name}_clean")
    )
    if out_root == src_root:
        sys.exit("❌ The output directory must be different from the input directory")
    if out_root.exists():
        if not args.force:
            sys.exit(f"❌ Output directory already exists: {out_root} (use --force)")
        shutil.rmtree(out_root)

    sk = SeqKit(resolve_seqkit(args.seqkit), args.threads, quiet=not args.verbose)

    fasta_files = sorted(p for p in src_root.rglob("*") if p.is_file() and is_fasta(p))
    if not fasta_files:
        sys.exit(f"❌ No FASTA file found under {src_root}")

    meta_dir = out_root / META_DIRNAME
    meta_dir.mkdir(parents=True, exist_ok=True)

    log(f"🧬 SeqKit: {sk.version()}")
    log(f"📂 Input : {src_root}")
    log(f"📂 Output: {out_root}")
    log(f"📄 {len(fasta_files)} FASTA file(s) found\n")

    report: dict = {
        "input_dir": str(src_root),
        "output_dir": str(out_root),
        "seqkit": sk.version(),
        "rmdup_by": args.rmdup_by,
        "common_by": args.by,
        "files": {},
        "groups": {},
    }

    # ---- (i) rmdup, file by file -------------------------------------------
    log("1️⃣  seqkit rmdup — removing duplicated sequences in each file")
    dedup_dir = Path(tempfile.mkdtemp(prefix="greac_clean_"))
    rel_paths = [f.relative_to(src_root) for f in fasta_files]

    def _dedup(rel: Path) -> tuple[Path, int, int]:
        src = src_root / rel
        tmp = dedup_dir / rel
        before = sk.count(src)
        step_rmdup(sk, src, tmp, args.rmdup_by)
        return rel, before, sk.count(tmp)

    try:
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
            for rel, before, after in pool.map(_dedup, rel_paths):
                log(f"   • {rel}: {before} → {after} sequences")
                report["files"][str(rel)] = {
                    "input": before, "after_rmdup": after, "removed_common": 0,
                }

        # ---- (ii)+(iii)+(iv) common sequences, group by group --------------
        groups: dict[Path, list[Path]] = defaultdict(list)
        for rel in rel_paths:
            groups[Path(".") if args.group == "all" else rel.parent].append(rel)

        for group_dir, rels in sorted(groups.items()):
            group_name = "<all>" if args.group == "all" else (
                str(group_dir) if str(group_dir) != "." else "<root>"
            )
            files = [dedup_dir / r for r in rels]
            keys: set[str] = set()

            if args.keep_common or len(files) < 2:
                if not args.keep_common:
                    log(f"\n2️⃣  {group_name}: single file, no common sequence to look for")
            else:
                log(f"\n2️⃣  seqkit common — group {group_name} ({len(files)} files, by {args.by})")
                pair_dir = Path(tempfile.mkdtemp(prefix="greac_pairs_", dir=dedup_dir))
                keys = step_common_keys(sk, files, args.by, pair_dir)
                shutil.rmtree(pair_dir, ignore_errors=True)
                log(f"   • {len(keys)} common sequence(s) found")

            safe = group_name.strip("<>").replace("/", "_").strip("_") or "root"
            report["groups"][group_name] = {"files": len(files), "common": len(keys)}

            # ids to drop, per file (a key matched by sequence may carry
            # different ids in different files)
            drop: dict[Path, set[str]] = {}
            if keys:
                for rel, f in zip(rels, files):
                    drop[rel] = {i for i, k in file_keys(sk, f, args.by) if k in keys}

            if any(drop.values()):
                ids_dir = meta_dir / "ids"
                per_file: list[tuple[Path, Path]] = []
                all_ids: set[str] = set()
                for rel, f in zip(rels, files):
                    ids = drop.get(rel, set())
                    all_ids |= ids
                    if ids:
                        per_file.append(
                            (f, write_ids(ids, ids_dir / f"{safe}__{str(rel).replace('/', '_')}.ids.txt"))
                        )
                write_ids(all_ids, meta_dir / f"common_{safe}.ids.txt")

                ref = meta_dir / f"common_{safe}.fasta"
                log("3️⃣  seqkit seq — extracting the common sequences to the reference file")
                step_extract_reference(sk, per_file, ref)
                log(f"   • reference: {ref.relative_to(out_root)}")
                report["groups"][group_name]["reference"] = str(ref.relative_to(out_root))

                log("4️⃣  seqkit grep — writing the filtered files")
                for rel, f in zip(rels, files):
                    dst = out_root / rel
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    ids = drop.get(rel, set())
                    entry = report["files"][str(rel)]
                    if ids:
                        step_filter(sk, f, dst, ids_dir / f"{safe}__{str(rel).replace('/', '_')}.ids.txt")
                        kept = sk.count(dst)
                    else:
                        shutil.copy2(f, dst)
                        kept = entry["after_rmdup"]
                    entry["removed_common"] = entry["after_rmdup"] - kept
                    entry["output"] = kept
                    warn = "  ⚠️  empty file!" if kept == 0 else ""
                    log(f"   • {rel}: {entry['after_rmdup']} → {kept} sequences{warn}")
            else:
                for rel, f in zip(rels, files):
                    dst = out_root / rel
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(f, dst)
                    report["files"][str(rel)]["output"] = report["files"][str(rel)]["after_rmdup"]
    finally:
        shutil.rmtree(dedup_dir, ignore_errors=True)

    # companion files (references, metadata, models, ...) are copied untouched
    extra = 0
    if not args.skip_extra:
        for src in src_root.rglob("*"):
            if not src.is_file() or is_fasta(src):
                continue
            dst = out_root / src.relative_to(src_root)
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            extra += 1
        if extra:
            log(f"\n\U0001f4ce {extra} non-FASTA file(s) copied as-is")
    report["extra_files_copied"] = extra

    (meta_dir / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    total_in = sum(v["input"] for v in report["files"].values())
    total_out = sum(v.get("output", 0) for v in report["files"].values())
    log("\n✅ Done")
    log(f"   sequences: {total_in} → {total_out}")
    log(f"   clean dataset: {out_root}")
    log(f"   report: {(meta_dir / 'report.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
