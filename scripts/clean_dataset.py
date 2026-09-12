#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Dataset cleaning pipeline for GREAC FASTA datasets, built on SeqKit2.

Given a directory holding FASTA files (searched recursively), this script
writes a sibling directory `<input_dir>_clean` with the very same tree and the
very same file names, so any downstream GREAC/GRAMEP workflow can be pointed at
the clean copy without any other change.

Cleaning rules, all performed by SeqKit2 commands:

  (i)   seqkit rmdup   -> remove duplicated sequences inside each file;
  (ii)  seqkit common  -> find the sequences shared between the files of a
                          group, compared pairwise;
  (iii) seqkit seq     -> extract those shared sequences to a reference file;
  (iv)  seqkit grep -v -> write the filtered files, without the shared
                          sequences, to avoid over-fitting.

A group is, by default, the set of FASTA files of a same directory; with
`--group all` the whole tree is a single group, which is what the GREAC layout
needs when every class lives in its own subdirectory. Non-FASTA files are
copied untouched. The reference files and the id lists are kept in the hidden
folder `.clean_meta/` of the output directory, so they never show up in a
`*.fasta` glob of the clean dataset.

Requires SeqKit2 in the PATH (https://bioinf.shenwei.me/seqkit/).
Compatible with Python 3.6+.

Usage:
    python clean_dataset.py <fasta_dir> [options]

Examples:
    python clean_dataset.py ~/datasets/virus_dataset
    python clean_dataset.py ./virus_dataset --group all --by seq -j 8
"""

import argparse
import itertools
import os
import shutil
import subprocess
import sys
import tempfile

FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".faa", ".ffn", ".frn")
COMPRESSED_EXTENSIONS = (".gz", ".xz", ".bz2", ".zst")
META_DIRNAME = ".clean_meta"

# comparison criteria of `seqkit common`
BY_FLAGS = {"id": [], "name": ["-n"], "seq": ["-s"]}

SEQKIT = "seqkit"
THREADS = "4"


def log(msg):
    print(msg)
    sys.stdout.flush()


def die(msg):
    sys.stderr.write(msg + "\n")
    sys.exit(1)


def is_fasta(path):
    name = path.lower()
    for comp in COMPRESSED_EXTENSIONS:
        if name.endswith(comp):
            name = name[: -len(comp)]
            break
    return name.endswith(FASTA_EXTENSIONS)


def find_fasta(root):
    """Relative paths of every FASTA file under root, recursively."""
    found = []
    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            full = os.path.join(dirpath, fn)
            if is_fasta(fn):
                found.append(os.path.relpath(full, root))
    return sorted(found)


def seqkit(args, stdout=None):
    """Run a seqkit command; return its stdout, or write it to a file."""
    cmd = [SEQKIT] + args + ["--threads", THREADS, "--quiet"]
    if stdout is not None:
        with open(stdout, "w") as fh:
            proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.PIPE,
                                    universal_newlines=True)
            _, err = proc.communicate()
        out = ""
    else:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        out, err = proc.communicate()
    if proc.returncode != 0:
        die("[ERRO] seqkit falhou (%d): %s\n%s" % (proc.returncode, " ".join(cmd), err.strip()))
    return out


def count(fasta):
    """Number of sequences in a file, via `seqkit stats`."""
    lines = seqkit(["stats", "-T", fasta]).splitlines()
    try:
        return int(lines[1].split("\t")[3])
    except (IndexError, ValueError):
        return 0


def keys_of(fasta, by):
    """[(id, comparison key)] of every record of a file, listed by `seqkit seq`.

    The key is what `seqkit common` compares (id, full name or sequence), while
    the id is what `seqkit grep` filters on, so a match found by sequence is
    still removed from the right records of every file, even when their ids
    differ.
    """
    if by == "seq":
        # sequences are only materialised when the comparison needs them
        pairs = []
        for line in seqkit(["fx2tab", fasta]).splitlines():
            if line.strip():
                name, _, seq = line.partition("\t")
                pairs.append((name.split()[0], seq.strip().upper()))
        return pairs

    flags = ["-n"] if by == "name" else ["-n", "-i"]
    return [(l.split()[0], l.strip()) for l in seqkit(["seq"] + flags + [fasta]).splitlines()
            if l.strip()]


def common_keys(files, by, workdir):
    """(ii) keys of the sequences shared between the files of a group.

    `seqkit common` returns what is present in *all* the given files, so the
    files are compared pairwise and the results joined: anything seen in more
    than one file of the group counts as common.
    """
    keys = set()
    for idx, (a, b) in enumerate(itertools.combinations(files, 2)):
        pair = os.path.join(workdir, "pair_%d.fasta" % idx)
        seqkit(["common"] + BY_FLAGS[by] + [a, b], stdout=pair)
        if os.path.getsize(pair) > 0:
            # (iii) the shared records are read back with `seqkit seq`
            keys.update(k for _, k in keys_of(pair, by))
        os.remove(pair)
    return keys


def write_lines(items, path):
    mkdir(os.path.dirname(path))
    with open(path, "w") as fh:
        fh.write("\n".join(sorted(items)))
        fh.write("\n")
    return path


def mkdir(path):
    if path and not os.path.isdir(path):
        os.makedirs(path)


def clean_group(name, rels, sources, out_root, meta_dir, by, counts):
    """Steps (ii) to (iv) for one group, writing the clean files."""
    safe = name.strip("<>").replace(os.sep, "_").strip("_") or "root"
    keys = set()

    if len(sources) < 2:
        log("\n[2] %s: um unico ficheiro, nada em comum a procurar" % name)
    else:
        log("\n[2] seqkit common - grupo %s (%d ficheiros, por %s)" % (name, len(sources), by))
        pair_dir = tempfile.mkdtemp(prefix="greac_pairs_")
        try:
            keys = common_keys(sources, by, pair_dir)
        finally:
            shutil.rmtree(pair_dir, ignore_errors=True)
        log("    %d sequencia(s) em comum" % len(keys))

    # ids to drop, per file: the same sequence may carry different ids
    drop = {}
    if keys:
        for rel, src in zip(rels, sources):
            drop[rel] = set(i for i, k in keys_of(src, by) if k in keys)

    if not any(drop.values()):
        for rel, src in zip(rels, sources):
            copy_to(src, os.path.join(out_root, rel))
            counts[rel]["output"] = counts[rel]["rmdup"]
        return

    ids_dir = os.path.join(meta_dir, "ids")
    ids_files = {}
    for rel in rels:
        if drop.get(rel):
            flat = rel.replace(os.sep, "_")
            ids_files[rel] = write_lines(drop[rel], os.path.join(ids_dir, "%s__%s.ids.txt" % (safe, flat)))
    write_lines(set().union(*drop.values()), os.path.join(meta_dir, "common_%s.ids.txt" % safe))

    # (iii) reference file with one copy of each shared sequence
    log("[3] seqkit seq - a extrair as sequencias comuns para o ficheiro de referencia")
    ref = os.path.join(meta_dir, "common_%s.fasta" % safe)
    tmp_ref = ref + ".tmp"
    with open(tmp_ref, "w") as fh:
        for rel, src in zip(rels, sources):
            if rel in ids_files:
                fh.write(seqkit(["grep", "-f", ids_files[rel], src]))
    seqkit(["rmdup", "-s", tmp_ref], stdout=ref)
    os.remove(tmp_ref)
    log("    referencia: %s" % os.path.join(META_DIRNAME, os.path.basename(ref)))

    # (iv) filtered files
    log("[4] seqkit grep - a escrever os ficheiros filtrados")
    for rel, src in zip(rels, sources):
        dst = os.path.join(out_root, rel)
        mkdir(os.path.dirname(dst))
        if rel in ids_files:
            seqkit(["grep", "-v", "-f", ids_files[rel], src], stdout=dst)
            kept = count(dst)
        else:
            copy_to(src, dst)
            kept = counts[rel]["rmdup"]
        counts[rel]["output"] = kept
        warn = "   [AVISO] ficheiro vazio!" if kept == 0 else ""
        log("    %s: %d -> %d sequencias%s" % (rel, counts[rel]["rmdup"], kept, warn))


def copy_to(src, dst):
    mkdir(os.path.dirname(dst))
    shutil.copy2(src, dst)


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Limpa um dataset FASTA com SeqKit2, espelhando a arvore de diretorios.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("input_dir", help="diretorio com os ficheiros FASTA")
    p.add_argument("-o", "--output-dir", default=None,
                   help="diretorio de saida (por omissao: <input_dir>_clean)")
    p.add_argument("--by", choices=("id", "name", "seq"), default="id",
                   help="criterio do seqkit common: id (omissao), name ou seq")
    p.add_argument("--group", choices=("dir", "all"), default="dir",
                   help="agrupamento para o seqkit common: dir (omissao, ficheiros da mesma "
                        "pasta) ou all (toda a arvore, util quando cada classe tem a sua subpasta)")
    p.add_argument("-j", "--threads", default="4", help="threads passadas ao seqkit (omissao: 4)")
    p.add_argument("--seqkit", default=None, help="caminho do binario seqkit")
    p.add_argument("-f", "--force", action="store_true",
                   help="sobrescreve o diretorio de saida se ja existir")
    return p.parse_args(argv)


def main(argv=None):
    global SEQKIT, THREADS
    args = parse_args(argv)
    THREADS = str(args.threads)

    for cand in ([args.seqkit] if args.seqkit else ["seqkit", "seqkit2"]):
        if shutil.which(cand):
            SEQKIT = cand
            break
    else:
        die("[ERRO] SeqKit2 nao encontrado. Instale-o (https://bioinf.shenwei.me/seqkit/) "
            "ou indique o binario com --seqkit /caminho/para/seqkit")

    src_root = os.path.abspath(os.path.expanduser(args.input_dir))
    if not os.path.isdir(src_root):
        die("[ERRO] nao e um diretorio: %s" % src_root)

    if args.output_dir:
        out_root = os.path.abspath(os.path.expanduser(args.output_dir))
    else:
        out_root = src_root.rstrip(os.sep) + "_clean"
    if out_root == src_root:
        die("[ERRO] o diretorio de saida tem de ser diferente do de entrada")
    if os.path.exists(out_root):
        if not args.force:
            die("[ERRO] o diretorio de saida ja existe: %s (use --force)" % out_root)
        shutil.rmtree(out_root)

    rels = find_fasta(src_root)
    if not rels:
        die("[ERRO] nenhum ficheiro FASTA encontrado em %s" % src_root)

    meta_dir = os.path.join(out_root, META_DIRNAME)
    mkdir(meta_dir)

    log("SeqKit : %s" % seqkit(["version"]).strip())
    log("Entrada: %s" % src_root)
    log("Saida  : %s" % out_root)
    log("%d ficheiro(s) FASTA encontrado(s)\n" % len(rels))

    # ---- (i) rmdup, file by file ----------------------------------------
    log("[1] seqkit rmdup - a remover sequencias duplicadas em cada ficheiro")
    tmp_root = tempfile.mkdtemp(prefix="greac_clean_")
    counts = {}
    try:
        for rel in rels:
            src = os.path.join(src_root, rel)
            tmp = os.path.join(tmp_root, rel)
            mkdir(os.path.dirname(tmp))
            before = count(src)
            seqkit(["rmdup", "-s", src], stdout=tmp)
            after = count(tmp)
            counts[rel] = {"input": before, "rmdup": after, "output": after}
            log("    %s: %d -> %d sequencias" % (rel, before, after))

        # ---- (ii) + (iii) + (iv), group by group -------------------------
        groups = {}
        for rel in rels:
            key = "<all>" if args.group == "all" else (os.path.dirname(rel) or "<root>")
            groups.setdefault(key, []).append(rel)

        for name in sorted(groups):
            group_rels = groups[name]
            sources = [os.path.join(tmp_root, r) for r in group_rels]
            clean_group(name, group_rels, sources, out_root, meta_dir, args.by, counts)
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)

    # ---- companion files (references, metadata, models, ...) -------------
    extra = 0
    for dirpath, _, filenames in os.walk(src_root):
        for fn in filenames:
            if is_fasta(fn):
                continue
            src = os.path.join(dirpath, fn)
            copy_to(src, os.path.join(out_root, os.path.relpath(src, src_root)))
            extra += 1
    if extra:
        log("\n%d ficheiro(s) nao-FASTA copiado(s) sem alteracao" % extra)

    total_in = sum(c["input"] for c in counts.values())
    total_out = sum(c["output"] for c in counts.values())
    log("\nConcluido")
    log("    sequencias: %d -> %d" % (total_in, total_out))
    log("    dataset limpo: %s" % out_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
