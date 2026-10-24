"""
Helpers used by cleaning.smk to locate and normalise the input reads of a sample.

A sample is a directory (``<Dataset>/<sample>``) that contains FASTQ files
(``.fastq``, ``.fq``, ``.fastq.gz`` or ``.fq.gz``).

Illumina data can be delivered in two layouts:

* ``paired``      one file per mate, e.g. ``S1_L001_R1_001.fastq.gz`` and
                  ``S1_L001_R2_001.fastq.gz`` (several lanes/chunks are allowed);
* ``interleaved`` one file in which forward and reverse reads alternate
                  (R1, R2, R1, R2, ...);
* ``single``      single-end reads: every file holds independent reads without mates
                  (several files/chunks are allowed and are concatenated).

Nanopore and PacBio reads are single files per chunk and only need to be
concatenated (see ``concat_fastq``).

The module only depends on the Python standard library so it can be tested
without Snakemake.
"""
import gzip
import os
import re
import shutil

FASTQ_NAME = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)

# <stem><sep><R?><1|2><_001?>.<fastq|fq>[.gz]
#   S1_L001_R1_001.fastq.gz | sample_R1.fq.gz | sample_1.fastq | sample.R2.fastq.gz
MATE_NAME = re.compile(
    r"^(?P<stem>.+?)(?P<sep>[._-])(?P<r>R?)(?P<mate>[12])(?P<tail>(?:_001)?)"
    r"(?P<ext>\.(?:fastq|fq)(?:\.gz)?)$",
    re.IGNORECASE,
)

LAYOUTS = ("auto", "paired", "interleaved", "single")


# ---------------------------------------------------------------------------
# Low level FASTQ access
# ---------------------------------------------------------------------------
def _open_binary(path):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def _read_record(handle, path):
    """Return the next FASTQ record as four byte lines, or None at end of file."""
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        return None
    if not lines[3]:
        raise ValueError(f"{path}: truncated FASTQ record (incomplete last entry)")
    if lines[0][:1] != b"@" or lines[2][:1] != b"+":
        raise ValueError(
            f"{path}: malformed FASTQ record near '{lines[0][:60].decode('ascii', 'replace').strip()}'"
        )
    if not lines[3].endswith(b"\n"):
        lines[3] += b"\n"
    return lines


def _record_id(header):
    """Read identifier shared by both mates of a pair (without /1, /2 or the CASAVA mate field)."""
    text = header.decode("ascii", "replace").strip()
    if text.startswith("@"):
        text = text[1:]
    token = text.split()[0] if text.split() else ""
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def _first_record(path):
    with _open_binary(path) as handle:
        return _read_record(handle, path)


# ---------------------------------------------------------------------------
# Locating files
# ---------------------------------------------------------------------------
def list_fastq(directory):
    """Sorted FASTQ files located directly inside ``directory``."""
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"reads directory not found: {directory}")
    files = sorted(
        os.path.join(directory, name)
        for name in os.listdir(directory)
        if FASTQ_NAME.search(name) and os.path.isfile(os.path.join(directory, name))
    )
    if not files:
        raise FileNotFoundError(
            f"no FASTQ files (.fastq, .fq, .fastq.gz, .fq.gz) found in {directory}"
        )
    return files


def is_interleaved(path, n_pairs=3):
    """True when the first records of ``path`` alternate between the two mates of a pair."""
    checked = 0
    with _open_binary(path) as handle:
        for _ in range(n_pairs):
            first = _read_record(handle, path)
            if first is None:
                break
            second = _read_record(handle, path)
            if second is None:
                return False
            if _record_id(first[0]) != _record_id(second[0]):
                return False
            checked += 1
    return checked > 0


def pair_by_name(files):
    """
    Pair forward/reverse files by their names.

    Returns two lists (forward, reverse) in matching order. Raises ValueError if a file
    has no recognisable mate marker or if a mate is missing.
    """
    groups = {}
    for path in files:
        match = MATE_NAME.match(os.path.basename(path))
        if match is None:
            raise ValueError(
                f"{os.path.basename(path)}: cannot identify the mate; paired files must contain "
                "R1/R2 or 1/2 (e.g. sample_R1.fastq.gz and sample_R2.fastq.gz)"
            )
        key = (match["stem"], match["sep"], match["r"].upper(), match["tail"])
        mate = int(match["mate"])
        if mate in groups.setdefault(key, {}):
            raise ValueError(
                f"{os.path.basename(path)}: more than one file for mate {mate} of the same pair"
            )
        groups[key][mate] = path

    forward, reverse = [], []
    for key in sorted(groups):
        mates = groups[key]
        if 1 not in mates or 2 not in mates:
            present = mates[1] if 1 in mates else mates[2]
            raise ValueError(f"{os.path.basename(present)}: its mate file was not found")
        forward.append(mates[1])
        reverse.append(mates[2])
    return forward, reverse


def _check_pairs_in_sync(forward, reverse):
    """Compare the first read of every forward/reverse file pair."""
    for f_path, r_path in zip(forward, reverse):
        f_rec, r_rec = _first_record(f_path), _first_record(r_path)
        if f_rec is None or r_rec is None:
            raise ValueError(f"{f_path} / {r_path}: empty file in a pair")
        if _record_id(f_rec[0]) != _record_id(r_rec[0]):
            raise ValueError(
                f"{os.path.basename(f_path)} and {os.path.basename(r_path)} do not start with the "
                "same read name; they are not mates of the same pair (if the files are interleaved "
                "set library_layout: interleaved)"
            )


def resolve_illumina_sample(directory, layout="auto"):
    """
    Describe the Illumina input of one sample.

    ``layout`` is ``auto`` (detect), ``paired``, ``interleaved`` or ``single``.
    Returns ``{"layout": "paired", "r1": [...], "r2": [...]}``,
    ``{"layout": "interleaved", "files": [...]}`` or ``{"layout": "single", "files": [...]}``.

    With ``auto`` the layout is decided as follows:
      1. every file carries an R1/R2 (or 1/2) marker           -> paired
      2. no file carries a marker and all alternate mates      -> interleaved
      3. no file carries a marker and none alternates mates    -> single
      4. anything else (mixed markers, mixed content)          -> error
    A single-end file whose name ends in _1 (e.g. SRR000001_1.fastq) is read as a mate by rule 1;
    in that case the error asks for ``library_layout: single``.
    """
    if layout not in LAYOUTS:
        raise ValueError(f"library_layout must be one of {', '.join(LAYOUTS)} (got '{layout}')")

    files = list_fastq(directory)

    if layout == "single":
        # explicit choice: the files are used as they are, without looking for mates
        return {"layout": "single", "files": files}

    with_mate = [f for f in files if MATE_NAME.match(os.path.basename(f))]
    without_mate = [f for f in files if f not in with_mate]

    if layout == "paired" or (layout == "auto" and with_mate and not without_mate):
        try:
            forward, reverse = pair_by_name(files)
            _check_pairs_in_sync(forward, reverse)
        except ValueError as error:
            if layout == "auto":
                raise ValueError(
                    f"{error} (if these files are single-end reads set library_layout: single)"
                ) from None
            raise
        return {"layout": "paired", "r1": forward, "r2": reverse}

    if layout == "auto" and with_mate and without_mate:
        raise ValueError(
            f"{directory}: files with and without R1/R2 markers are mixed; cannot decide the layout "
            "(separate the files, or set library_layout to paired, interleaved or single)"
        )

    # no mate markers in the names: interleaved (explicit, or detected) or single-end (detected)
    interleaved_flags = {f: is_interleaved(f) for f in files}
    if layout == "interleaved":
        not_interleaved = [os.path.basename(f) for f, flag in interleaved_flags.items() if not flag]
        if not_interleaved:
            raise ValueError(
                f"{directory}: not interleaved paired reads: {', '.join(not_interleaved)}. "
                "The files do not alternate forward/reverse reads "
                "(for single-end reads set library_layout: single)"
            )
        return {"layout": "interleaved", "files": files}

    # auto
    if all(interleaved_flags.values()):
        return {"layout": "interleaved", "files": files}
    if not any(interleaved_flags.values()):
        return {"layout": "single", "files": files}
    raise ValueError(
        f"{directory}: some files alternate forward/reverse reads and others do not "
        f"(interleaved: {', '.join(os.path.basename(f) for f, flag in interleaved_flags.items() if flag)}); "
        "separate them or set library_layout explicitly"
    )


# ---------------------------------------------------------------------------
# Writing normalised files
# ---------------------------------------------------------------------------
def concat_fastq(files, output, compresslevel=6):
    """Concatenate FASTQ files (plain or gzip) into one gzip file."""
    with open(output, "wb") as destination:
        for path in files:
            with open(path, "rb") as source:
                if path.lower().endswith(".gz"):
                    shutil.copyfileobj(source, destination, 1 << 20)  # gzip members can be chained
                else:
                    with gzip.GzipFile(
                        fileobj=destination, mode="wb", compresslevel=compresslevel
                    ) as compressed:
                        shutil.copyfileobj(source, compressed, 1 << 20)


def deinterleave(files, output_r1, output_r2, compresslevel=4):
    """
    Split interleaved FASTQ files into forward and reverse gzip files.

    Checks that every pair carries the same read name and returns the number of pairs.
    """
    pairs = 0
    with gzip.open(output_r1, "wb", compresslevel) as out1, gzip.open(
        output_r2, "wb", compresslevel
    ) as out2:
        for path in files:
            with _open_binary(path) as handle:
                while True:
                    first = _read_record(handle, path)
                    if first is None:
                        break
                    second = _read_record(handle, path)
                    if second is None:
                        raise ValueError(
                            f"{path}: odd number of reads; the last read has no mate "
                            f"({_record_id(first[0])})"
                        )
                    if _record_id(first[0]) != _record_id(second[0]):
                        raise ValueError(
                            f"{path}: consecutive reads are not mates of the same pair: "
                            f"'{first[0].decode('ascii', 'replace').strip()}' / "
                            f"'{second[0].decode('ascii', 'replace').strip()}'"
                        )
                    out1.writelines(first)
                    out2.writelines(second)
                    pairs += 1
    return pairs
