# rst_caller.py
# Michael J. Foster
# github.com/mjfos2r
# 2025-May-27
# Version 0.1.0

# Requires Seqkit: https://bioinf.shenwei.me/seqkit/usage/#amplicon
# Refer to the associated markdown files in the vault.

from __future__ import annotations

import argparse
import io
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import Restriction, Seq, SeqIO, SeqRecord
from Bio import __version__ as BioPythonVersion

from . import __about__ as about

__version__ = about.__version__

# patterns, primers, and type definitions are from Liveris et al. (1995)
RFLP_PATTERNS: dict[str, list[int]] = {
    "H1": [1078, 372, 310],
    "H2": [1078, 310, 241, 131],
    "M1": [258, 149, 136, 128, 102],
    "M2": [364, 258, 136, 102],
}

RST_TYPES: dict[int : list[str]] = {
    1: ["H1", "M1"],
    2: ["H2", "M2"],
    3: ["H2", "M1"],
}


class _FullVersion(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print(f"{parser.prog} {__version__}")
        print(f"Biopython {BioPythonVersion}")
        print(_ensure_deps())
        parser.exit(0)


def _ensure_deps() -> str:
    prog = shutil.which("seqkit")
    if prog is None:
        sys.exit(
            "ERROR: Cannot locate SeqKit binary. Please download the latest release for your OS from github!\nhttps://github.com/shenwei356/seqkit/releases"
        )
    else:
        version = subprocess.run(
            ["seqkit", "version"],
            text=True,
            capture_output=True,
        )
        return version.stdout


def _existing_file(path_str: str) -> Path:
    p = Path(path_str)
    if not p.is_file():
        raise argparse.ArgumentTypeError(f"Input FASTA not found: {p}")
    return p


def fuzzy_match(observed, expected, tolerance=10) -> bool:
    unmatched = expected[:]
    for obs in observed:
        for exp in unmatched:
            if abs(obs - exp) <= tolerance:
                unmatched.remove(exp)
                break
    return len(unmatched) == 0


def get_best_pattern(obs_frags, enzyme, tolerance=10) -> list | None:
    candidates = None
    for pat_name, pat_frags in RFLP_PATTERNS.items():
        if pat_name.startswith(enzyme[0]) and fuzzy_match(
            obs_frags, sorted(pat_frags), tolerance
        ):
            candidates: str = pat_name
    return candidates


def amplify_and_cut(input_fa_file, output: Path, quiet=False, **kwargs):
    # primers from Liveris et al. (1995)
    fwd = "GGTATGTTTAGTGAGGG"
    rev = "CAGGCTCTACACTTCTG"

    file_id = input_fa_file.stem

    with open(input_fa_file, "r") as f:
        proc_out = subprocess.run(
            [
                "seqkit",
                "amplicon",
                "--bed",
                "-F",
                fwd,
                "-R",
                rev,
                "-m",
                "2",
                "-M",
                "-I",
            ],
            stdin=f,
            text=True,
            capture_output=True,
        )
        output_file = output / f"{file_id}_AMPLICON.fna"
        with open(output_file, "w") as outhandle:
            # convert to BED6+1 to add coordinates to header?
            # let's try it.
            # (0, 'chrom') (1, 'start') (2, 'end') (3, 'name') (4, 'score') (5, 'strand') (6, 'sequence') (7, 'mismatches') (8, 'mismatches_5') (9, 'mismatches_3')
            # keys: list[str] = ["chrom", "start", "end", "name", "score", "strand", "sequence", "mismatches", "mismatches_5", "mismatches_3"]
            lines = proc_out.stdout.split("\n")
            if len(lines) == 0:
                if not quiet:
                    print("ERROR: No amplicon returned!")
                return (
                    'NoAMP',
                    [],
                    [],
                    'NaN',
                    'NaN',
                    amplicons,
                    amplicon,
                )
            amplicons = []
            for line in lines:
                if not line.strip():
                    continue
                fields: list[str] = line.split("\t")
                rec_id: str = fields[0]
                rec_name: str = "rRNA_ITS_amplicon"
                rec_description: str = f"[{fields[1]}:{fields[2]}] {file_id}_{rec_name} mismatches={fields[7]}({fields[8]}+{fields[9]})"
                amplicons.append(
                    SeqRecord.SeqRecord(
                        seq=Seq.Seq(fields[6]),
                        id=rec_id,
                        name=rec_name,
                        description=rec_description,
                    )
                )
                amplicon: SeqRecord = amplicons[0]
                amplicon_out: int = SeqIO.write(amplicons, outhandle, "fasta")
                if not quiet:
                    print(f"{amplicon_out} amplicons written.")
        # run the digestion but toss out the MseI fragments below 100. (since we're fuzzy matching, drop this to 90bp)
        HinfI_fragments: list[int] = sorted(
            {
                len(frag)
                for frag in Restriction.HinfI.catalyze(amplicon.seq, linear=True)
            }
        )
        MseI_fragments: list[int] = sorted(
            {
                len(frag)
                for frag in Restriction.MseI.catalyze(amplicon.seq, linear=True)
                if len(frag) >= 90
            }
        )

        # this is not optimal but this is the only way I can get it to return types that have been experimentally determined.
        # bump the tolerance for M_match to 20.
        H_match: list | None = get_best_pattern(HinfI_fragments, "HinfI", 45)
        M_match: list | None = get_best_pattern(MseI_fragments, "MseI", 20)
        if not quiet:
            print(amplicon.description.split(" ")[-1], len(amplicon.seq))
            print(f"HinfI fragments: {HinfI_fragments}")
            print(f"MseI fragments:  {MseI_fragments}")
            print(f"Matched HinfI Pattern: {H_match}")
            print(f"Matched MseI Pattern:  {M_match}")

        if not H_match or not M_match:
            # I want to return 0 if one of the patterns doesn't match
            # but isn't RST3 anything not 1 or 2? 0 for now to make debugging easier.
            if not quiet:
                print("Error: Novel fragment pattern detected!")
                # should make it always output all fragments
                # TODO: see above
            return (
                0,
                HinfI_fragments,
                MseI_fragments,
                H_match,
                M_match,
                amplicons,
                amplicon,
            )

        for rst, (hpat, mpat) in RST_TYPES.items():
            if hpat in H_match and mpat in M_match:
                called_RST: int = rst
        if not quiet:
            print(f"Matched RST: {called_RST}")

        return (
            called_RST,
            HinfI_fragments,
            MseI_fragments,
            H_match,
            M_match,
            amplicons,
            amplicon,
        )


def main():
    parser = argparse.ArgumentParser(
        prog="rst_caller",
        description="RST type classifier for a single FASTA file containing B.burgdorferi contigs. This does insilico PCR and RFLP to determine the type",
    )
    parser.add_argument(
        "-v",
        "--version",
        nargs=0,
        action=_FullVersion,
        help="Show program and BLAST versions, then exit.",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=_existing_file,
        required=True,
        help="Input FASTA file path",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Directory where results will be output. It will be created if nonexistant.",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        required=False,
        default=False,
        help="Suppress verbose output to stdout.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    results = amplify_and_cut(args.input, args.output, args.quiet)

    with open(args.output / f"{args.input.stem}_RST_TYPE.txt", "w") as outf:
        outf.write(f"{results[0]}\n")

    with open(args.output / f"{args.input.stem}_FRAGMENT_LENGTHS.txt", "w") as outf:
        # results = (0, 'called_RST') (1, 'HinfI_fragments') (2, 'MseI_fragments') (3, 'H_match') (4, 'M_match') (5, 'amplicons') (6, 'amplicon')
        header = [
            "enzyme",
            "rflp_pattern_name",
            "rflp_pattern_fragments",
            "observed_fragments",
        ]
        H_match_frags = RFLP_PATTERNS.get(results[3], "None")
        M_match_frags = RFLP_PATTERNS.get(results[4], "None")  # results[4]
        outf.write(f"{'\t'.join(header)}\n")
        outf.write(f"HinfI\t{results[3]}\t{H_match_frags}\t{results[1]}\n")
        outf.write(f"MseI\t{results[4]}\t{M_match_frags}\t{results[2]}\n")
    sys.exit(0)


if __name__ == "__main__":
    main()
