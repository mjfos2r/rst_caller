# rst_caller.py
# Michael J. Foster
# github.com/mjfos2r
# 2025-May-27
# Version 0.1.0

# Requires Seqkit: https://bioinf.shenwei.me/seqkit/usage/#amplicon
import subprocess
import io
import sys
from Bio import Restriction
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
import argparse

def _existing_file(path_str: str) -> Path:
    p = Path(path_str)
    if not p.is_file():
        raise argparse.ArgumentTypeError(f"Input FASTA not found: {p}")
    return p

def fuzzy_match(observed, expected, tolerance=10):
    unmatched = expected[:]
    for obs in observed:
        for exp in unmatched:
            if abs(obs-exp) <= tolerance:
                unmatched.remove(exp)
                break
    return len(unmatched) == 0

def get_best_pattern(obs_frags, enzyme, tolerance=10):
    candidates = []
    for pat_name, pat_frags in RFLP_PATTERNS.items():
        if pat_name.startswith(enzyme[0]) and fuzzy_match(obs_frags, sorted(pat_frags), tolerance):
            candidates.append(pat_name)
    return candidates

def amplify_and_cut(input_fa_file, quiet=False):
    # primers from Liveris et al. (1995)
    fwd = "GGTATGTTTAGTGAGGG"
    rev = "CAGGCTCTACACTTCTG"

    with open(input_fa_file, 'r') as f:
        proc_out = subprocess.run(
            ["seqkit", "amplicon",
             "-F", fwd, "-R", rev,
             "-m", "2", "-M", "-I"],
            stdin=f,
            text=True,
            capture_output=True,
        )

        amplicon = [rec for rec in SeqIO.parse(io.StringIO(proc_out.stdout), 'fasta')][0]
        # run the digestion but toss out the MseI fragments below 100.
        HinfI_fragments = sorted({len(frag) for frag in Restriction.HinfI.catalyze(amplicon.seq, linear=True)})
        MseI_fragments = sorted({len(frag) for frag in Restriction.MseI.catalyze(amplicon.seq, linear=True) if len(frag) >= 100})

        # this is not optimal but this is the only way I can get it to return types that have been experimentally determined.
        # fuzzy matching +/- 45bp doesn't feel right but it's working.
        # given that most HinfI fragments are 1037 instead of 1078 let's give it some wiggle room to account for indels.
        # and MseI is usually pretty close. let's use a semi-arbitrary cutoff of 17bp 
        # since Jones et al. 2006 use 18bp as their minimum length to distinguish bands. (for ospC typing via RFLP)
        # This needs review and work.
        H_match = get_best_pattern(HinfI_fragments, 'HinfI', 45)
        M_match = get_best_pattern(MseI_fragments, 'MseI', 17)
        if not quiet:
            print(amplicon.description.split(' ')[-1], len(amplicon.seq))
            print(f"HinfI fragments: {HinfI_fragments}")
            print(f" MseI fragments: {MseI_fragments}")
            print(f"Matched HinfI Pattern: {H_match}")
            print(f" Matched MseI Pattern: {M_match}")

        if not H_match or not M_match:
            # I want to return 0 if one of the patterns doesn't match
            # but isn't RST3 anything not 1 or 2? 0 for now to make debugging easier.
            if not quiet:
               print("Error: Novel fragment pattern detected!")
            return 0,  HinfI_fragments, MseI_fragments

        for rst, (hpat, mpat) in RST_TYPES.items():
            if hpat in H_match and mpat in M_match:
                called_RST = rst

        if not quiet:
            print(f"Matched RST: {called_RST}")
                
        return called_RST, HinfI_fragments, MseI_fragments, hpat, mpat

# patterns, primers, and type definitions are from Liveris et al. (1995)
RFLP_PATTERNS = {
    'H1':[1078,372,310],
    'H2':[1078,310,241,131],
    'M1':[258,149,136,128,102],
    'M2':[364,258,136,102],
}

RST_TYPES = {
    1: ['H1', 'M1'],
    2: ['H2', 'M2'],
    3: ['H2', 'M1'],
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="rst_caller",
        description="RST type classifier for a single FASTA file containing B.burgdorferi contigs. This does insilico PCR and RFLP to determine the type",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=_existing_file,
        required=True,
        help="Input FASTA file path"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output directory"
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        required=False,
        default=False,
        help="Suppress verbose output to stdout."
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    results = amplify_and_cut(args.input, args.quiet)

    with open(args.output / "RST_TYPE.txt", 'w') as outf:
        outf.write(f"{results[0]}\n")
    with open(args.output / "FRAGMENT_LENGTHS.txt", 'w') as outf:
        header = ["enzyme", "pattern_name", "pattern_fragments", "observed_fragments"]
        outf.write(f"HinfI\t{results[3]}\t{RFLP_PATTERNS[results[3]]}\t{results[1]}\n")
        outf.write(f"MseI\t{results[4]}\t{RFLP_PATTERNS[results[4]]}\t{results[2]}\n")
    sys.exit(0)

