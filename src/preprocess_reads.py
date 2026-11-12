import csv
import gzip
import os
import re
import subprocess
import sys


sample_file = "results/sample_summary.tsv"
out_folder = "results/preprocessed"
tmp_folder = "results/preprocessed/tmp"
summary_file = "results/preprocessing_summary.tsv"

r1_adapter = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT"
r2_adapter = "TACCCTTCGCTTCACACACAAGGGGAAAGAGTGTAGATCTCGGTGGTCGC"
min_length = 14
cutadapt_exe = os.path.join(os.environ["USERPROFILE"],
                            "miniconda3", "envs", "bme415",
                            "Scripts", "cutadapt.exe")
umi_pattern = re.compile(r"([A-Z]{3})rand=([ACGT]{7})\b")


def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def load_samples(only_samples=None):
    if not os.path.exists(sample_file):
        print("Sample summary not found. Run Module 1 first.")
        sys.exit(1)
    rows = []
    with open(sample_file, newline="") as f:
        rows = [r for r in csv.DictReader(f, delimiter="\t")
                if r.get("sample_id")]
    if only_samples:
        wanted = [s.strip() for s in only_samples.split(",")]
        rows = [r for r in rows if r["sample_id"] in wanted]
        if len(rows) != len(wanted):
            print("ERROR: --samples contained unknown IDs")
            sys.exit(1)
    for r in rows:
        if r["status"] != "OK":
            print("ERROR: sample not OK:", r["sample_id"])
            sys.exit(1)
    return rows


def raw_paths(row):
    ""
    layout = row["layout"]
    if layout == "interleaved":
        for ext in (".fastq", ".fastq.gz"):
            p = os.path.join("raw_data", row["sample_id"] + ext)
            if os.path.exists(p):
                return p, p, True
        return None, None, True
    files = []
    for mate in ("_1", "_2"):
        found = None
        for ext in (".fastq", ".fastq.gz"):
            p = os.path.join("raw_data", row["sample_id"] + mate + ext)
            if os.path.exists(p):
                found = p
                break
        if found is None:
            return None, None, False
        files.append(found)
    return files[0], files[1], False


def umi_of(header):
    ""
    m = umi_pattern.search(header)
    if m:
        return m.group(2)
    return None


def rename_and_split(row, out_r1, out_r2):
    ""






    r1_raw, r2_raw, interleaved = raw_paths(row)
    if r1_raw is None:
        raise RuntimeError("raw FASTQ files missing")
    f1 = open_fastq(r1_raw)
    f2 = f1 if interleaved else open_fastq(r2_raw)
    out1 = open(out_r1, "w")
    out2 = open(out_r2, "w")
    pairs_in = 0
    pairs_renamed = 0
    invalid_umi = 0
    mate_mismatch = 0
    try:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            r1 = [h1] + [f1.readline() for _ in range(3)]
            r2 = [f2.readline() for _ in range(4)]
            if not r1[3] or not r2[3]:
                raise RuntimeError("file ends in the middle of pair %d"
                                   % (pairs_in + 1))
            pairs_in += 1
            id1 = r1[0].split()[0].lstrip("@")
            id2 = r2[0].split()[0].lstrip("@")
            if id1 != id2:
                mate_mismatch += 1
                continue
            umi = umi_of(r1[0])
            if umi is None or umi_of(r2[0]) != umi:
                invalid_umi += 1
                continue
            new_id = id1 + "_" + umi
            
            
            
            out1.write("@%s UMI=%s\n" % (new_id, umi))
            out2.write("@%s UMI=%s\n" % (new_id, umi))
            out1.write(r1[1])
            out1.write("+\n")
            out1.write(r1[3])
            out2.write(r2[1])
            out2.write("+\n")
            out2.write(r2[3])
            pairs_renamed += 1
    finally:
        f1.close()
        if f2 is not f1:
            f2.close()
        out1.close()
        out2.close()
    
    for p in (out_r1, out_r2):
        with open(p) as f:
            first = f.readline()
        if not first.startswith("@"):
            raise RuntimeError("renamed file %s is not valid FASTQ" % p)
    return pairs_in, pairs_renamed, invalid_umi, mate_mismatch


def run_cutadapt(row, in_r1, in_r2, out_r1, out_r2, too_short_r1,
                 too_short_r2):
    ""
    cmd = [cutadapt_exe,
           "-a", r1_adapter, "-A", r2_adapter,
           "--pair-filter=any",
           "--minimum-length", str(min_length),
           "--too-short-output", to_norm(too_short_r1),
           "--too-short-paired-output", to_norm(too_short_r2),
           "-o", to_norm(out_r1), "-p", to_norm(out_r2),
           to_norm(in_r1), to_norm(in_r2)]
    log_path = in_r1 + ".cutadapt.log"
    with open(log_path, "w") as logf:
        result = subprocess.run(cmd, stdout=logf, stderr=logf)
    err = open(log_path).read()
    if result.returncode != 0:
        raise RuntimeError("cutadapt failed: " + err[-400:])

    def num_after_colon(line):
        
        after_colon = line.split(":")[1]
        first_part = after_colon.strip().split(" ")[0]
        no_commas = first_part.replace(",", "")
        return int(no_commas)

    
    stats = {"pairs_in": 0, "r1_adapter": 0, "r2_adapter": 0,
             "too_short": 0, "pairs_out": 0}
    for line in err.splitlines():
        if "Total read pairs processed" in line:
            stats["pairs_in"] = num_after_colon(line)
        elif "Read 1 with adapter" in line:
            stats["r1_adapter"] = num_after_colon(line)
        elif "Read 2 with adapter" in line:
            stats["r2_adapter"] = num_after_colon(line)
        elif "Pairs that were too short" in line:
            stats["too_short"] = num_after_colon(line)
        elif "Pairs written (passing filters)" in line:
            stats["pairs_out"] = num_after_colon(line)
    if stats["pairs_out"] == 0:
        raise RuntimeError("cutadapt produced 0 pairs - log: %s" % log_path)
    return stats


def to_norm(path):
    return os.path.abspath(path)


def validate_pairing(out_r1, out_r2, n=10000):
    ""
    mismatched = 0
    checked = 0
    with open_fastq(out_r1) as f1, open_fastq(out_r2) as f2:
        for _ in range(n):
            h1 = f1.readline()
            h2 = f2.readline()
            if not h1 or not h2:
                break
            for _ in range(3):
                f1.readline()
                f2.readline()
            checked += 1
            if h1.split()[0] != h2.split()[0]:
                mismatched += 1
    return checked, mismatched


def process_sample(row):
    name = row["sample_id"]
    out_r1 = os.path.join(out_folder, name + "_R1.fastq.gz")
    out_r2 = os.path.join(out_folder, name + "_R2.fastq.gz")
    tmp_r1 = os.path.join(tmp_folder, name + "_renamed_R1.fastq")
    tmp_r2 = os.path.join(tmp_folder, name + "_renamed_R2.fastq")
    short1 = os.path.join(tmp_folder, name + "_tooshort_R1.fastq")
    short2 = os.path.join(tmp_folder, name + "_tooshort_R2.fastq")

    os.makedirs(tmp_folder, exist_ok=True)

    print("  renaming (split + UMI to read name) ...")
    pairs_in, pairs_renamed, invalid_umi, mate_mismatch = \
        rename_and_split(row, tmp_r1, tmp_r2)
    print("  pairs in: %d, renamed: %d, invalid UMI: %d, mate mismatch: %d"
          % (pairs_in, pairs_renamed, invalid_umi, mate_mismatch))
    if mate_mismatch:
        raise RuntimeError("mate IDs did not match - aborting sample")

    print("  cutadapt adapter trimming (easyCLIP adapters, min length %d) ..."
          % min_length)
    stats = run_cutadapt(row, tmp_r1, tmp_r2, out_r1, out_r2, short1, short2)
    print("  cutadapt: pairs in %d, R1 adapter %d, R2 adapter %d, "
          "too short %d, pairs out %d"
          % (stats["pairs_in"], stats["r1_adapter"], stats["r2_adapter"],
             stats["too_short"], stats["pairs_out"]))

    checked, mismatched = validate_pairing(out_r1, out_r2)
    print("  post-preprocessing pairing: %d/%d mismatched"
          % (mismatched, checked))
    if mismatched:
        raise RuntimeError("pairing broken by preprocessing")

    
    for p in (tmp_r1, tmp_r2, short1, short2):
        if os.path.exists(p):
            os.remove(p)

    return {"pairs_in": pairs_in, "pairs_renamed": pairs_renamed,
            "invalid_umi": invalid_umi, "stats": stats,
            "out_r1": out_r1, "out_r2": out_r2,
            "pairing_checked": checked, "pairing_mismatched": mismatched}


def main():
    only = None
    if "--samples" in sys.argv:
        i = sys.argv.index("--samples")
        if i + 1 >= len(sys.argv):
            print("Give a sample ID after --samples")
            sys.exit(1)
        only = sys.argv[i + 1]

    print("=" * 70)
    print("Module 2a: easyCLIP read preprocessing (cutadapt adapter trimming,")
    print("UMI into read name, pairing validation)")
    print("=" * 70)

    if not os.path.exists(cutadapt_exe):
        print("ERROR: cutadapt not found at", cutadapt_exe)
        sys.exit(1)

    rows = load_samples(only)
    os.makedirs(out_folder, exist_ok=True)
    os.makedirs(tmp_folder, exist_ok=True)

    
    existing = {}
    if os.path.exists(summary_file):
        with open(summary_file, newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                existing[r["sample_id"]] = [r[c] for c in
                                            r.keys()]

    summary_rows = []
    failures = []
    for row in rows:
        name = row["sample_id"]
        out_r1 = os.path.join(out_folder, name + "_R1.fastq.gz")
        out_r2 = os.path.join(out_folder, name + "_R2.fastq.gz")
        if os.path.exists(out_r1) and os.path.exists(out_r2) \
                and name in existing:
            print("\nSample %s: outputs and summary row exist, skipping"
                  % name)
            summary_rows.append(existing[name])
            continue
        print("\nSample %s (%s, replicate %s, layout %s)"
              % (name, row["condition"], row["replicate"], row["layout"]))
        try:
            result = process_sample(row)
        except (RuntimeError, OSError) as exc:
            print("  ERROR:", exc)
            failures.append(name)
            summary_rows.append([name, row["condition"], row["replicate"],
                                 row["reads_r1"], "", "", "", "FAILED: "
                                 + str(exc)[:80], "", "FAILED"])
            continue
        if result is None:
            continue

        input_pairs = int(row["reads_r1"])
        retained = result["stats"]["pairs_out"]
        adapter_reads = (result["stats"]["r1_adapter"]
                         + result["stats"]["r2_adapter"])
        status = "OK"
        if retained == 0:
            status = "FAILED"
            failures.append(name)
        summary_rows.append([
            name, row["condition"], row["replicate"],
            input_pairs,                                   
            retained,                                      
            adapter_reads,                                 
            result["pairs_renamed"],                       
            "UMI (7 nt) in headers, verified, appended to read name",
            result["out_r1"] + " ; " + result["out_r2"],
            status])

    with open(summary_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "condition", "biological_replicate",
                         "input_fragments", "fragments_retained",
                         "adapter_trimmed_reads", "barcode_processed_fragments",
                         "umi_status", "output_paths", "status"])
        writer.writerows(summary_rows)
    print("\nPreprocessing summary written to", summary_file)

    if failures:
        print("Samples that failed:", ", ".join(failures))
        sys.exit(1)
    print("Finished preprocessing")


if __name__ == "__main__":
    main()
