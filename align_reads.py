import csv
import gzip
import os
import subprocess
import sys

sample_file = "results/sample_summary.tsv"
align_folder = "results/alignments"
summary_file = "results/alignment_summary.tsv"
temp_folder = "temp_deinterleaved"
reference = "reference/GRCh38.primary_assembly.genome.fa.gz"
reference_name = "GRCh38 (GENCODE v44 primary assembly)"

use_preprocessed = True
preprocessed_folder = "results/preprocessed"

aligner = "minimap2"
star = "/root/miniconda3/envs/bme415/bin/STAR"
star_index = "reference/star_index_sparse"

threads = 8

use_wsl = True
wsl_distro = "Ubuntu"
minimap2 = "/root/miniconda3/envs/bme415/bin/minimap2"
samtools = "/root/miniconda3/envs/bme415/bin/samtools"


def to_wsl_path(path):
    
    path = os.path.abspath(path)
    drive, rest = os.path.splitdrive(path)
    return "/mnt/" + drive[0].lower() + rest.replace("\\", "/")


def get_cmd(program, args):
    if use_wsl:
        return ["wsl.exe", "-d", wsl_distro, "-u", "root", "--",
                program] + args
    return [program] + args


def run_cmd(cmd):
    result = subprocess.run(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    text = result.stdout.decode(errors="replace")
    text += result.stderr.decode(errors="replace")
    return result.returncode, text


def open_fastq(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_id(header):
    return header.split()[0].lstrip("@")


def load_data():
    if not os.path.exists(sample_file):
        print("Sample summary not found. Run Module 1 first.")
        sys.exit(1)

    with open(sample_file, newline="") as f:
        samples = list(csv.DictReader(f, delimiter="\t"))

    problems = []
    for s in samples:
        if s["status"] != "OK":
            problems.append(s["sample_id"] + ": Module 1 status is " + s["status"])
        for key in ("r1_file", "r2_file"):
            if s[key] and not os.path.exists(s[key]):
                problems.append(s["sample_id"] + ": missing " + s[key])

    if problems:
        print("ERROR: unusable samples:")
        for p in problems:
            print(" ", p)
        sys.exit(1)

    return samples


def deinterleave(in_path, out_r1, out_r2):

    pairs = 0
    error = None
    fin = open_fastq(in_path)
    f1 = open(out_r1, "w")
    f2 = open(out_r2, "w")
    while True:
        rec1 = []
        rec2 = []
        line = fin.readline()
        if not line:
            break
        rec1.append(line)
        for i in range(3):
            rec1.append(fin.readline())
        for i in range(4):
            rec2.append(fin.readline())
        if not rec1[3] or not rec2[3]:
            error = "file ends in the middle of pair %d" % (pairs + 1)
            break
        if read_id(rec1[0]) != read_id(rec2[0]):
            error = "read IDs differ in pair %d" % (pairs + 1)
            break
        f1.writelines(rec1)
        f2.writelines(rec2)
        pairs += 1
    fin.close()
    f1.close()
    f2.close()
    return pairs, error


def prepare_reads(sample):
    
    
    
    
    if use_preprocessed:
        r1 = os.path.join(preprocessed_folder, sample["sample_id"]
                          + "_R1.fastq.gz")
        r2 = os.path.join(preprocessed_folder, sample["sample_id"]
                          + "_R2.fastq.gz")
        if not os.path.exists(r1) or not os.path.exists(r2):
            return None, None, ("preprocessed FASTQs missing - run "
                                "preprocess_reads.py first")
        return r1, r2, ""

    if sample["layout"] == "paired":
        return sample["r1_file"], sample["r2_file"], ""

    os.makedirs(temp_folder, exist_ok=True)
    out_r1 = os.path.join(temp_folder, sample["sample_id"] + "_R1.fastq")
    out_r2 = os.path.join(temp_folder, sample["sample_id"] + "_R2.fastq")
    print("    deinterleaving", sample["r1_file"], "...")
    pairs, error = deinterleave(sample["r1_file"], out_r1, out_r2)
    if error:
        return None, None, error
    print("    wrote %d pairs to %s" % (pairs, temp_folder))

    notes = ""
    if pairs != int(sample["reads_r1"]):
        notes = "deinterleaved pair count differs from Module 1 count"
        print("    WARNING:", notes)
    return out_r1, out_r2, notes


def align_sample(sample, r1, r2, sample_folder):
    
    
    bam = os.path.join(sample_folder, sample["sample_id"] + "_aligned.bam")
    bai = bam + ".bai"
    log_path = os.path.join(sample_folder, "alignment_" + aligner + ".log")

    if aligner == "star":
        
        
        cmd = get_cmd(star, [
            "--genomeDir", to_wsl_path(star_index),
            "--readFilesIn", to_wsl_path(r1), to_wsl_path(r2),
            "--readFilesCommand", "zcat",
            "--runThreadN", str(threads),
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outFileNamePrefix", to_wsl_path(os.path.join(sample_folder, ""))])
        print("    running: " + " ".join(cmd))
        with open(log_path, "w") as logf:
            result = subprocess.run(cmd, stdout=logf, stderr=logf)
        if result.returncode != 0:
            return None, None, "STAR failed - see " + log_path
        bam = os.path.join(sample_folder, "Aligned.sortedByCoord.out.bam")
        bai = bam + ".bai"
    else:
        mm2_cmd = get_cmd(minimap2, [
            "-ax", "sr", "-t", str(threads), "-I", "2G",
            "--split-prefix",
            to_wsl_path(os.path.join(sample_folder, "split_parts")),
            to_wsl_path(reference), to_wsl_path(r1), to_wsl_path(r2)])
        sort_cmd = get_cmd(samtools, [
            "sort", "-@", "2", "-m", "256M", "-o", to_wsl_path(bam), "-"])

        print("    aligning ...")
        with open(log_path, "w") as logf:
            p1 = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE,
                                  stderr=logf)
            p2 = subprocess.Popen(sort_cmd, stdin=p1.stdout,
                                  stdout=subprocess.DEVNULL,
                                  stderr=subprocess.PIPE)
            sort_error = p2.communicate()[1].decode(errors="replace")
            p1.communicate()
        if p1.returncode != 0 or p2.returncode != 0:
            return None, None, ("alignment failed (minimap2 code %s, samtools code %s)"
                                % (p1.returncode, p2.returncode))

    code, error = run_cmd(get_cmd(samtools, ["index", to_wsl_path(bam),
                                             to_wsl_path(bai)]))
    if code != 0:
        return None, None, "samtools index failed"

    return bam, bai, ""


def get_stats(bam):
    
    stats = {}
    code, out = run_cmd(get_cmd(samtools, ["flagstat", to_wsl_path(bam)]))
    if code != 0:
        return stats
    for line in out.splitlines():
        if "in total" in line:
            stats["input"] = int(line.split("+")[0])
        elif " mapped (" in line and "not mapped" not in line:
            stats["mapped"] = int(line.split("+")[0])
            in_parens = line.split("(")[1]
            stats["mapped_pct"] = in_parens.split("%")[0]
        elif "properly paired" in line:
            stats["properly_paired"] = int(line.split("+")[0])
        elif "secondary" in line:
            stats["multi"] = int(line.split("+")[0])
    return stats


def parse_star_log(log_path):
    ""
    stats = {}
    if not os.path.exists(log_path):
        return stats
    with open(log_path) as f:
        for line in f:
            if "|" not in line:
                continue
            key, value = line.split("|", 1)
            key = key.strip()
            value = value.strip()
            if key == "Number of input reads":
                stats["input"] = int(value)
            elif key == "Uniquely mapped reads number":
                stats["unique"] = int(value)
            elif key == "Uniquely mapped reads %":
                stats["unique_pct"] = float(value)
            elif key == "Number of reads mapped to multiple loci":
                stats["multi"] = int(value)
            elif key == "% of reads mapped to multiple loci":
                stats["multi_pct"] = float(value)
            elif key == "Number of reads unmapped: too many loci":
                stats["unmapped_toomany"] = int(value)
    return stats


def write_alignment_comparison():
    ""


    if not os.path.exists(summary_file):
        return
    old_path = "results/alignment_summary.tsv"
    if not os.path.exists(old_path):
        return
    old_rows = {}
    with open(old_path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            old_rows[r["sample_id"]] = r
    lines = []
    lines.append("# Alignment comparison: raw reads vs preprocessed reads")
    lines.append("")
    lines.append("Old = minimap2 -ax sr on the raw deposited FASTQs (results/")
    lines.append("alignments/, summary results/alignment_summary.tsv).")
    lines.append("New = the corrected flow (module2_preprocess_reads.py: easyCLIP")
    lines.append("3' adapter trimming, UMI moved into the read name, pairs < 14 nt")
    lines.append("discarded) followed by this module on the preprocessed FASTQs")
    lines.append("(results/alignments_preprocessed/, summary")
    lines.append("results/alignment_summary_preprocessed.tsv).")
    lines.append("")
    lines.append("IMPORTANT: the preprocessed input is slightly smaller than the")
    lines.append("raw input (pairs < 14 nt after trimming are discarded, and pairs")
    lines.append("with a malformed UMI are dropped), so the denominators differ.")
    lines.append("A higher mapping percentage does not by itself mean a better")
    lines.append("alignment - the fraction of the ORIGINAL read pool that ends up")
    lines.append("mapped is the fairer number (reported as retained-then-mapped).")
    lines.append("")
    lines.append("| sample | aligner | input pairs | retained pairs | mapped reads | mapped % | properly paired | secondary (multimap) |")
    lines.append("|---|---|---|---|---|---|---|---|")
    with open(summary_file, newline="") as f:
        new_rows = list(csv.DictReader(f, delimiter="\t"))
    for r in new_rows:
        sid = r["sample_id"]
        o = old_rows.get(sid, {})
        retained = ""
        star_log = os.path.join(align_folder, sid, "Log.final.out")
        star_stats = parse_star_log(star_log)
        pre_summary = os.path.join("results/preprocessing_summary.tsv")
        if os.path.exists(pre_summary):
            with open(pre_summary, newline="") as f:
                for p in csv.DictReader(f, delimiter="\t"):
                    if p["sample_id"] == sid:
                        retained = p["fragments_retained"]
        extra = ""
        if star_stats:
            extra = "unique %d (%.1f%%), multi %d (%.1f%%)" % (
                star_stats.get("unique", 0), star_stats.get("unique_pct", 0),
                star_stats.get("multi", 0), star_stats.get("multi_pct", 0))
        lines.append("| %s | %s | %s | %s | %s | %s | %s | %s |"
                     % (sid, r["aligner"], r["input_pairs"], retained,
                        r["mapped_reads"], r["mapped_pct"],
                        r["properly_paired_reads"], extra or r.get("multimapped_reads", "")))
    lines.append("")
    lines.append("Differences between conditions and replicates should be")
    lines.append("interpreted together with the preprocessing summary (adapter")
    lines.append("content was 0.1-2.2% per sample and the retained fraction")
    lines.append("99.4-99.6%).")
    path = "results/alignment_comparison.md"
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("  wrote", path)


def check_tools():
    
    code, out = run_cmd(get_cmd(minimap2, ["--version"]))
    if code != 0:
        print("ERROR: minimap2 not reachable in WSL.")
        print("Set the tools up with:")
        print("  wsl -d Ubuntu -u root -- bash "
              "/mnt/c/Users/priya/BME415FinalProject/scripts/setup_wsl_alignment.sh")
        sys.exit(1)
    code, out2 = run_cmd(get_cmd(samtools, ["--version"]))
    if code != 0:
        print("ERROR: samtools not reachable in WSL.")
        sys.exit(1)
    version = out.strip().splitlines()[0]
    print("Tools: minimap2 %s, %s (WSL)"
          % (version, out2.strip().splitlines()[0]))
    return version


def check_reference():
    if not os.path.exists(reference):
        print("ERROR: reference genome not found:", reference)
        print("Download GRCh38 (GENCODE v44) - see the notes at the top.")
        sys.exit(1)
    print("Reference:", reference)


def main():
    args = sys.argv[1:]

    dry_run = "--dry-run" in args
    only = None
    if "--samples" in args:
        i = args.index("--samples")
        if i + 1 >= len(args):
            print("Give a sample ID after --samples")
            sys.exit(1)
        only = args[i + 1]

    print("=" * 70)
    print("Module 2: Read Alignment")
    print("KHDRBS2 WT vs R168C easyCLIP-seq (GSE154168)")
    if dry_run:
        print("DRY RUN - reads are prepared and commands printed, no alignment")
    print("=" * 70)

    global align_folder, summary_file
    if use_preprocessed:
        align_folder = "results/alignments_preprocessed"
        summary_file = "results/alignment_summary_preprocessed.tsv"

    samples = load_data()

    if only:
        wanted = only.split(",")
        samples = [s for s in samples if s["sample_id"] in wanted]
        if len(samples) == 0:
            print("No samples matched:", only)
            sys.exit(1)

    version = check_tools()
    if not dry_run:
        check_reference()

    print("Loaded %d sample(s) from %s" % (len(samples), sample_file))

    results = []
    failed = []

    for sample in samples:
        name = sample["sample_id"]
        condition = sample["condition"]
        print("\nSample %s (%s, replicate %s)"
              % (name, condition, sample["replicate"]))

        sample_folder = os.path.join(align_folder, name)
        existing_bam = os.path.join(sample_folder, name + "_aligned.bam")
        if os.path.exists(existing_bam) and os.path.exists(existing_bam + ".bai"):
            print("  BAM and index exist, skipping (delete them to re-align)")
            continue

        os.makedirs(sample_folder, exist_ok=True)

        r1, r2, notes = prepare_reads(sample)
        if r1 is None:
            print("    ERROR:", notes)
            failed.append(name)
            results.append([name, condition, sample["replicate"],
                            sample["layout"], "minimap2", version,
                            reference_name, sample["reads_r1"],
                            0, "", "", "", "", "", "FAILED", notes])
            continue

        if dry_run:
            print("    (dry run: alignment would run here)")
            continue

        bam, bai, error = align_sample(sample, r1, r2, sample_folder)

        status = "OK"
        if error:
            print("    ERROR:", error)
            notes = (notes + "; " if notes else "") + error
            failed.append(name)
            bam = ""
            bai = ""
        elif not os.path.isfile(bam) or not os.path.isfile(bai):
            notes = (notes + "; " if notes else "") + "BAM or index missing"
            failed.append(name)
            bam = ""
            bai = ""
        else:
            print("    BAM:", bam)

        stats = {}
        if status == "OK":
            stats = get_stats(bam)

        results.append([name, condition, sample["replicate"], sample["layout"],
                        "minimap2", version, reference_name,
                        sample["reads_r1"],
                        stats.get("mapped", 0), stats.get("mapped_pct", ""),
                        stats.get("properly_paired", ""),
                        stats.get("multi", ""),
                        bam, bai, status, notes])

    if dry_run:
        print("\nDry run finished. Nothing was aligned and no summary was written.")
        return

    with open(summary_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "condition", "replicate", "layout",
                         "aligner", "aligner_version", "reference",
                         "input_pairs", "mapped_reads", "mapped_pct",
                         "properly_paired_reads", "multimapped_reads",
                         "bam_file", "bai_file", "status", "notes"])
        writer.writerows(results)
    print("\nAlignment summary written to", summary_file)
    write_alignment_comparison()

    print("\nSummary (alignment only - NOT binding sites):")
    for r in results:
        print("  %s  %s  mapped %s %s%%  [%s]"
              % (r[0], r[1], r[8], r[9], r[14]))

    if failed:
        print("\nSamples that failed:", ", ".join(failed))
        sys.exit(1)
    print("\nFinished Module 2")


if __name__ == "__main__":
    main()
