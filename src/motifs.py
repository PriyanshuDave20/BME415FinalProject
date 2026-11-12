import csv
import glob
import gzip
import math
import os
import re
import subprocess
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import load_prepare as m1
import find_binding_regions as m3

streme = "/root/miniconda3/envs/bme415/bin/streme"
winlen = 100
minw = 4
maxw = 12
ncond = 5
nrep = 1
sseed = 154168
simcut = 0.85
npar = 3
ref = "reference/GRCh38.primary_assembly.genome.fa.gz"
outdir = "results/motifs"
figdir = "results/motifs/figures"
fadir = "results/motifs/fasta"
dsdir = "results/downsampling/regions_downsampled"
condfile = "results/comparison/region_fold_changes.tsv"


def wslp(p):
    p = os.path.abspath(p)
    d, rest = os.path.splitdrive(p)
    return "/mnt/" + d[0].lower() + rest.replace("\\", "/")


def loadcond():
    wt = []
    r = []
    f = open(condfile, newline="")
    for x in csv.DictReader(f, delimiter="\t"):
        d = {"id": x["region_id"], "chrom": x["chrom"],
             "start": int(x["start"]), "end": int(x["end"])}
        if x["region_id"].startswith("WT_"):
            wt.append(d)
        else:
            r.append(d)
    f.close()
    return wt, r


def getbr():
    brs = {
        "WT-BR1": "SRR12192530", "WT-BR2": "SRR12192532",
        "WT-BR3": "SRR12192534", "WT-BR4": "SRR12192535",
        "R168C-BR1": "SRR12192536", "R168C-BR2": "SRR12192538",
        "R168C-BR3": "SRR12192540", "R168C-BR4": "SRR12192541",
    }
    out = {}
    for nm in sorted(brs):
        bam = os.path.join("results", "dedup", nm + ".bam")
        cnt, strands, tot = m3.count_reads(bam)
        thr, lam = m3.get_threshold(cnt)
        found = m3.find_regions(cnt, strands, thr)
        cond = "WT" if nm.startswith("WT") else "R168C"
        regs = []
        for x in found:
            regs.append({"id": "%s_%05d" % (nm, len(regs) + 1),
                         "sample_id": nm, "condition": cond,
                         "chrom": x["chrom"], "start": x["start"],
                         "end": x["end"]})
        out[nm] = regs
        print("  %s: %d regions" % (nm, len(regs)))
    return out


def getds():
    out = {}
    for path in glob.glob(os.path.join(dsdir, "*_regions.tsv")):
        nm = os.path.basename(path).replace("_regions.tsv", "")
        nm = nm.replace("_", "-", 1)
        cond = "WT" if nm.startswith("WT") else "R168C"
        regs = []
        f = open(path, newline="")
        for x in csv.DictReader(f, delimiter="\t"):
            regs.append({"id": "%s_%05d" % (nm, len(regs) + 1),
                         "sample_id": nm, "condition": cond,
                         "chrom": x["chrom"], "start": int(x["start"]),
                         "end": int(x["end"])})
        f.close()
        out[nm] = regs
    return out


def mergeit(regions, cond):
    cov = {}
    for r in regions:
        if r["condition"] != cond:
            continue
        f = r["start"] // 4000
        l = (r["end"] - 1) // 4000
        for b in range(f, l + 1):
            cov[(r["chrom"], b)] = True
    bychrom = {}
    for chrom, b in cov:
        bychrom.setdefault(chrom, []).append(b)
    out = []
    for chrom in sorted(bychrom):
        bl = sorted(bychrom[chrom])
        groups = []
        for b in bl:
            if groups and b == groups[-1][-1] + 1:
                groups[-1].append(b)
            else:
                groups.append([b])
        for g in groups:
            out.append({"id": "%s_%05d" % (cond, len(out) + 1),
                        "condition": cond, "chrom": chrom,
                        "start": g[0] * 4000, "end": (g[-1] + 1) * 4000,
                        "bins": g})
    return out


def extract(sets):
    bychrom = {}
    for nm in sets:
        for r in sets[nm]:
            bychrom.setdefault(r["chrom"], []).append((nm, r))
    seqs = {}
    stats = {}
    for nm in sets:
        seqs[nm] = {}
        stats[nm] = {"ex": 0, "dup": 0}

    def dochrom(chrom, parts):
        if chrom is None or chrom not in bychrom:
            return
        seq = "".join(parts)
        half = winlen // 2
        for nm, r in bychrom[chrom]:
            c = r["start"] + (r["end"] - r["start"]) // 2
            s = seq[c - half:c - half + winlen].upper()
            if len(s) < winlen or "N" in s:
                stats[nm]["ex"] += 1
                continue
            if s in seqs[nm]:
                stats[nm]["dup"] += 1
                continue
            seqs[nm][r["id"]] = s

    cur = None
    parts = []
    f = gzip.open(ref, "rt")
    for line in f:
        if line.startswith(">"):
            dochrom(cur, parts)
            cur = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    dochrom(cur, parts)
    f.close()
    return seqs, stats


def writefa(path, seqs, label):
    f = open(path, "w")
    for rid in seqs:
        f.write(">%s\n" % rid)
        s = seqs[rid]
        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + "\n")
    f.close()
    print("  wrote %s (%d %s sequences)" % (path, len(seqs), label))


def countfa(path):
    n = 0
    f = open(path)
    for line in f:
        if line.startswith(">"):
            n = n + 1
    f.close()
    return n


def run1(job, sd=True):
    os.makedirs(job["outdir"], exist_ok=True)
    cmd = ["wsl.exe", "-d", "Ubuntu", "-u", "root", "--", streme,
           "--p", wslp(job["fasta"]), "--oc", wslp(job["outdir"]),
           "--nmotifs", str(job["n"]), "--minw", str(minw),
           "--maxw", str(maxw)]
    if sd:
        cmd = cmd + ["--seed", str(sseed)]
    e = open(os.path.join(job["outdir"], "streme_stderr.log"), "w")
    r = subprocess.run(cmd, stdout=e, stderr=e)
    e.close()
    if r.returncode != 0 and sd:
        t = open(os.path.join(job["outdir"], "streme_stderr.log")).read()[-300:]
        if "seed" in t or "unrecognized" in t:
            return run1(job, False)
        print("  FAILED: %s: %s" % (job["name"], t))
        return False
    return True


def runjobs(jobs):
    bad = []
    todo = []
    for j in jobs:
        j["nseq"] = countfa(j["fasta"])
        if os.path.exists(os.path.join(j["outdir"], "streme.txt")):
            print("  %-18s already done, skipping" % j["name"])
            continue
        todo.append(j)
    running = []
    while todo or running:
        while todo and len(running) < npar:
            j = todo.pop(0)
            print("  starting STREME: %s (%d sequences)" % (j["name"], j["nseq"]))
            os.makedirs(j["outdir"], exist_ok=True)
            cmd = ["wsl.exe", "-d", "Ubuntu", "-u", "root", "--", streme,
                   "--p", wslp(j["fasta"]), "--oc", wslp(j["outdir"]),
                   "--nmotifs", str(j["n"]), "--minw", str(minw),
                   "--maxw", str(maxw), "--seed", str(sseed)]
            e = open(os.path.join(j["outdir"], "streme_stderr.log"), "w")
            p = subprocess.Popen(cmd, stdout=e, stderr=e)
            running.append((j, p, e))
        still = []
        for j, p, e in running:
            if p.poll() is None:
                still.append((j, p, e))
                continue
            e.close()
            if p.returncode == 0:
                print("  finished STREME: %s" % j["name"])
            else:
                t = open(os.path.join(j["outdir"],
                                      "streme_stderr.log")).read()[-300:]
                if "seed" in t or "unrecognized" in t:
                    ok = run1(j, False)
                    if not ok:
                        bad.append(j["name"])
                else:
                    bad.append(j["name"])
                    print("  FAILED: %s: %s" % (j["name"], t))
        running = still
        if running:
            time.sleep(5)
    return bad


def parsestreme(path):
    motifs = []
    cur = None
    f = open(path)
    for line in f:
        if line.startswith("MOTIF"):
            x = line.split()
            con = ""
            mid = x[-1] if len(x) > 1 else "motif"
            if len(x) >= 3:
                con = x[1]
                mid = x[2]
            if "-" in con and con.split("-")[0].isdigit():
                con = con.split("-", 1)[1]
            cur = {"id": mid, "con": con, "w": 0, "nsites": 0,
                   "e": None, "pwm": []}
            motifs.append(cur)
        elif cur is not None and line.startswith("letter-probability"):
            mw = re.search(r"w=\s*(\d+)", line)
            mn = re.search(r"nsites=\s*(\d+)", line)
            me = re.search(r"E=\s*(\S+)", line)
            if mw:
                cur["w"] = int(mw.group(1))
            if mn:
                cur["nsites"] = int(mn.group(1))
            if me:
                cur["e"] = float(me.group(1))
        elif cur is not None:
            x = line.split()
            if len(x) == 4:
                try:
                    cur["pwm"].append([float(v) for v in x])
                except ValueError:
                    continue
                if len(cur["pwm"]) == cur["w"]:
                    cur = None
    f.close()
    return motifs


def seqcounts(tag):
    cnt = {}
    path = os.path.join(outdir, "streme_" + tag, "sequences.tsv")
    if os.path.exists(path):
        pairs = set()
        f = open(path)
        for line in f:
            if line.startswith("motif_ID"):
                continue
            x = line.split("\t")
            if len(x) >= 4:
                pairs.add((x[1], x[3]))
        f.close()
        for p in pairs:
            cnt[p[0]] = cnt.get(p[0], 0) + 1
    return cnt


def rcpwm(pwm):
    out = []
    for row in reversed(pwm):
        out.append([row[3], row[2], row[1], row[0]])
    return out


def pwmsim(a, b):
    best = 0.0
    for ob in (b, rcpwm(b)):
        for off in range(-(len(ob) - 1), len(a)):
            sims = []
            for i in range(len(a)):
                j = i - off
                if 0 <= j < len(ob):
                    ca = a[i]
                    cb = ob[j]
                    dot = ca[0]*cb[0] + ca[1]*cb[1] + ca[2]*cb[2] + ca[3]*cb[3]
                    na = math.sqrt(ca[0]**2 + ca[1]**2 + ca[2]**2 + ca[3]**2)
                    nb = math.sqrt(cb[0]**2 + cb[1]**2 + cb[2]**2 + cb[3]**2)
                    if na > 0 and nb > 0:
                        sims.append(dot / (na * nb))
            if sims:
                m = 0.0
                for s in sims:
                    if s > m:
                        m = s
                if m > best:
                    best = m
    return best


def render(ax, pwm, title):
    bases = ["A", "C", "G", "T"]
    cols = {"A": "green", "C": "blue", "G": "orange", "T": "red"}
    for i in range(len(pwm)):
        y = 0.0
        row = pwm[i]
        order = []
        
        remaining = [0, 1, 2, 3]
        while remaining:
            best_j = remaining[0]
            for j in remaining:
                if row[j] > row[best_j]:
                    best_j = j
            order.append(best_j)
            remaining.remove(best_j)
        for j in order:
            p = row[j]
            if p <= 0.02:
                continue
            ax.text(i + 0.5, y + p / 2, bases[j], fontsize=5 + 13 * p,
                    color=cols[bases[j]], ha="center", va="center")
            y = y + p
    ax.set_xlim(0, len(pwm))
    ax.set_ylim(0, 1.02)
    ax.set_xticks([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=8)


def main():
    print("=" * 70)
    print("Module 5: motif discovery and comparison (STREME)")
    print("=" * 70)

    r = subprocess.run(["wsl.exe", "-d", "Ubuntu", "-u", "root", "--",
                        streme, "--version"], stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    if r.returncode != 0:
        print("ERROR: STREME is not reachable in WSL.")
        sys.exit(1)
    print("STREME is available in WSL.")

    os.makedirs(outdir, exist_ok=True)
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(fadir, exist_ok=True)

    print("Loading condition regions from Module 4 ...")
    wt, r = loadcond()
    print("  WT: %d, R168C: %d" % (len(wt), len(r)))
    print("Per-BR regions from the deduped bams (Module 3 functions) ...")
    brs = getbr()
    print("Downsampled BR regions from Module 4 ...")
    ds = getds()
    for nm in sorted(ds):
        print("  %s (downsampled): %d regions" % (nm, len(ds[nm])))

    ds_wt = []
    ds_r = []
    for nm in brs:
        use = ds.get(nm, brs[nm])
        for x in use:
            if nm.startswith("WT"):
                ds_wt.append(x)
            else:
                ds_r.append(x)
    ds_wtu = mergeit(ds_wt, "WT")
    ds_ru = mergeit(ds_r, "R168C")
    print("  downsampled unions: WT %d, R168C %d" % (len(ds_wtu), len(ds_ru)))

    sets = {"WT": wt, "R168C": r, "DS-WT": ds_wtu, "DS-R168C": ds_ru}
    for nm in brs:
        sets[nm] = brs[nm]

    print("Extracting %d nt windows ..." % winlen)
    seqs, stats = extract(sets)
    print("Sequence extraction results:")
    for nm in sorted(sets):
        print("  %-10s kept %6d (N excluded: %d, duplicates: %d)"
              % (nm, len(seqs[nm]), stats[nm]["ex"], stats[nm]["dup"]))

    writefa(os.path.join(outdir, "wt_sequences.fasta"), seqs["WT"], "WT")
    writefa(os.path.join(outdir, "r168c_sequences.fasta"), seqs["R168C"],
            "R168C")
    writefa(os.path.join(fadir, "ds_wt_sequences.fasta"), seqs["DS-WT"],
            "downsampled WT")
    writefa(os.path.join(fadir, "ds_r168c_sequences.fasta"), seqs["DS-R168C"],
            "downsampled R168C")
    for nm in sorted(brs):
        fnm = nm.replace("-", "_").lower()
        writefa(os.path.join(fadir, fnm + "_sequences.fasta"), seqs[nm], nm)

    f = open(os.path.join(outdir, "extraction_stats.tsv"), "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["sequence_set", "regions", "sequences_kept",
                "excluded_N_or_short", "duplicate_sequences"])
    for nm in sorted(sets):
        w.writerow([nm, len(sets[nm]), len(seqs[nm]), stats[nm]["ex"],
                    stats[nm]["dup"]])
    f.close()
    print("  wrote extraction_stats.tsv")

    jobs = []
    jobs.append({"name": "wt", "fasta": os.path.join(outdir, "wt_sequences.fasta"),
                 "outdir": os.path.join(outdir, "streme_wt"), "n": ncond})
    jobs.append({"name": "r168c",
                 "fasta": os.path.join(outdir, "r168c_sequences.fasta"),
                 "outdir": os.path.join(outdir, "streme_r168c"), "n": ncond})
    jobs.append({"name": "ds_wt", "fasta": os.path.join(fadir, "ds_wt_sequences.fasta"),
                 "outdir": os.path.join(outdir, "streme_ds_wt"), "n": ncond})
    jobs.append({"name": "ds_r168c",
                 "fasta": os.path.join(fadir, "ds_r168c_sequences.fasta"),
                 "outdir": os.path.join(outdir, "streme_ds_r168c"), "n": ncond})
    for nm, tag in (("wt_br1", "WT-BR1"), ("wt_br2", "WT-BR2"),
                    ("wt_br3", "WT-BR3"), ("wt_br4", "WT-BR4"),
                    ("r168c_br1", "R168C-BR1"), ("r168c_br2", "R168C-BR2"),
                    ("r168c_br3", "R168C-BR3"), ("r168c_br4", "R168C-BR4")):
        jobs.append({"name": nm,
                     "fasta": os.path.join(fadir, tag.replace("-", "_").lower()
                                           + "_sequences.fasta"),
                     "outdir": os.path.join(outdir, "streme_" + nm),
                     "n": nrep})
    print("Running %d STREME jobs ..." % len(jobs))
    bad = runjobs(jobs)

    wtm = parsestreme(os.path.join(outdir, "streme_wt", "streme.txt"))
    rm = parsestreme(os.path.join(outdir, "streme_r168c", "streme.txt"))
    wseq = seqcounts("wt")
    rseq = seqcounts("r168c")
    print("WT motifs: %d, R168C motifs: %d" % (len(wtm), len(rm)))

    pairs = []
    for a in wtm:
        for b in rm:
            s = pwmsim(a["pwm"], b["pwm"])
            pairs.append((s, a, b))
    
    pairs.sort(reverse=True)
    usedw = []
    usedr = []
    matched = []
    for p in pairs:
        if p[0] >= simcut and p[1] not in usedw and p[2] not in usedr:
            usedw.append(p[1])
            usedr.append(p[2])
            matched.append(p)

    f = open(os.path.join(outdir, "motif_comparison.tsv"), "w", newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["wt_motif_id", "r168c_motif_id", "wt_consensus",
                "r168c_consensus", "width_wt", "width_r168c", "similarity",
                "wt_evalue", "r168c_evalue", "wt_seqs_with_motif",
                "r168c_seqs_with_motif", "relationship"])
    for p in matched:
        a = p[1]
        b = p[2]
        w.writerow([a["id"], b["id"], a["con"], b["con"], a["w"], b["w"],
                    "%.3f" % p[0], "%.2e" % a["e"], "%.2e" % b["e"],
                    wseq.get(a["id"], 0), rseq.get(b["id"], 0), "shared"])
    for a in wtm:
        if a in usedw:
            continue
        best = 0.0
        for b in rm:
            s = pwmsim(a["pwm"], b["pwm"])
            if s > best:
                best = s
        w.writerow([a["id"], "", a["con"], "", a["w"], "", "%.2e" % a["e"],
                    "", wseq.get(a["id"], 0), "",
                    "WT-associated (best sim to R168C: %.2f)" % best])
    for b in rm:
        if b in usedr:
            continue
        best = 0.0
        for a in wtm:
            s = pwmsim(a["pwm"], b["pwm"])
            if s > best:
                best = s
        w.writerow(["", b["id"], "", b["con"], "", b["w"], "",
                    "%.2e" % b["e"], "", rseq.get(b["id"], 0),
                    "R168C-associated (best sim to WT: %.2f)" % best])
    f.close()
    print("  wrote motif_comparison.tsv")

    f = open(os.path.join(outdir, "replicate_motif_summary.tsv"), "w",
             newline="")
    w = csv.writer(f, delimiter="\t")
    w.writerow(["condition", "motif_id", "consensus", "width", "evalue",
                "seqs_with_motif", "brs_with_match", "br_similarities",
                "recurring"])
    rows = []
    for cond, motifs, tag, seqs2 in (("WT", wtm, "wt", wseq),
                                     ("R168C", rm, "r168c", rseq)):
        sims = {}
        for m in motifs:
            sims[m["id"]] = []
        for i in (1, 2, 3, 4):
            brf = os.path.join(outdir, "streme_%s_br%d" % (tag, i),
                               "streme.txt")
            if not os.path.exists(brf):
                continue
            brm = parsestreme(brf)
            for m in motifs:
                best = 0.0
                for b in brm:
                    s = pwmsim(m["pwm"], b["pwm"])
                    if s > best:
                        best = s
                sims[m["id"]].append("%.2f" % best)
        for m in motifs:
            ss = sims[m["id"]]
            n = 0
            for s in ss:
                if float(s) >= simcut:
                    n = n + 1
            rec = "no"
            if n >= 2:
                rec = "yes"
            rows.append({"condition": cond, "motif_id": m["id"],
                         "consensus": m["con"], "brs_with_match": n,
                         "br_similarities": "/".join(ss) if ss else ""})
            w.writerow([cond, m["id"], m["con"], m["w"],
                        "%.2e" % m["e"], seqs2.get(m["id"], 0), n,
                        rows[-1]["br_similarities"], rec])
    f.close()
    print("  wrote replicate_motif_summary.tsv")

    for m in wtm:
        plt.figure(figsize=(0.7 * m["w"] + 1, 2.2))
        render(plt.gca(), m["pwm"], "WT %s" % m["con"])
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "wt_%s_logo.png" % m["id"]), dpi=150)
        plt.close()
    for m in rm:
        plt.figure(figsize=(0.7 * m["w"] + 1, 2.2))
        render(plt.gca(), m["pwm"], "R168C %s" % m["con"])
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "r168c_%s_logo.png" % m["id"]),
                    dpi=150)
        plt.close()
    for i, p in enumerate(matched[:4], start=1):
        plt.figure(figsize=(1.4 * max(p[1]["w"], p[2]["w"]) + 2, 3))
        render(plt.subplot(1, 2, 1), p[1]["pwm"], "WT: %s" % p[1]["con"])
        render(plt.subplot(1, 2, 2), p[2]["pwm"], "R168C: %s" % p[2]["con"])
        plt.suptitle("Shared motif pair %d (PWM similarity %.2f)" % (i, p[0]))
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "shared_motif_pair_%d.png" % i),
                    dpi=150)
        plt.close()
    if rows:
        plt.figure(figsize=(8, 4.5))
        cols = []
        for x in rows:
            if x["condition"] == "WT":
                cols.append("seagreen")
            else:
                cols.append("firebrick")
        plt.bar(range(len(rows)),
                [int(x["brs_with_match"]) for x in rows], color=cols)
        plt.axhline(2, color="grey", linestyle="--", linewidth=1,
                    label="recurring (>= 2 BRs)")
        plt.xticks(range(len(rows)),
                   ["%s\n%s" % (x["condition"], x["consensus"]) for x in rows],
                   fontsize=7)
        plt.ylabel("biological replicates with the motif")
        plt.title("Motif recurrence across biological replicates")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "replicate_recurrence.png"), dpi=150)
        plt.close()

    rec = []
    for cond, motifs, tag in (("WT", wtm, "ds_wt"),
                              ("R168C", rm, "ds_r168c")):
        dsf = os.path.join(outdir, "streme_" + tag, "streme.txt")
        if not os.path.exists(dsf):
            continue
        dsm = parsestreme(dsf)
        for m in motifs:
            best = 0.0
            for b in dsm:
                s = pwmsim(m["pwm"], b["pwm"])
                if s > best:
                    best = s
            rec.append((cond, m["con"], best))
    if rec:
        plt.figure(figsize=(8, 4.5))
        plt.bar(range(len(rec)), [x[2] for x in rec],
                color=["seagreen" if x[0] == "WT" else "firebrick"
                       for x in rec])
        plt.axhline(simcut, color="grey", linestyle="--", linewidth=1,
                    label="recovered (>= %.2f)" % simcut)
        plt.xticks(range(len(rec)),
                   ["%s\n%s" % (x[0], x[1]) for x in rec], fontsize=7)
        plt.ylabel("best PWM similarity to downsampled motifs")
        plt.title("Original vs downsampled motif recovery")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(figdir, "downsampled_motif_recovery.png"),
                    dpi=150)
        plt.close()

    lines = []
    lines.append("# Module 5: motif discovery summary (STREME)")
    lines.append("")
    lines.append("STREME (MEME Suite 5.5.9) was run on %d nt windows centred on"
                 % winlen)
    lines.append("the condition-level candidate regions from the Module 4")
    lines.append("biological-replicate analysis (deduplicated BAMs, deeper")
    lines.append("HiSeq run per library). Without control sequences STREME")
    lines.append("shuffles the primary sequences preserving dinucleotide")
    lines.append("frequencies, so the background comes from each condition's")
    lines.append("own regions. Settings: --minw %d --maxw %d --nmotifs %d"
                 % (minw, maxw, ncond))
    lines.append("(condition and downsampled runs) / %d (replicate runs),"
                 % nrep)
    lines.append("--seed %d." % sseed)
    lines.append("")
    lines.append("## Sequences analysed")
    lines.append("")
    for nm in sorted(sets):
        lines.append("- %s: %d sequences kept (%d excluded with N, %d duplicates)"
                     % (nm, len(seqs[nm]), stats[nm]["ex"], stats[nm]["dup"]))
    lines.append("")
    lines.append("Windows are region-centred (Module 3 stored no per-base")
    lines.append("summits), so motif signal is diluted by the 4 kb bin width.")
    lines.append("Sequences are genomic DNA; candidate regions are not")
    lines.append("reliably stranded, so no per-strand orientation was used.")
    lines.append("")
    lines.append("## Motifs found")
    lines.append("")
    lines.append("### WT (%d motifs)" % len(wtm))
    lines.append("")
    for m in wtm:
        lines.append("- %s (width %d, %d sites, %d sequences, E = %.2e)"
                     % (m["con"], m["w"], m["nsites"],
                        wseq.get(m["id"], 0), m["e"]))
    lines.append("")
    lines.append("### R168C (%d motifs)" % len(rm))
    lines.append("")
    for m in rm:
        lines.append("- %s (width %d, %d sites, %d sequences, E = %.2e)"
                     % (m["con"], m["w"], m["nsites"],
                        rseq.get(m["id"], 0), m["e"]))
    lines.append("")
    lines.append("## WT vs R168C comparison")
    lines.append("")
    lines.append("Motifs were matched by best-offset PWM cosine similarity in")
    lines.append("both orientations; similarity >= %.2f counts as the same"
                 % simcut)
    lines.append("motif (see motif_comparison.tsv).")
    lines.append("")
    lines.append("- shared motifs: %d" % len(matched))
    for p in matched:
        lines.append("  - %s ~ %s (similarity %.2f)"
                     % (p[1]["con"], p[2]["con"], p[0]))
    for a in wtm:
        if a not in usedw:
            lines.append("- %s listed as WT-associated, but it is a close")
            lines.append("  relative of the shared motif family")
    for b in rm:
        if b not in usedr:
            lines.append("- %s listed as R168C-associated, but it is a close")
            lines.append("  relative of the shared motif family")
    lines.append("")
    lines.append("## Motif recurrence across biological replicates")
    lines.append("")
    lines.append("STREME was run separately on each of the 4 biological")
    lines.append("replicates per condition (1 motif each). A condition motif")
    lines.append("counts as detected in a replicate when the PWM similarity is")
    lines.append(">= %.2f (see replicate_motif_summary.tsv):" % simcut)
    lines.append("")
    for x in rows:
        lines.append("- %s %s: detected in %s of 4 replicates (%s)"
                     % (x["condition"], x["consensus"], x["brs_with_match"],
                        x["br_similarities"].replace("/", ", ")))
    lines.append("")
    lines.append("## Downsampled sensitivity check")
    lines.append("")
    lines.append("STREME was re-run on the Module 4 downsampled region windows")
    lines.append("(equalized sequencing depth within matched replicate pairs):")
    lines.append("")
    for x in rec:
        lines.append("- %s %s: similarity %.2f (%s)"
                     % (x[0], x[1], x[2],
                        "recovered" if x[2] >= simcut else "NOT recovered"))
    lines.append("")
    lines.append("## Interpretation limits")
    lines.append("")
    lines.append("- An enriched motif in candidate regions is not a confirmed")
    lines.append("  KHDRBS2 binding motif.")
    lines.append("- A motif seen in only one condition can reflect candidate-")
    lines.append("  region selection and sequencing depth, not binding")
    lines.append("  specificity.")
    lines.append("- Module 4 showed the two conditions' region sets overlap")
    lines.append("  strongly and their overall k-mer composition is nearly")
    lines.append("  identical, so large motif differences would be surprising.")
    lines.append("- STREME finds enriched patterns; it does not establish")
    lines.append("  biological function or a direct biochemical effect of the")
    lines.append("  R168C mutation.")

    f = open(os.path.join(outdir, "motif_summary.md"), "w")
    f.write("\n".join(lines) + "\n")
    f.close()
    print("  wrote motif_summary.md")

    print("\nFinished Module 5")
    if bad:
        print("STREME jobs that failed:", ", ".join(bad))
        sys.exit(1)


if __name__ == "__main__":
    main()
