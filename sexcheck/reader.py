import os
import json
import numpy as np
import subprocess
import tempfile
import logging
from typing import Tuple, Dict, List

logger = logging.getLogger(__name__)

class VCFData:
    def __init__(self, samples: List[str], n_samples: int, n_snps: int):
        self.samples = samples
        self.n_samples = n_samples
        self.n_snps = n_snps
        self.variants_id = []
        self.variants_ref = []

def _ensure_processor_compiled() -> str:
    native_dir = os.path.join(os.path.dirname(__file__), "native")
    os.makedirs(native_dir, exist_ok=True)
    source = os.path.join(native_dir, "vcf_reader.c")
    binary = os.path.join(native_dir, "vcf_reader")

    if os.path.exists(binary) and os.path.exists(source):
        if os.path.getmtime(binary) >= os.path.getmtime(source):
            return binary

    if not os.path.exists(source):
        raise RuntimeError(f"No se encuentra {source}")

    cmd = ["gcc", "-O3", "-march=native", "-flto", "-funroll-loops", source, "-o", binary, "-lz"]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Fallo compilando vcf_reader.c: {e.stderr}") from e
    return binary

def read_vcf(vcf_path: str) -> Tuple[np.ndarray, VCFData, Dict]:
    """
    Fast VCF processing using C implementation.
    
    Directly processes VCF file into a normalized zygosity histogram.
    
    Parameters:
        vcf_path: Path to VCF or VCF.gz file
    
    Returns:
        histogram: Histogram data (n_samples, len(mode))
        vcf_data: VCFData object with sample info
        stats: Processing statistics
    """
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(vcf_path)

    processor = _ensure_processor_compiled()

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as fcsv:
        data_file = fcsv.name
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as fjs:
        stats_file = fjs.name

    try:
        mode = "zot"
        normalize = "1"
        if vcf_path.endswith(".gz"):
            cmd = [processor, mode, normalize, data_file, stats_file, vcf_path]
            proc = subprocess.run(cmd, capture_output=True)
        else:
            cmd = [processor, mode, normalize, data_file, stats_file]
            with open(vcf_path, "r") as fin:
                proc = subprocess.run(cmd, stdin=fin, capture_output=True)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr.decode("utf-8") if proc.stderr else "Error ejecutando processor")

        with open(data_file, "r") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if len(lines) < 2:
            raise ValueError("Sin datos tras procesar el VCF")

        header = lines[0].split(",")
        samples, rows = [], []
        for ln in lines[1:]:
            parts = ln.split(",")
            if len(parts) != len(header):
                continue
            samples.append(parts[0])
            rows.append([float(x) for x in parts[1:]])
        hist = np.array(rows, dtype=np.float32)

        # Swap logic:
        # If haploid (#2==0): if #1>#0, swap #0 and #1
        # Else (diploid): if #2>#0, swap #0 and #2
        for h in hist:
            if h[2] == 0:
                if h[1] > h[0]:
                    h[0], h[1] = h[1], h[0]
            else:
                if h[2] > h[0] and h[0] != 0:
                    h[0], h[2] = h[2], h[0]

        # Stats JSON
        try:
            with open(stats_file, "r") as f:
                stats = json.load(f)
        except Exception:
            stats = {"num_samples": len(samples), "num_snps_used": 0, "mode": mode}

        vcf_data = VCFData(samples=samples, n_samples=len(samples), n_snps=stats.get("num_snps_used", 0))
        return hist, vcf_data, stats

    finally:
        try: os.unlink(stats_file)
        except: pass
