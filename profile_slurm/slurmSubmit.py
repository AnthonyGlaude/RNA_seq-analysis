#!/usr/bin/env python3

import subprocess
import sys
from snakemake.utils import read_job_properties

jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)

cmd = ["sbatch"]

for param, val in job_properties["cluster"].items():
    cmd.extend([f"--{param}", str(val)])

dependencies = [d for d in sys.argv[1:-1] if d]
if dependencies:
    cmd.append(f"--dependency=afterok:{':'.join(dependencies)}")

cmd.append(jobscript)

res = subprocess.run(cmd, capture_output=True, text=True)

if res.returncode != 0:
    sys.stderr.write(res.stderr)
    sys.exit(res.returncode)

# Ex: "Submitted batch job 12345678"
jobid = res.stdout.strip().split()[-1]
print(jobid)