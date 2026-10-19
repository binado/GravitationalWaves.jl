"""Run independent native Enzyme probes with bounded compilation time.

Example: JULIA='julia +1.12' python3 integration/GWSEOBNR/probe_enzyme.py
Logs and TOML records are saved under gitignored reference/seobnr/enzyme.
"""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time
import tomllib

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'reference/seobnr/enzyme'
DEST.mkdir(parents=True, exist_ok=True)
JULIA = shlex.split(os.environ.get('JULIA', 'julia +1.12'))
CASES = [('calibration', 'default'), ('hamiltonian', 'default'),
         ('flux', 'default'), ('fft', 'default'), ('tail', 'default'),
         ('timing', 'default'), ('timing', 'relaxed_aliasing'),
         ('initial', 'default'), ('initial', 'runtime_activity'),
         ('evolution', 'default'), ('waveform', 'default'),
         ('waveform', 'runtime_activity'), ('likelihood', 'default')]
summary = []
for case, mode in CASES:
    name = f'{case}_{mode}'
    record, log = DEST / (name + '.toml'), DEST / (name + '.log')
    if record.exists():
        record.unlink()
    started = time.monotonic()
    print(f'Starting {name}', flush=True)
    with log.open('w') as stream:
        try:
            completed = subprocess.run(
                JULIA + ['--startup-file=no', '--project=integration/GWSEOBNR',
                         'integration/GWSEOBNR/probe_enzyme.jl', case, mode, str(record)],
                cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, timeout=180)
            status = {'exit_code': completed.returncode}
        except subprocess.TimeoutExpired:
            status = {'status': 'timeout', 'timeout_seconds': 180}
    status.update(case=case, mode=mode, wall_seconds=time.monotonic()-started)
    if record.exists():
        data = tomllib.loads(record.read_text())
        status['status'] = data['status']
        if 'enzyme_forwarddiff_error' in data:
            status['enzyme_forwarddiff_error'] = data['enzyme_forwarddiff_error']
        if 'error_type' in data:
            status['error_type'] = data['error_type']
    else:
        status.setdefault('status', 'process_failure')
    summary.append(status)
    (DEST / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(status, flush=True)

sys.exit(0 if all(item['status'] == 'pass' for item in summary) else 1)
