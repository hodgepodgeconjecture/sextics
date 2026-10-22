"""Mirror scripts/proposals/sandbox.verify without the Docker container.

Same steps and flags: build the trusted imports with Lake, generate the
protected target from instance.rec, compile every submitted module in import
order with `lean -j 2 -M <memoryMb>`, replay each module with a fresh
`leanchecker --verbose` process, then run the trusted Audit.lean.
Reports wall time and peak RSS (VmHWM) per step.

usage: python3 verify_local.py <instance-id> <proof-folder> <lean-project> <out-dir>
(<lean-project> is a copy of formal/lean with the Mathlib cache; the default
formal/lean also works.)
"""
import json, os, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]  # repository root
sys.path.insert(0, str(ROOT/'scripts/proposals'))
import yaml  # noqa: E402
from check import module_order, target_source  # noqa: E402

instance_id, folder, lean_root, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
out.mkdir(parents=True, exist_ok=True)
policy = json.loads((ROOT/'.github/solve-policy.json').read_text())
instance = yaml.safe_load((ROOT/'data'/instance_id/'instance.rec').read_text())
files = {p.name: p.read_bytes() for p in folder.iterdir() if p.suffix == '.lean' or p.name in ('proof.json', 'solution.md')}
manifest = json.loads(files['proof.json'])
order, imports = module_order(files)
elan = Path.home()/'.elan/bin'
env0 = dict(os.environ, PATH=f'{elan}:{os.environ["PATH"]}')


def run(args, label, cwd, env):
    start = time.time(); peak = 0
    with (out/f'{label}.log').open('wb') as log:
        p = subprocess.Popen(args, cwd=cwd, env=env, stdout=log, stderr=subprocess.STDOUT)
        while p.poll() is None:
            try:
                for line in Path(f'/proc/{p.pid}/status').read_text().splitlines():
                    if line.startswith('VmHWM:'):
                        peak = max(peak, int(line.split()[1]))
            except OSError:
                pass
            time.sleep(0.2)
    wall = time.time()-start
    print(f'{label:40s} rc={p.returncode} {wall:8.1f}s peak={peak/1048576:5.2f}GiB', flush=True)
    if p.returncode != 0:
        print((out/f'{label}.log').read_text()[-3000:])
        raise SystemExit(f'{label} failed')
    return wall, peak


subprocess.run(['lake', 'build', 'HodgeComputer.Verify', *imports], cwd=lean_root, env=env0, check=True,
               stdout=subprocess.DEVNULL)
lean_path = subprocess.check_output(['lake', 'env', 'printenv', 'LEAN_PATH'], cwd=lean_root, env=env0).decode().strip()
prefix = Path(subprocess.check_output(['lake', 'env', 'lean', '--print-prefix'], cwd=lean_root, env=env0).decode().strip())
lean, checker = prefix/'bin/lean', prefix/'bin/leanchecker'

work = out/'tree'
sources, artifacts, protected, checked = work/'source', work/'artifacts', work/'protected', work/'checked'
for d in (sources/'Proposal', artifacts/'Proposal', protected, checked):
    d.mkdir(parents=True, exist_ok=True)
for name, raw in files.items():
    if name.endswith('.lean'):
        (sources/'Proposal'/name).write_bytes(raw)
(sources/'HodgeTarget.lean').write_text(target_source(instance))
(sources/'HodgeChecked.lean').write_text(
    'import HodgeTarget\nimport Proposal.Solution\n'
    'theorem HodgeSubmission.checked : HodgeComputer.HasCertificate '
    'HodgeSubmission.F HodgeSubmission.A HodgeSubmission.q := '+manifest['theorem']+'\n')
env = {'PATH': str(prefix/'bin')+':/usr/bin:/bin', 'HOME': '/nonexistent', 'LEAN_SYSROOT': str(prefix),
       'LEAN_PATH': lean_path+':'+':'.join(map(str, (protected, checked, artifacts))), 'LANG': 'en_US.UTF-8'}
protected_names = ['HodgeTarget.lean']
sequence = [*protected_names, *order, 'HodgeChecked.lean']
totals = {'compile': [0, 0], 'replay': [0, 0]}
for file in sequence:
    dest = protected if file in protected_names else checked if file == 'HodgeChecked.lean' else artifacts
    rel = 'Proposal/'+file if dest == artifacts else file
    w, pk = run([str(lean), '-j', '2', '-M', str(policy['memoryMb']), '-o', str(dest/rel.replace('.lean', '.olean')), rel],
                'compile-'+file.removesuffix('.lean'), sources, env)
    totals['compile'][0] += w; totals['compile'][1] = max(totals['compile'][1], pk)
for file in sequence:
    stem = file.removesuffix('.lean')
    module = stem if file in protected_names or file == 'HodgeChecked.lean' else 'Proposal.'+stem
    w, pk = run([str(checker), '--verbose', module], 'kernel-'+stem, sources, env)
    totals['replay'][0] += w; totals['replay'][1] = max(totals['replay'][1], pk)
run([str(lean), '--run', str(ROOT/'scripts/proposals/Audit.lean')], 'theorem-audit', sources, env)
audit = json.loads((out/'theorem-audit.log').read_text())
src = sum(len(raw) for name, raw in files.items() if name.endswith('.lean'))
summary = {'instance': instance_id, 'modules': len(order), 'leanSourceBytes': src,
           'compileSeconds': round(totals['compile'][0], 1), 'compilePeakGiB': round(totals['compile'][1]/1048576, 2),
           'replaySeconds': round(totals['replay'][0], 1), 'replayPeakGiB': round(totals['replay'][1]/1048576, 2),
           'audit': audit}
(out/'summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
