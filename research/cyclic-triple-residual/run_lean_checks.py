"""Compile and independently kernel-replay the research reduction.

This is not the protected instance verifier: no HasCertificate theorem is
claimed. It uses the repository's pinned toolchain, manifest, and options.
Only the source and textual results belong in git; build files are temporary.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import tempfile
import time
import tomllib

ROOT=Path(__file__).resolve().parent


def main(formal_root,lake,output):
    formal_root=formal_root.resolve()
    lake=Path(lake).resolve()
    env=os.environ.copy()
    env['PATH']=str(lake.parent)+os.pathsep+env.get('PATH','')
    def run(args,extra_env=None):
        child_env=env.copy()
        if extra_env: child_env.update(extra_env)
        result=subprocess.run([str(x) for x in args],cwd=formal_root,env=child_env,
                              capture_output=True,text=True)
        if result.returncode:
            raise RuntimeError(result.stdout+'\n'+result.stderr)
        return result.stdout+result.stderr
    version=run([lake.parent/'lean','--version']).strip()
    assert 'version 4.36.0,' in version
    assert (formal_root/'lean-toolchain').read_text().strip()=='leanprover/lean4:v4.36.0'
    manifest=json.loads((formal_root/'lake-manifest.json').read_text())
    mathlib=next(p['rev'] for p in manifest['packages'] if p['name']=='mathlib')
    assert mathlib=='fc4f59e577feb1a4956671bc04ac18aaaf8cd842'
    options=tomllib.loads((formal_root/'lakefile.toml').read_text())['leanOptions']
    flags=[f'-D{k}={str(v).lower() if isinstance(v,bool) else v}' for k,v in options.items()]
    lean_path=run([lake,'env','printenv','LEAN_PATH']).strip()
    source=ROOT/'lean/CyclicCubicReduction.lean'
    with tempfile.TemporaryDirectory(prefix='cyclic-lean-',dir=ROOT) as scratch:
        scratch=Path(scratch)
        start=time.monotonic()
        compile_log=run([lake.parent/'lean',*flags,'-R',source.parent,'-o',
                         scratch/'CyclicCubicReduction.olean',source],{'LEAN_PATH':lean_path})
        compile_seconds=time.monotonic()-start
        assert 'warning:' not in compile_log and 'sorryAx' not in compile_log, compile_log
        axiom_rows=re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]",compile_log)
        assert len(axiom_rows)==9
        allowed={'propext','Classical.choice','Quot.sound'}
        for theorem,axioms in axiom_rows:
            assert set(x.strip() for x in axioms.split(',')) <= allowed, theorem
        start=time.monotonic()
        replay_log=run([lake.parent/'leanchecker','--verbose','CyclicCubicReduction'],
                       {'LEAN_PATH':str(scratch)+os.pathsep+lean_path})
        replay_seconds=time.monotonic()-start
    report={'toolchain':version,'mathlib_revision':mathlib,
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'compile':'PASS; no warnings','kernel_replay':'PASS',
            'axioms':{name:[x.strip() for x in axioms.split(',')] for name,axioms in axiom_rows},
            'compile_seconds':round(compile_seconds,3),'kernel_replay_seconds':round(replay_seconds,3),
            'peak_child_rss_kib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
            'scope':'Universal cubic compact reduction, sixth-root connection trivialization, and conditional current reconstruction only. No target incidence binding, rational scalar primitive, HasCertificate theorem, or protected instance verification.'}
    output.write_text(json.dumps(report,indent=2)+'\n')
    (output.parent/'lean-compile.log').write_text(compile_log)
    (output.parent/'lean-kernel-replay.log').write_text(replay_log)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--formal-root',type=Path,default=ROOT.parents[1]/'formal/lean')
    ap.add_argument('--lake',default=shutil.which('lake'))
    ap.add_argument('--output',type=Path,default=ROOT/'lean-checks.json')
    args=ap.parse_args()
    if not args.lake: ap.error('Provide --lake pointing to the pinned Lean 4.36.0 toolchain')
    main(args.formal_root,args.lake,args.output)
