"""Paired per-wrap bootstrap (score_ab.paired_bootstrap, unchanged) of arm B in one score JSON against arm A in another
(same truth, same scorer). Also mean per-wrap modal share and hop rate. usage: compare_arms.py <A json> <B json> [boot]"""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'vendor'))
import score_ab as sa
ja, jb = json.load(open(sys.argv[1])), json.load(open(sys.argv[2]))
boot = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
A, B = ja['arms']['A'], jb['arms']['B']
ra, rb = A['per_wrap'], B['per_wrap']
for k in ('within_2', 'within_4'):
    r = sa.paired_bootstrap(ra, rb, k, boot, 12345)
    print(k, f"A {A[k]:.3f} B {B[k]:.3f}  B-A {r['point']*100:+.1f} points, 95% CI {r['ci95'][0]*100:+.1f} to {r['ci95'][1]*100:+.1f}")
for name, X in (('A', A), ('B', B)):
    pw = X['per_wrap']
    n4 = sum(v['within_4'] for v in pw.values())
    print(name, f"adjacent exact {X['adjacent_pairs']['exact_share']:.3f} same {X['adjacent_pairs']['same_winding_share']:.3f} "
          f"hops {X['winding_hops_total']} hop_rate {X['winding_hops_total']/max(n4,1):.4f} "
          f"modal {np.mean([v['modal_share'] for v in pw.values()]):.3f} dup(scorer) {X['duplicate_sheet']['shared_share']:.3f}")
