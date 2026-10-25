import sys, re, math
sys.path.insert(0, '.')
from replay import kernel_calls

res = sys.argv[1] if len(sys.argv) > 1 else 'results.txt'
t = {}; flips = {}
for line in open(res):
    m = re.match(r'(\S+) (?:load=(\S+) )?real (\S+) user (\S+)\s+flips=(\d+)', line)
    if not m:
        continue
    name, ld, real, user, fl = m.groups()
    t[name] = float(user); flips[name] = int(fl)

UNIFORM = re.compile(r'^d\d+$')
depthsets = {'default': {4: 6, 5: 5, 7: 5, 8: 5}, 'legacy': {4: 12, 5: 8, 7: 7, 8: 7}}
rows = []
for v in ['d1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd8', 'd10', 'd12', 'default', 'legacy', 'frozen']:
    t1, t2 = t.get(v + '_n1000'), t.get(v + '_n2000')
    if t1 is None or t2 is None:
        continue
    rate = (t2 - t1) / 1000.0
    if v == 'frozen':
        K = 960; corners = 0
    else:
        dd = int(v[1:]) if UNIFORM.match(v) else None
        ds = depthsets.get(v) or {4: dd, 5: dd, 7: dd, 8: dd}
        Ka, _, ca = kernel_calls(v + '_n1000', ds); Kb, _, cb = kernel_calls(v + '_n2000', ds)
        K = (Ka + Kb) / 2; corners = (ca + cb) / 2
    rows.append(dict(v=v, t1=t1, t2=t2, ms=rate * 1e3, K=K, corners=corners,
                     f1=flips.get(v + '_n1000'), f2=flips.get(v + '_n2000')))
    print('%-8s t1000=%7.2f t2000=%7.2f  ms/step(1000-2000)=%7.2f  steps/s=%5.1f  flips=%d/%d  K/step=%6.0f  EVcorners=%.0f'
          % (v, t1, t2, rate * 1e3, 1 / rate, flips.get(v + '_n1000'), flips.get(v + '_n2000'), K, corners))


def fit(pts):
    n = len(pts); sx = sum(x for x, y in pts); sy = sum(y for x, y in pts)
    sxx = sum(x * x for x, y in pts); sxy = sum(x * y for x, y in pts)
    b = (n * sxy - sx * sy) / (n * sxx - sx * sx); a = (sy - b * sx) / n
    ss_res = sum((y - a - b * x) ** 2 for x, y in pts); ymean = sy / n
    ss_tot = sum((y - ymean) ** 2 for x, y in pts)
    return a, b, 1 - ss_res / ss_tot, math.sqrt(ss_res / n), max(abs(y - a - b * x) / y for x, y in pts)


uni = [r for r in rows if UNIFORM.match(r['v'])]
for dmin in (1, 2, 3):
    pts = [(int(r['v'][1:]), r['ms']) for r in uni if int(r['v'][1:]) >= dmin]
    a, b, r2, rms, maxrel = fit(pts)
    print('fit depth>=%2d: ms/step = %6.2f + %6.3f*depth   R2=%.5f  rms=%.2f ms  max|res|/t=%.1f%%  (n=%d)'
          % (dmin, a, b, r2, rms, 100 * maxrel, len(pts)))
pts = [(int(r['v'][1:]), r['ms']) for r in uni if int(r['v'][1:]) >= 3]
a3, b3, r2, rms, maxrel = fit(pts)
# cost per force-pass kernel call from the slope: each level adds 3 calls per EV corner
corners3 = sum(r['corners'] for r in uni if int(r['v'][1:]) >= 3) / len(pts)
print('mean EV corners (depth>=3 runs) = %.0f -> %.0f extra force-kernel calls per level per step -> %.0f ns per call (upper bound: area pass and flip trials add calls)'
      % (corners3, 3 * corners3, b3 * 1e6 / (3 * corners3)))
ptsK = [(r['K'], r['ms']) for r in uni if int(r['v'][1:]) >= 3]
aK, bK, r2K, rmsK, maxrelK = fit(ptsK)
print('fit vs K (depth>=3): ms/step = %.2f + %.1f ns/call * K   R2=%.5f  max|res|/t=%.1f%%' % (aK, bK * 1e6, r2K, 100 * maxrelK))
for r in rows:
    if r['v'] in ('default', 'legacy'):
        print('   %-7s predicted from K: %.2f ms, measured %.2f ms (%.1f%%)' % (r['v'], aK + bK * r['K'], r['ms'], 100 * (aK + bK * r['K'] - r['ms']) / r['ms']))
tf = [r for r in rows if r['v'] == 'frozen'][0]['ms']
print('intercept a=%.2f ms/step (depth->0), frozen=%.2f ms/step (regular lattice, no flips), slope b=%.3f ms/step/level' % (a3, tf, b3))
print('%-8s %8s %12s %12s %9s %10s' % ('variant', 'ms/step', 'vs depth->0', 'vs frozen', 'x frozen', 'x depth->0'))
for r in rows:
    if r['v'] == 'frozen':
        continue
    tt = r['ms']
    print('%-8s %8.2f %11.1f%% %11.1f%% %9.2f %10.2f' % (r['v'], tt, 100 * (tt - a3) / tt, 100 * (tt - tf) / tt, tt / tf, tt / a3))
