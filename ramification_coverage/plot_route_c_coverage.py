"""Regenerates route_c_coverage.pdf/.png from route_c_coverage_data.csv.

Columns of the CSV, for 4 <= d <= 2000:
  eps_alpha_fundamental  1 if O_F^x = <-1, eps_alpha>            (computed, place_ramification.py)
  in_RD_family           1 if d = 2m+1, m even, m^2-1 squarefree  (Thm. "unconditional family")
  cyclic_2Sylow          1 if the 2-Sylow of Cl(d' inf_-) is cyclic
                         (the 24 dimensions quoted from Sec. "class-field form"; NOT recomputed here)
Ramification at inf_- (condition (C_d)) holds for every d in the range (Table "ram"), so the
set of all d is drawn as the reference curve.

Panel (a): counting functions #{d' <= d} on log-log axes. A set of positive density is a line
of slope 1, parallel to the reference; a set of density zero bends away from it.
Panel (b): the same sets as a share of all dimensions <= d, with the predicted density of the
Richaud-Degert family, (1/4) prod_{p>2} (1 - 2/p^2) = 0.16132, as a dotted line.
"""
import csv
import matplotlib.pyplot as plt

rows = list(csv.DictReader(open('route_c_coverage_data.csv')))
d = [int(r['d']) for r in rows]
n_all = [x - 3 for x in d]                                    # all d' in [4, d]
eps = [int(r['cum_eps_alpha_fundamental']) for r in rows]
rd = [int(r['cum_RD_family']) for r in rows]
syl = [int(r['cum_cyclic_2Sylow']) for r in rows]
RD_DENS = 0.16132

INK, MUTED, GRID = '#1f1f1f', '#6b6b6b', '#e6e6e6'
BLUE, ORANGE, AQUA, REF = '#2a78d6', '#eb6834', '#1baf7a', '#9a9a9a'   # validated, first 3 slots
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['cmr10'], 'axes.formatter.use_mathtext': True, 'axes.unicode_minus': False, 'mathtext.fontset': 'cm', 'font.size': 10,
    'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
})
series = [  # (values, color, style, width, label)
    (n_all, REF, (0, (1, 1.6)), 1.4, 'all $d$:  $|\\mathrm{Cl}(d^{\\prime}\\infty_-)|=2\\,|\\mathrm{Cl}(d^{\\prime})|$  (verified)'),
    (eps, BLUE, '-', 1.8, '$\\mathcal{O}_F^{\\times}=\\langle -1,\\varepsilon_\\alpha\\rangle$,  $\\varepsilon_\\alpha=\\frac{1}{2}(d-1+\\sqrt{\\Delta})$  (Thm. VI.37)'),
    (rd, AQUA, (0, (5, 1.6)), 1.8, '$d=2m+1$,  $m$ even,  $m^2-1$ squarefree  (Thm. VI.38, unconditional)'),
    (syl, ORANGE, (0, (3, 1.2, 1, 1.2)), 1.8, '$\\mathrm{Cl}(d^{\\prime}\\infty_-)[2^\\infty]$ cyclic  (Prop. VI.30, unique involution)'),
]

fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.6, 4.2), gridspec_kw={'wspace': 0.55})

# (a) counting functions, log-log
for v, c, ls, lw, lab in series:
    ax.plot(d, v, color=c, ls=ls, lw=lw, solid_capstyle='round', label=lab)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlim(4, 2000); ax.set_ylim(0.8, 3000)
ax.set_xticks([5, 10, 50, 100, 500, 2000]); ax.set_xticklabels(['5', '10', '50', '100', '500', '2000'])
ax.set_yticks([1, 10, 100, 1000]); ax.set_yticklabels(['1', '10', '100', '1000'])
ax.set_xlabel('dimension $d$'); ax.set_ylabel('number of dimensions $\\leq d$')
ax.grid(True, which='major', color=GRID, lw=0.6); ax.set_axisbelow(True)
for (v, c, ls, lw, lab), dy in zip(series, (5, -5, 0, 0)):   # end values, in text ink
    ax.annotate(f'{v[-1]}', (2000, v[-1]), xytext=(4, dy), textcoords='offset points',
                va='center', ha='left', fontsize=9, color=INK)
ax.set_title('(a)  counting functions', loc='left', fontsize=10, color=INK)

# (b) share of all dimensions <= d
lo = 49   # start at d = 50 so that small-d noise does not dominate
bx.axhline(100 * RD_DENS, color=AQUA, lw=0.8, ls=':', zorder=1)
bx.annotate('$\\frac{1}{4}\\,\\Pi_{p>2}(1-2p^{-2})\\approx 16.13\\%$', (1950, 100 * RD_DENS), xytext=(0, 6), textcoords='offset points',
            ha='right', va='bottom', fontsize=9, color=INK)
for (v, c, ls, lw, lab), dy in zip(series, (5, -5, 0, 0)):
    share = [100 * v[i] / n_all[i] for i in range(lo, len(d))]
    bx.plot(d[lo:], share, color=c, ls=ls, lw=lw, solid_capstyle='round')
    bx.annotate(f'{share[-1]:.1f}%', (2000, share[-1]), xytext=(4, dy), textcoords='offset points',
                va='center', ha='left', fontsize=9, color=INK)
bx.set_xlim(50, 2000); bx.set_ylim(0, 105)
bx.set_xticks([50, 500, 1000, 1500, 2000])
bx.set_yticks([0, 20, 40, 60, 80, 100]); bx.set_yticklabels(['0', '20', '40', '60', '80', '100%'])
bx.set_xlabel('dimension $d$'); bx.set_ylabel('share of dimensions $\\leq d$')
bx.grid(True, axis='y', color=GRID, lw=0.6); bx.set_axisbelow(True)
bx.set_title('(b)  share of all dimensions', loc='left', fontsize=10, color=INK)

for a in (ax, bx):
    for sp in a.spines.values():
        sp.set_visible(True); sp.set_color(MUTED); sp.set_linewidth(0.8)
    a.tick_params(which='both', top=True, right=True, direction='in')

fig.legend(*ax.get_legend_handles_labels(), loc='upper center', ncol=1, frameon=False,
           fontsize=9, handlelength=2.8, columnspacing=1.4, bbox_to_anchor=(0.5, 1.0))
fig.subplots_adjust(left=0.11, right=0.89, bottom=0.14, top=0.71)
fig.savefig('route_c_coverage.pdf')
fig.savefig('route_c_coverage.png', dpi=300)
print('endpoints at d=2000:', eps[-1], rd[-1], syl[-1], 'of', n_all[-1])
