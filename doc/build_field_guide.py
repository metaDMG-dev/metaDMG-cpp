#!/usr/bin/env python3
"""Build offline field-guide.html from README.md + local CSS/JS/figures.
Requires Python 3 and Pandoc on PATH. Run from any directory.
Optional --repo points to a documentation staging tree containing README.md.
The guide uses no CDN, analytics, remote fonts, or runtime dependencies.
"""
from pathlib import Path
import argparse,base64,html,re,subprocess
parser=argparse.ArgumentParser();parser.add_argument('--repo',type=Path);args=parser.parse_args()
doc=Path(__file__).resolve().parent;repo=args.repo.resolve() if args.repo else doc.parent
md=(repo/'README.md').read_text()
def convert(s):
 out=subprocess.run(['pandoc','--from=gfm','--to=html5','--wrap=none','--no-highlight'],input=s,text=True,capture_output=True,check=True).stdout
 out=re.sub(r'<table>', '<div class="table-wrap"><table>',out).replace('</table>','</table></div>')
 def embed(m):
  rel=m.group(1);path=repo/rel
  if rel.startswith('doc/figures/') and path.exists():return 'src="data:image/svg+xml;base64,'+base64.b64encode(path.read_bytes()).decode()+'"'
  return m.group(0)
 out=re.sub(r'src="([^"]+)"',embed,out)
 out=out.replace('href="doc/','href="').replace('href="README.md"','href="../README.md"')
 return out
model='''<div class="lab" aria-label="Interactive damage model illustration">
<div class="plate"><span class="eyebrow">Plate 03 · Try the model</span>
<svg viewBox="0 0 520 330" role="img" aria-labelledby="plot-title plot-description"><title id="plot-title">Terminal damage mean curve</title><desc id="plot-description">Illustrative curve, no measured data.</desc>
<g font-family="sans-serif" font-size="12" fill="#526057">
<text x="18" y="22">Expected damage frequency</text>
<path d="M62 48H462 M62 121.33H462 M62 194.67H462 M62 268H462" fill="none" stroke="#dce2d8"/>
<text x="25" y="52">60%</text><text x="25" y="125">40%</text><text x="25" y="199">20%</text><text x="32" y="272">0%</text>
<text x="58" y="290">0</text><text x="170" y="290">4</text><text x="284" y="290">8</text><text x="394" y="290">12</text><text x="454" y="290">14</text>
<text x="260" y="318" text-anchor="middle">Distance inward (terminal base = 0)</text></g>
<path d="M62 48V268H462" fill="none" stroke="#526057"/>
<line id="background-line" x1="62" x2="462" y1="260.67" y2="260.67" stroke="#285e45" stroke-dasharray="6 5" stroke-width="2"/>
<path id="model-line" d="M62 187.33 L90.57 209.33 L119.14 224.73 L147.71 235.51 L176.29 243.06 L204.86 248.34 L233.43 252.04 L262 254.63 L290.57 256.44 L319.14 257.71 L347.71 258.60 L376.29 259.22 L404.86 259.65 L433.43 259.96 L462 260.17" fill="none" stroke="#ba481f" stroke-width="3"/>
</svg>
<p class="caption">Orange: expected mean D(x). Dashed green: background c. Illustrative parameters only; no observations or confidence band.</p>
<div class="controls">
<div class="control"><label for="amplitude">A · amplitude <output id="amplitude-value" for="amplitude">20%</output></label><input id="amplitude" type="range" min="0" max="0.4" step="0.01" value="0.20"></div>
<div class="control"><label for="decay">q · decay <output id="decay-value" for="decay">0.30</output></label><input id="decay" type="range" min="0.01" max="0.95" step="0.01" value="0.30"></div>
<div class="control"><label for="background">c · background <output id="background-value" for="background">2.0%</output></label><input id="background" type="range" min="0" max="0.1" step="0.005" value="0.02"></div>
<div class="control"><label for="concentration">phi · concentration <output id="concentration-value" for="concentration">100</output></label><input id="concentration" type="range" min="1" max="1000" step="1" value="100"></div>
</div><div class="readout" aria-live="polite">First base, A + c: <strong id="terminal">22.0%</strong><br>At distance 14: <strong id="interior">2.14%</strong><br>Illustrative frequency SD at the first base, N = 100: <strong id="dispersion">5.83 percentage points</strong></div>
<p class="caption">SD = √[p(1−p)(N+phi)/(N(1+phi))], p = A+c. This is beta-binomial count variation at a fixed mean, not the fitted uncertainty of A or a Zfit calculation. Slider limits are demonstration ranges, not recommended filters.</p>
<button class="reset" id="reset-model" type="button">Reset illustration</button><noscript><p>Enable JavaScript to move the controls. The initial curve and all instructions remain readable.</p></noscript></div>
<div class="marks">
<div class="mark"><span class="kind">01 · Mean model</span><h3>D(x) = A(1 − q)<sup>x</sup> + c</h3><p>Start at the read end, x = 0, and move towards its interior. The first-base prediction is A + c.</p></div>
<div class="mark"><span class="kind">02 · Terminal excess</span><h3>A raises the terminal signal</h3><p>Amplitude measures damage above the background component. Moving A changes the height of the curve near the read end.</p></div>
<div class="mark"><span class="kind">03 · Shape</span><h3>q sets how quickly it falls</h3><p>A larger decay rate concentrates the signal at the terminus. A smaller value extends it further into the read.</p></div>
<div class="mark"><span class="kind">04 · Baseline</span><h3>c shifts the background</h3><p>The background can include several mismatch sources. It is not a direct measure of sequencing error alone.</p></div>
<div class="mark"><span class="kind">05 · Variation around the mean</span><h3>phi changes dispersion</h3><p>Move phi: the mean curve stays fixed, but the illustrative standard deviation changes. Larger values approach binomial variation; phi alone does not diagnose reference quality.</p></div>
</div></div>'''
terms=[
('MD:Z','Input tag','Reference mismatch description used to reconstruct aligned reference bases. Required for damage counting.'),
('NM:i / AS:i','Input tags','Keep NM for edit-distance/ANI filtering. AS is required when --best_as 1 is enabled.'),
('--run_mode','Counting choice','0 pools alignments; 1 groups by BAM reference; 2 groups references by mapped taxid without LCA.'),
('--weight_type','LCA option','Source default is 0 (full alignment weight). Specify 1 for fractional alignment contributions to damage counts.'),
('--lca_rank','LCA option','Controls rank acceptance for damage/statistics. An ambiguous genus assignment is not resolved to species.'),
('--lib ds / ss','Library model','ds: 5′ C→T and 3′ G→A. ss: C→T at both ends. Match the library preparation.'),
('--nbootstrap','Fit option','1: default fit without bootstrap summaries. Values >1 enable resampling. 0 is rejected in this checkout.'),
('A / q / c / phi','Fit outputs','Amplitude above background / positional decay / background mismatch component / beta-binomial concentration.'),
('fwK{i} / bwK{i}','Observed counts','Numerator counts for the informative substitutions at each end; interpreted with the selected library model.'),
('fwN{i} / bwN{i}','Reference counts','Denominators for substitution frequencies. Support varies by position and group; zero support needs care.'),
('fwf{i} / bwf{i}','Observed frequencies','K/N at terminal index i. Index 0 denotes the first base. Available with --showfits 2.'),
('fwdx{i} / bwdx{i}','Fitted frequencies','Model predictions at index i, rather than raw observed frequencies. Available with --showfits 1 or 2.'),
('sigmaD_old / sigmaD_new','Fit uncertainty','Legacy and updated uncertainty estimates. Pair each with the corresponding Zfit column.'),
('Zfit_old / Zfit_new','Damage scores','Not probabilities that a taxon is ancient, correctly mapped, or correctly classified.'),
('A_b / q_b / c_b / phi_b','Bootstrap summaries','In this C++ checkout these summarize bootstrap estimates; they are not Bayesian MAP outputs.'),
('*_CI_l / *_CI_h','Bootstrap intervals','Lower and upper confidence bounds under the chosen resampling and CI method; do not label them Bayesian credible intervals.'),
('nalign / nreads','Count fields','Keep alignment support distinct from read-group counts and weighted mismatch counts. Ancestor and descendant rows overlap.'),
('.bdamage.gz / .dfit.gz','Output files','Binary substitution matrix (export with print) / compressed tab-separated fitted results (open with gzip -dc).')]
rows=''.join('<tr>'+''.join('<td>'+html.escape(v)+'</td>' for v in row)+'</tr>' for row in terms)
glossary='''<section class="chapter" id="parameter-reference"><h2>7. Parameter field notes</h2><div class="body"><p>Look up the command options and output columns used in this guide. Definitions refer to this C++ checkout; downstream pipelines can rename fields.</p><div class="lookup"><label for="term-search">Find a parameter</label><input type="search" id="term-search" placeholder="Try bootstrap, weight, or frequency"></div><p id="term-status" role="status">18 entries</p><div class="table-wrap"><table id="parameter-table"><thead><tr><th scope="col">Symbol</th><th scope="col">Role</th><th scope="col">What it tells you</th></tr></thead><tbody>'''+rows+'''</tbody></table></div><p>Detailed options: <a href="commands.md">advanced command reference</a>. This dictionary intentionally distinguishes this checkout's bootstrap fields from the terminology in the Holi reference guide.</p></div></section>'''
parts=re.split(r'^## (.+)$',md,flags=re.M)
intro=parts[0]
# Exclude GitHub navigation/badge and the link to this generated companion.
intro=re.sub(r'^# metaDMG-cpp\n','',intro)
intro='\n'.join(l for l in intro.splitlines() if not l.startswith(('[![tests]','[Install]','**Visual field guide:**')))
chapters=[];ids=[]
for i in range(1,len(parts),2):
 title,body=parts[i:i+2];slug=re.sub(r'[^a-z0-9 -]','',title.lower()).replace(' ','-');ids.append(slug)
 body=convert(body)
 if title.startswith('5.'):
  # Replace the first static model figure/caption with the interactive plate.
  body=re.sub(r'<(?:figure|p)>\s*<img.*?</(?:figure|p)>',lambda m:model,body,count=1,flags=re.S)
  body=re.sub(r'<p><em>Illustration only:.*?</em></p>','',body,flags=re.S)
 chapters.append(f'<section class="chapter" id="{slug}"><h2>{html.escape(title)}</h2><div class="body">{body}</div></section>')
nav=''.join(f'<a href="#{id}">{i+1:02} {label}</a>' for i,(id,label) in enumerate(zip(ids,['Install','Test','Prepare','Run','Interpret','Troubleshoot'])))+'<a href="#parameter-reference">07 Field notes</a>'
css=(doc/'field-guide.css').read_text();js=(doc/'field-guide.js').read_text()
page='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="Install, test, run, and interpret metaDMG-cpp. Offline visual field guide with an interactive damage model."><title>metaDMG · Field guide</title><style>'''+css+'''</style></head><body><a class="skip" href="#main">Skip to guide</a><header><div class="mast"><a class="brand" href="#main">metaDMG <small>A field guide to damage analysis</small></a><div class="header-links"><a href="../README.md">README</a><a href="commands.md">Command reference</a></div></div><nav aria-label="Guide chapters">'''+nav+'''</nav></header><main id="main"><div class="intro"><div class="eyebrow">Mapped reads → damage counts → fitted evidence</div><h1>Read the damage.<br>Understand the evidence.</h1><p class="lead">A practical route through installation, your first successful run, and the meaning of the results.</p><div class="cards"><a class="card" href="#1-install"><span>START HERE</span><b>New to metaDMG?</b><span>Install the tools, then run a bundled example.</span></a><a class="card" href="#4-choose-and-run-an-analysis"><span>BRING YOUR DATA</span><b>Choose an analysis</b><span>Whole sample, individual references, or taxonomic groups.</span></a><a class="card" href="#5-understand-the-results"><span>READ THE RESULTS</span><b>Explore the model</b><span>See what amplitude, decay, background, and dispersion do.</span></a></div>'''+convert(intro)+'''</div>'''+''.join(chapters)+glossary+'''<footer class="footer"><p>Documentation for the dev2026 checkout · Visual structure inspired by the supplied Holi Field Guide. The README is the source for installation and workflow instructions; the parameter notes and model explorer are maintained alongside the generator.</p><p>All model illustrations are synthetic. This guide loads offline; links to repository Markdown/PDF references require the accompanying repository files. No external scripts or fonts are loaded.</p></footer></main><script>'''+js+'''</script></body></html>'''
assert 'id="model-line"' in page
(doc/'field-guide.html').write_text(page)
print('Generated',doc/'field-guide.html',len(page.encode()),'bytes')
