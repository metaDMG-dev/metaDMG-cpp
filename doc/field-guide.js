'use strict';
const $=id=>document.getElementById(id);
function draw(){
 const A=+$('amplitude').value,q=+$('decay').value,c=+$('background').value,phi=+$('concentration').value;
 $('amplitude-value').textContent=(100*A).toFixed(0)+'%';$('decay-value').textContent=q.toFixed(2);$('background-value').textContent=(100*c).toFixed(1)+'%';$('concentration-value').textContent=phi;
 const x=i=>62+i*400/14,y=v=>268-v*220/.6;
 const coords=Array.from({length:15},(_,i)=>[x(i),y(A*Math.pow(1-q,i)+c)]);
 $('model-line').setAttribute('d',coords.map((p,i)=>(i?'L':'M')+p.join(',')).join(' '));
 $('background-line').setAttribute('y1',y(c));$('background-line').setAttribute('y2',y(c));
 $('terminal').textContent=(100*(A+c)).toFixed(1)+'%';$('interior').textContent=(100*(A*Math.pow(1-q,14)+c)).toFixed(2)+'%';
 const p=A+c,N=100;const sd=Math.sqrt(N*p*(1-p)*(N+phi)/(1+phi))/N;
 $('dispersion').textContent=(100*sd).toFixed(2)+' percentage points';
 $('plot-description').textContent=`Illustrative mean damage curve. First-base frequency ${(100*(A+c)).toFixed(1)} percent. Distance 14 frequency ${(100*(A*Math.pow(1-q,14)+c)).toFixed(2)} percent. No observations or fitted confidence intervals.`;
}
for(const id of ['amplitude','decay','background','concentration'])$(id).addEventListener('input',draw);
$('reset-model').addEventListener('click',()=>{for(const [id,v] of Object.entries({amplitude:.20,decay:.30,background:.02,concentration:100}))$(id).value=v;draw();});
draw();
for(const pre of document.querySelectorAll('pre')){
 const code=pre.querySelector('code');if(!code)continue;
 const button=document.createElement('button');button.className='copy';button.textContent='Copy';button.type='button';button.setAttribute('aria-label','Copy command block');
 button.addEventListener('click',async()=>{try{await navigator.clipboard.writeText(code.textContent);button.textContent='Copied';}catch{const range=document.createRange();range.selectNodeContents(code);const selection=window.getSelection();selection.removeAllRanges();selection.addRange(range);button.textContent='Selected — copy manually';}setTimeout(()=>button.textContent='Copy',2500);});pre.append(button);
}
const termRows=[...document.querySelectorAll('#parameter-table tbody tr')];
$('term-search').addEventListener('input',()=>{const q=$('term-search').value.trim().toLowerCase();let count=0;for(const r of termRows){r.hidden=!r.textContent.toLowerCase().includes(q);if(!r.hidden)count++;}$('term-status').textContent=count+' matching entries';});
const nav=[...document.querySelectorAll('nav a')];
if('IntersectionObserver' in window){const observer=new IntersectionObserver(entries=>{for(const entry of entries){if(entry.isIntersecting){for(const a of nav)a.setAttribute('aria-current',a.hash==='#'+entry.target.id?'true':'false');}}},{rootMargin:'-125px 0px -55% 0px'});document.querySelectorAll('.chapter').forEach(s=>observer.observe(s));}
