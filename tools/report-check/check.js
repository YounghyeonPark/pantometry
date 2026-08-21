'use strict';
// What the report must do when it is actually run.
//
//     node tools/report-check/check.js out1.html out2.html ...
//
// Every assertion here is about something that came out of the viewer, not about the HTML that
// went in. A string check cannot tell "the heatmap drew 2,834 cells" from "the heatmap threw on
// line 3 and the page is blank", and both of those have shipped from this repository before in
// other renderers.

const { load } = require('./harness.js');

let failures = 0;
let checks = 0;
let current = '';

function ok(cond, what, detail) {
  checks++;
  if (cond) return true;
  failures++;
  console.log('  FAIL  ' + what + (detail ? '\n          ' + detail : ''));
  return false;
}
function note(s) { console.log('        ' + s); }

// --- CIELAB, so the page's own colour table can be checked where it will be used ----------------
function lightness(r, g, b) {
  const lin = (c) => {
    c /= 255;
    return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
  };
  const [R, G, B] = [lin(r), lin(g), lin(b)];
  const y = (0.2126729 * R + 0.7151522 * G + 0.0721750 * B) / 1.0000001;
  return y > 216 / 24389 ? 116 * Math.cbrt(y) - 16 : 24389 / 27 * y;
}

function checkReport(path) {
  current = path;
  console.log('\n=== ' + path);
  let r;
  try {
    r = load(path);
  } catch (e) {
    failures++;
    console.log('  FAIL  the viewer did not run: ' + e.message);
    return;
  }

  // --- the colour table that reached the page -------------------------------------------------
  const maps = JSON.parse(r.byId.maps.textContent);
  ok(maps.seq && maps.seq.length === 256 * 6, 'the sequential table is 256 entries',
     'got ' + (maps.seq || '').length + ' hex characters');
  ok(maps.div && maps.div.length === 256 * 6, 'the diverging table is 256 entries');
  let back = 0;
  let prev = -1;
  for (let i = 0; i < 256; i++) {
    const L = lightness(parseInt(maps.seq.substr(6 * i, 2), 16),
                        parseInt(maps.seq.substr(6 * i + 2, 2), 16),
                        parseInt(maps.seq.substr(6 * i + 4, 2), 16));
    if (prev >= 0 && L < prev - 0.1) back++;
    prev = L;
  }
  ok(back === 0, 'no step of the sequential scale is darker than the one below it',
     back + ' of 255 steps run backwards, which is what the four-stop ramp did');

  // --- every view drew something ----------------------------------------------------------------
  const slots = r.canvases.map((c) => c.dataset.slot);
  ok(slots.length > 0, 'the page has at least one view');
  for (const c of r.canvases) {
    const slot = c.dataset.slot;
    const ctx = c._ctx;
    if (!ok(ctx, slot + ': the view asked for a drawing context')) continue;
    const drew = ctx._calls.fillRect + ctx._calls.drawImage + ctx._calls.stroke + ctx._calls.fill;
    ok(drew > 2, slot + ': the view drew something',
       'fillRect ' + ctx._calls.fillRect + ', drawImage ' + ctx._calls.drawImage
       + ', stroke ' + ctx._calls.stroke);
    // The backing store must follow the element, not a number baked into the markup: the report
    // ships width="1400" and the element is 1148 CSS pixels wide.
    ok(c.width !== 1400 || c.clientWidth === 1400,
       slot + ': the canvas was sized from the element and not from the markup',
       'canvas.width is still ' + c.width + ' for a ' + c.clientWidth + ' px element');
    ok(r.caption(slot).length > 10, slot + ': the view wrote a caption',
       JSON.stringify(r.caption(slot)));
  }

  // --- numbers on the page ------------------------------------------------------------------------
  for (const c of r.canvases) {
    const slot = c.dataset.slot;
    const kind = c.dataset.kind;
    const text = (c._ctx ? c._ctx._seen.text : []).join(' | ');
    if (kind === 'heatmap' || kind === 'slices' || kind === 'profile') {
      // Ticks are bare numbers and the unit is named once on the axis — "x / mm" — so this looks
      // for the unit token wherever it is, not for a unit glued to every tick. An axis whose
      // extent collapsed to nothing says "cell index" instead, and that is also an answer.
      ok(/(^|[\s/])(nm|µm|mm|m)\b/.test(text) || /cell index/.test(text),
         slot + ': a spatial axis is labelled in real units',
         'the labels drawn were ' + JSON.stringify((c._ctx ? c._ctx._seen.text : []).slice(0, 14)));
    }
    if (kind === 'heatmap' || kind === 'slices' || kind === 'volume'
        || kind === 'scene' || kind === 'layout') {
      // The colour bar draws at least three tick labels plus the unit.
      const labels = (c._ctx ? c._ctx._seen.text : []).length;
      ok(labels >= 3, slot + ': the view drew a colour bar with numbers on it',
         labels + ' text runs in total');
    }
    if (kind === 'series') {
      const t = c._ctx ? c._ctx._seen.text : [];
      const numeric = t.filter((s) => /^-?[\d.]/.test(s) || /e[+-]?\d/.test(s));
      ok(numeric.length >= 6, slot + ': the scalar chart has axis numbers',
         'only ' + numeric.length + ' numeric labels: ' + JSON.stringify(t.slice(0, 10)));
    }
  }

  // --- taking a figure away ------------------------------------------------------------------
  // A report ends up in a document. Without this the only way to get a figure out of one was a
  // screenshot at whatever size the window happened to be.
  for (const c of r.canvases) {
    const slot = c.dataset.slot;
    const a = r.savePng(slot);
    if (!ok(a, slot + ': the PNG button produced a download')) continue;
    ok(/^data:image\/png/.test(a.href || ''), slot + ': and it is a data URL, not a fetch',
       'href was ' + JSON.stringify((a.href || '').slice(0, 40)));
    ok(/\.png$/.test(a.download || ''), slot + ': with a .png filename',
       'filename was ' + JSON.stringify(a.download));
  }

  // --- pointing at things ---------------------------------------------------------------------
  for (const c of r.canvases) {
    const slot = c.dataset.slot;
    if (!['heatmap', 'slices', 'profile', 'scene', 'series'].includes(c.dataset.kind)) continue;
    // The element is 1148 CSS pixels wide in the stub; aim at the middle of it.
    const rect = c.getBoundingClientRect();
    let got = '';
    for (const [fx, fy] of [[0.5, 0.5], [0.4, 0.45], [0.6, 0.55], [0.5, 0.35], [0.35, 0.6]]) {
      got = r.hover(slot, rect.width * fx, rect.height * fy);
      if (got) break;
    }
    ok(got.length > 0, slot + ': hovering reads a value back',
       'the readout stayed empty over five points in the middle of the view');
  }

  // --- the animation runs, and how fast --------------------------------------------------------
  // **Before the transport checks, not after.** Pressing a key pauses playback — that is what a
  // reader wants and it is also what made the first version of this measurement vacuous: every
  // report reported a median frame time of 0.0 ms, because by then the loop was paused and
  // `step` was timing an early return. The frame counter is asserted to have moved, so a future
  // change that pauses here fails instead of reporting a very fast nothing.
  const seen = new Set([r.byId.tick.textContent]);
  const times = [];
  for (let i = 0; i < 12; i++) {
    const one = r.step(1);
    if (!one.length) break;
    times.push(one[0]);
    seen.add(r.byId.tick.textContent);
  }
  ok(times.length === 12, 'the loop kept scheduling for twelve frames',
     'it stopped after ' + times.length);
  // Distinct ticks, not "the tick differs from the first one": twelve steps of a twelve-frame run
  // wrap exactly back to where they started, and comparing the ends said the loop had not moved.
  ok(seen.size > 1, 'the loop actually advanced the frame while timing',
     'the tick never changed, so the times below measure an early return');
  if (times.length) {
    const sorted = times.slice().sort((a, b) => a - b);
    const median = sorted[sorted.length >> 1];
    note('median frame ' + median.toFixed(1) + ' ms, worst ' + sorted[sorted.length - 1].toFixed(1)
         + ' ms, over ' + slots.length + ' view(s)');
    // 62 ms is the loop's own interval at 1x. A frame that takes longer than that cannot keep up
    // with its own clock, and this stub does no compositing, so a browser has strictly less room.
    ok(median < 62, 'a frame is drawn inside the interval the loop asks for',
       'median ' + median.toFixed(1) + ' ms against a 62 ms budget');
  }

  // --- dragging a 3D view, and a display that is not 1x ------------------------------------------
  const spin = r.canvases.find((c) => ['scene', 'volume', 'layout'].includes(c.dataset.kind));
  if (spin) {
    const listeners = spin._listeners;
    const before = spin._ctx._calls.fillRect + spin._ctx._calls.drawImage;
    (listeners.pointerdown || []).forEach((f) => f({ clientX: 100, clientY: 100, pointerId: 1 }));
    (listeners.pointermove || []).forEach((f) => f({ clientX: 160, clientY: 130, pointerId: 1 }));
    (listeners.pointerup || []).forEach((f) => f({ pointerId: 1 }));
    const after = spin._ctx._calls.fillRect + spin._ctx._calls.drawImage;
    ok(after > before, spin.dataset.slot + ': dragging redrew the view',
       'the camera moved and nothing was painted');
  }

  // A high-DPI display gets a bigger backing store, not a blurrier one. The report ships
  // `width="1400"`; a viewer that ignores the device pixel ratio hands a 2x screen half the
  // resolution it can show, and draws 13 px labels at about 10.
  const hi = load(path, { dpr: 2 });
  for (const c of hi.canvases) {
    const one = r.canvases.find((q) => q.dataset.slot === c.dataset.slot);
    ok(c.width === one.width * 2, c.dataset.slot + ': a 2x display gets a 2x backing store',
       'it is ' + c.width + ' wide at dpr 2 and ' + one.width + ' at dpr 1');
    ok((c.style.height || '') === (one.style.height || ''),
       c.dataset.slot + ': and the same size on the page',
       'CSS height went from ' + one.style.height + ' to ' + c.style.height);
  }

  // --- the transport ----------------------------------------------------------------------------
  const before = r.byId.tick.textContent;
  r.press('ArrowRight');
  const after = r.byId.tick.textContent;
  ok(before !== after, 'the right arrow advances a frame',
     'the tick said ' + JSON.stringify(before) + ' both times');
  r.press('Home');
  ok(/frame 1 \//.test(r.byId.tick.textContent), 'Home returns to the first frame',
     r.byId.tick.textContent);

}

const paths = process.argv.slice(2);
if (!paths.length) {
  console.error('usage: node check.js <report.html> [...]');
  process.exit(2);
}
paths.forEach(checkReport);

console.log('\n' + checks + ' checks, ' + failures + ' failed');
process.exit(failures ? 1 : 0);
