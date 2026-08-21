'use strict';
// A stub browser, enough of one to run the viewer that `pantometry-view::report` inlines.
//
// Every test of that report has been a test of the HTML as a *string*: it contains a canvas with
// this data-kind, it mentions that unit. None of them has ever run the JavaScript, so a syntax
// error, a renamed element or a view that silently draws nothing would ship and every assertion
// would still pass. This runs it and reports what it drew.
//
// **Main realm, not `vm.runInContext`.** Measured on one identical hot loop: node runs it in
// 133 ms in the main realm and 3964 ms inside a vm context — a 30x penalty the optimiser never
// pays back. The first version of this file used a vm and measured the volume renderer at 118 ms
// a frame; in the main realm the same code is 4. Any performance number taken through a vm is a
// number about node. `new Function` with the globals as parameters stubs just as completely.

const fs = require('fs');

function makeCtx(w, h) {
  const calls = {
    fillRect: 0, strokeRect: 0, arc: 0, stroke: 0, fill: 0,
    fillText: 0, drawImage: 0, putImageData: 0, moveTo: 0, lineTo: 0,
  };
  const seen = { fill: new Set(), stroke: new Set(), text: [] };
  let painted = 0;
  const ctx = {
    canvas: { width: w, height: h },
    set fillStyle(v) { this._fs = v; seen.fill.add(String(v)); },
    get fillStyle() { return this._fs; },
    set strokeStyle(v) { this._ss = v; seen.stroke.add(String(v)); },
    get strokeStyle() { return this._ss; },
    lineWidth: 1, globalAlpha: 1, font: '', textAlign: '', textBaseline: '',
    imageSmoothingEnabled: true, lineJoin: '', lineCap: '', filter: '',
    shadowBlur: 0, shadowColor: '', globalCompositeOperation: '',
    fillRect(x, y, ww, hh) { calls.fillRect++; painted += Math.max(0, ww) * Math.max(0, hh); },
    strokeRect() { calls.strokeRect++; }, clearRect() {},
    beginPath() {}, moveTo() { calls.moveTo++; }, lineTo() { calls.lineTo++; }, closePath() {},
    arc() { calls.arc++; }, ellipse() {}, rect() {}, quadraticCurveTo() {}, bezierCurveTo() {},
    fill() { calls.fill++; }, stroke() { calls.stroke++; },
    save() {}, restore() {}, translate() {}, rotate() {}, scale() {}, setTransform() {}, clip() {},
    setLineDash() {}, getLineDash() { return []; },
    measureText(t) { return { width: String(t).length * 6 }; },
    fillText(t) { calls.fillText++; seen.text.push(String(t)); },
    strokeText() {},
    createImageData(ww, hh) {
      return { width: ww, height: hh, data: new Uint8ClampedArray(ww * hh * 4) };
    },
    putImageData() { calls.putImageData++; },
    drawImage() { calls.drawImage++; },
    createLinearGradient() { return { addColorStop() {} }; },
    createRadialGradient() { return { addColorStop() {} }; },
  };
  ctx._calls = calls;
  ctx._seen = seen;
  ctx._painted = () => painted;
  return ctx;
}

function makeEl(tag) {
  const el = {
    tagName: (tag || 'div').toUpperCase(), dataset: {}, style: {}, children: [],
    width: 300, height: 150, clientWidth: 1148, clientHeight: 0,
    textContent: '', innerHTML: '', value: '', max: '', min: '',
    className: '', attrs: {}, _listeners: {},
    getContext() {
      if (!this._ctx) this._ctx = makeCtx(this.width, this.height);
      // The context is created once and the canvas is resized afterwards, exactly as a browser
      // does it: `canvas.width = n` resets the bitmap without replacing the context object.
      this._ctx.canvas.width = this.width;
      this._ctx.canvas.height = this.height;
      return this._ctx;
    },
    addEventListener(k, f) { (this._listeners[k] = this._listeners[k] || []).push(f); },
    removeEventListener() {}, setPointerCapture() {}, releasePointerCapture() {},
    getBoundingClientRect() {
      // The viewer sets `style.height` after sizing the backing store, exactly as it must for a
      // canvas whose bitmap is `width * devicePixelRatio`. Read it back the way layout would, or
      // the hover mapping is tested against the wrong height on every high-DPI run.
      const h = parseFloat(this.style.height) || this.clientHeight || this.height;
      return { left: 0, top: 0, width: this.clientWidth, height: h,
               right: this.clientWidth, bottom: h };
    },
    appendChild(c) { this.children.push(c); return c; },
    remove() {}, querySelector() { return null; }, querySelectorAll() { return []; },
    setAttribute(k, v) { this.attrs[k] = String(v); }, getAttribute(k) { return this.attrs[k]; },
    focus() {},
    click() {
      (this._listeners.click || []).forEach((f) => f({}));
      if (this.onclick) this.onclick({});
    },
    classList: { add() {}, remove() {}, toggle() {}, contains() { return false; } },
    toDataURL() { return 'data:image/png;base64,'; },
  };
  Object.defineProperty(el, 'innerHTML', {
    get() { return this._html || ''; },
    set(v) {
      this._html = v;
      // Only what the viewer's legend actually writes: `<i></i><span></span>`.
      this.children = (String(v).match(/<(\w+)/g) || []).map((m) => makeEl(m.slice(1)));
    },
  });
  return el;
}

// --- the elements the viewer will ask for, read out of the report's own markup ------------------
function parse(html) {
  const byId = {};
  const canvases = [];
  const idRe = /<(\w+)\b([^>]*?)\bid="([^"]+)"([^>]*)>/g;
  let m;
  while ((m = idRe.exec(html))) {
    const el = makeEl(m[1]);
    el.id = m[3];
    byId[m[3]] = el;
  }
  const canRe = /<canvas\b([^>]*)>/g;
  while ((m = canRe.exec(html))) {
    const a = m[1];
    if (!/class="view"/.test(a)) continue;
    const el = makeEl('canvas');
    const g = (k) => {
      const r = new RegExp(k + '="([^"]*)"').exec(a);
      return r ? r[1] : undefined;
    };
    el.width = Number(g('width') || 300);
    el.height = Number(g('height') || 150);
    el.className = 'view';
    for (const k of ['data-panel', 'data-kind', 'data-slot', 'data-aspect']) {
      const v = g(k);
      if (v !== undefined) el.dataset[k.replace('data-', '')] = v;
    }
    canvases.push(el);
  }
  const grab = (id) => {
    const r = new RegExp('<script id="' + id + '"[^>]*>([\\s\\S]*?)</script>').exec(html);
    const el = makeEl('script');
    el.id = id;
    el.textContent = r ? r[1] : '{}';
    byId[id] = el;
    return el;
  };
  grab('run');
  grab('maps');
  return { byId, canvases };
}

/** Load a report, run its viewer, and hand back everything it touched. */
function load(path, opts) {
  opts = opts || {};
  const html = fs.readFileSync(path, 'utf8');
  const { byId, canvases } = parse(html);
  const scripts = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].map((x) => x[1]);
  const code = scripts[scripts.length - 1];
  if (!code || code.trim().length < 200) throw new Error('no viewer script in ' + path);

  const rafQueue = [];
  const keys = [];
  const requestAnimationFrame = (f) => rafQueue.push(f);
  const created = [];
  const document = {
    getElementById: (id) => byId[id] || null,
    querySelectorAll: (sel) => (sel.includes('canvas') ? canvases : []),
    querySelector: (sel) => (sel.includes('canvas') ? canvases[0] : null),
    createElement: (t) => { const e = makeEl(t); created.push(e); return e; },
    addEventListener() {},
    body: makeEl('body'),
    documentElement: makeEl('html'),
  };
  const win = {
    matchMedia: () => ({ matches: !!opts.reducedMotion, addEventListener() {}, addListener() {} }),
    addEventListener(k, f) { if (k === 'keydown') keys.push(f); },
    devicePixelRatio: opts.dpr || 1,
    requestAnimationFrame,
    innerWidth: 1200, innerHeight: 900, location: { href: '' },
    getComputedStyle: () => ({ getPropertyValue: () => '' }),
    document,
    ResizeObserver: opts.noResizeObserver ? undefined : Observer,
  };
  function Observer(cb) { this.observe = () => {}; this.disconnect = () => {}; this._cb = cb; }

  const fn = new Function(
    'document', 'window', 'requestAnimationFrame', 'cancelAnimationFrame',
    'setTimeout', 'clearTimeout', 'ResizeObserver', code,
  );
  fn(document, win, requestAnimationFrame, () => {}, () => 0, () => {}, win.ResizeObserver);

  let now = 0;
  const api = {
    byId, canvases, rafQueue, html, window: win, path,
    /** Advance the animation loop `n` times, returning the milliseconds each frame took.
     *
     * The clock **persists across calls**. It did not, and the second `step(1)` handed the loop
     * the same timestamp as the first, so `now - last` was zero, the loop returned early and
     * eleven of twelve measurements timed nothing. A monotonic clock is not an implementation
     * detail of a harness; it is the thing being simulated. */
    step(n) {
      const out = [];
      for (let i = 0; i < n; i++) {
        const f = rafQueue.shift();
        if (!f) break;
        now += 200;
        const t0 = process.hrtime.bigint();
        f(now);
        out.push(Number(process.hrtime.bigint() - t0) / 1e6);
      }
      return out;
    },
    /** Press a key, the way a reader would. */
    press(key) { keys.forEach((f) => f({ key, target: { tagName: 'BODY' }, preventDefault() {} })); },
    /** Move the pointer over a view and read back what it says. */
    hover(slot, x, y) {
      const c = canvases.find((q) => q.dataset.slot === slot);
      if (!c) throw new Error('no view ' + slot);
      (c._listeners.pointermove || []).forEach((f) => f({ clientX: x, clientY: y }));
      const out = byId['read-' + slot];
      return out ? out.textContent : '';
    },
    view(slot) { return canvases.find((q) => q.dataset.slot === slot); },
    /** Every element the viewer built for itself — the offscreen volume buffer, a save anchor. */
    created,
    /** Press a card's PNG button and hand back the anchor it made, if it made one. */
    savePng(slot) {
      const b = byId['png-' + slot];
      if (!b) return null;
      const before = created.length;
      b.click();
      return created.slice(before).find((e) => e.tagName === 'A') || null;
    },
    text(slot) {
      const c = api.view(slot);
      return c && c._ctx ? c._ctx._seen.text : [];
    },
    caption(slot) {
      const el = byId['cap-' + slot];
      return el ? el.textContent : '';
    },
  };
  return api;
}

module.exports = { load, makeEl, makeCtx };
