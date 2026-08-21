//! A self-contained HTML report, with the view for each domain chosen from the data's shape.
//!
//! # The problem this is for
//!
//! A researcher who can state a simulation and cannot draw one. Today this crate offers them a
//! filmstrip — one fixed shape, whatever the physics is — or a table and a data file with the
//! drawing left as an exercise. For a 2D room the filmstrip is right; for a winding it is empty;
//! for a bar it is a colour strip where a profile would say more. None of that is a choice the
//! researcher should have to make, and all of it is a choice the *data* can make.
//!
//! So: one file, opened in a browser, nothing installed. No Python, no plotting library, no
//! toolchain beyond the one that already ran the simulation.
//!
//! # How the view is chosen
//!
//! By shape, not by domain. A new domain gets a sensible picture without this module learning
//! its name, which is the same reason `Domain::as_field` exists.
//!
//! | what the data is | what it becomes |
//! | --- | --- |
//! | scalars over time | a line chart, one axis per unit |
//! | a 1D field | a profile that animates, over a faint ghost of the whole run |
//! | a 2D field | a heatmap that animates, on one colour scale throughout |
//! | a 3D field | a rotatable render **and** every z-slice as a montage — see below |
//! | points in space | a rotatable 3D scene, depth-sorted, that animates |
//! | paths in space | a rotatable 3D layout: rays through a lens train, a trajectory |
//!
//! A volume gets **both**, and that is not indecision. A raycast composites values along each
//! ray, so it shows the shape of a field and a reader cannot get a number back out of it; the
//! montage puts every sample on screen and is quantitative and unreadable as a shape. Offering
//! one would be choosing which half of the question a researcher is allowed to ask.
//!
//! The montage is one slice per tile rather than one slice behind a slider, because a viewer who
//! never touches a slider sees a picture of a solid that is really a picture of one plane.
//!
//! The scale is fixed across the run in every case. A frame that rescales makes a quantity
//! *look* constant while it changes by orders of magnitude, which is the one thing a picture of
//! a simulation must never do.
//!
//! # A picture you can read numbers off
//!
//! Every view here was, until recently, a picture with a colour bar and a cell count. That is
//! enough to see *that* something happened and not enough to say *what*, and the gap showed up in
//! four separate places:
//!
//! - **The axes were in cells.** A room was "61 x 43" and never "6.1 m by 4.3 m", so a hot spot
//!   could be seen and not located. Fields now carry the box they were sampled over — the
//!   `extent_m` of [`PanelData::Field`] — and every spatial axis is labelled in metres.
//! - **The scalar chart had no y-axis at all.** Sixty-four pixels were reserved for labels that
//!   were never drawn, every series was normalised to its own range so all of them filled the
//!   frame, and the only numbers on screen were the current values in the legend. Series are now
//!   grouped by unit, one axis per group, with ticks.
//! - **Nothing could be pointed at.** Hovering now reads back the sample under the cursor: its
//!   index, its position in metres, and its value in the panel's unit.
//! - **The bodies view had no colour bar**, alone among the six. It does.
//!
//! # Why the viewer is a string in this file
//!
//! It is a few hundred lines of JavaScript, inlined into the output. That is not elegant and the
//! alternatives are worse: a separate asset means the report is no longer one file, and a library
//! from a CDN means it does not open on a machine without a network. The whole promise is *open
//! it and it works*.
//!
//! It is not, however, unexecuted. `tools/report-check` runs this viewer against real reports
//! with a stub canvas and asserts on what it drew — see that directory's README for what a
//! measurement taken through a `vm` context is worth, which is nothing.

use crate::data::{compact as num, numbers as nums};
use crate::ramp;
use pantometry_scene::{Frame, PanelData};

/// Build the report for a finished run.
pub fn html(title: &str, frames: &[Frame]) -> String {
    let mut out = String::with_capacity(1 << 16);
    out.push_str("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n");
    out.push_str("<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\n");
    out.push_str(&format!("<title>{}</title>\n", escape(title)));
    out.push_str("<style>\n");
    out.push_str(STYLE);
    out.push_str("</style>\n</head>\n<body>\n<div class=\"wrap\">\n");

    out.push_str(&format!(
        "<header><p class=\"eyebrow\">pantometry-world</p><h1>{}</h1>\
         <p class=\"sub\">{} frames &middot; {:.4} s &middot; conservation held throughout</p>\
         </header>\n",
        escape(title),
        frames.len(),
        frames.last().map_or(0.0, |f| f.time_s)
    ));

    // The transport. `step` either way, a speed, and a scrub — and every one of them also has a
    // key, because a reader comparing two frames presses the same key twice and does not want to
    // hunt a 12 px button to do it.
    out.push_str(
        "<div class=\"bar\">\
         <button id=\"back\" title=\"Previous frame (left arrow)\" aria-label=\"Previous frame\">\
         &#9664;|</button>\
         <button id=\"play\" title=\"Play or pause (space)\">Pause</button>\
         <button id=\"fwd\" title=\"Next frame (right arrow)\" aria-label=\"Next frame\">\
         |&#9654;</button>\
         <input id=\"scrub\" type=\"range\" min=\"0\" value=\"0\" aria-label=\"Frame\">\
         <label class=\"spd\">speed\
         <select id=\"speed\" aria-label=\"Playback speed\">\
         <option value=\"0.25\">0.25x</option><option value=\"0.5\">0.5x</option>\
         <option value=\"1\" selected>1x</option><option value=\"2\">2x</option>\
         <option value=\"4\">4x</option></select></label>\
         <span class=\"tick\" id=\"tick\"></span></div>\n",
    );

    // One card per drawable domain, plus one for the readings, which every domain has.
    //
    // **A volume gets two.** A raycast is what a three-dimensional field looks like and a slice
    // montage is what it *is*, and neither substitutes: the render composites values along a ray,
    // so a reader cannot get a number back out of it, while the montage puts every sample on
    // screen and is unreadable as a shape. Offering only one would be choosing which half of the
    // question a researcher is allowed to ask.
    if let Some(first) = frames.first() {
        for panel in &first.panels {
            let kinds: &[&str] = match &panel.data {
                PanelData::Paths { .. } => &["layout"],
                PanelData::Field { nz, .. } if *nz > 1 => &["volume", "slices"],
                PanelData::Field { ny, .. } if *ny <= 1 => &["profile"],
                PanelData::Field { .. } => &["heatmap"],
                PanelData::Points { .. } => &["scene"],
            };
            for kind in kinds {
                let slot = format!("{}-{}", escape(&panel.name), kind);
                out.push_str(&format!(
                    "<section class=\"card\"><div class=\"head\"><h2>{}</h2>\
                     <span class=\"kind\">{}</span>\
                     <span class=\"read\" id=\"read-{slot}\"></span>\
                     <button class=\"png\" id=\"png-{slot}\" title=\"Save this view as a PNG\">\
                     PNG</button></div>\
                     <canvas class=\"view\" data-panel=\"{}\" data-kind=\"{}\" \
                     data-slot=\"{slot}\" data-aspect=\"{}\" width=\"1400\" height=\"620\" \
                     role=\"img\" aria-label=\"{} of {}\"></canvas>\
                     <p class=\"cap\" id=\"cap-{slot}\"></p></section>\n",
                    escape(&panel.name),
                    match *kind {
                        "profile" => "1D field &middot; profile",
                        "heatmap" => "2D field &middot; heatmap",
                        "layout" => "paths &middot; 3D, drag to rotate",
                        "volume" => "3D field &middot; rendered, drag to rotate",
                        "slices" => "3D field &middot; every z-slice, and the numbers",
                        _ => "bodies &middot; 3D, drag to rotate",
                    },
                    escape(&panel.name),
                    kind,
                    // Wider than tall for a chart, closer to square for anything spatial.
                    match *kind {
                        "profile" => "2.26",
                        _ => "1.9",
                    },
                    kind,
                    escape(&panel.name),
                ));
            }
        }
        if !first.readings.is_empty() {
            out.push_str(
                "<section class=\"card\"><div class=\"head\"><h2>Readings</h2>\
                 <span class=\"kind\">scalars &middot; over time</span>\
                 <span class=\"read\" id=\"read-series\"></span>\
                 <button class=\"png\" id=\"png-series\" title=\"Save this view as a PNG\">\
                 PNG</button></div>\
                 <canvas class=\"view\" data-kind=\"series\" data-slot=\"series\" \
                 data-aspect=\"2.5\" width=\"1400\" height=\"560\" role=\"img\" \
                 aria-label=\"Scalar readings over time\"></canvas>\
                 <div class=\"legend\" id=\"legend\"></div>\
                 <p class=\"cap\" id=\"cap-series\"></p></section>\n",
            );
        }
    }

    out.push_str("</div>\n<script id=\"run\" type=\"application/json\">");
    out.push_str(&json(frames));
    out.push_str("</script>\n<script id=\"maps\" type=\"application/json\">{\"seq\":\"");
    // The colour tables are built by `crate::ramp`, whose tests pin what they are: lightness
    // that never falls, no step that stands out, and two arms that mirror. A second table
    // written in JavaScript would agree on the day it was written.
    out.push_str(&ramp::hex_table(ramp::sequential));
    out.push_str("\",\"div\":\"");
    out.push_str(&ramp::hex_table(ramp::diverging));
    out.push_str("\"}</script>\n<script>\n");
    out.push_str(VIEWER);
    out.push_str("</script>\n</body>\n</html>\n");
    out
}

/// The run as JSON, for the viewer above it.
fn json(frames: &[Frame]) -> String {
    let mut out = String::from("{\"frames\":[");
    for (fi, f) in frames.iter().enumerate() {
        if fi > 0 {
            out.push(',');
        }
        out.push_str(&format!("{{\"t\":{},\"panels\":[", num(f.time_s)));
        for (pi, p) in f.panels.iter().enumerate() {
            if pi > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"name\":{},\"unit\":{},",
                quote(&p.name),
                quote(p.unit)
            ));
            match &p.data {
                PanelData::Field {
                    nx,
                    ny,
                    nz,
                    extent_m,
                    values,
                } => out.push_str(&format!(
                    "\"kind\":\"field\",\"nx\":{nx},\"ny\":{ny},\"nz\":{nz},\"e\":{},\"v\":{}",
                    nums(extent_m),
                    nums(values)
                )),
                PanelData::Paths {
                    vertices,
                    starts,
                    values,
                    bounds,
                } => {
                    let flat: Vec<f64> = vertices.iter().flatten().copied().collect();
                    let heads: Vec<f64> = starts.iter().map(|k| *k as f64).collect();
                    out.push_str(&format!(
                        "\"kind\":\"paths\",\"b\":{},\"s\":{},\"p\":{},\"v\":{}",
                        nums(bounds),
                        nums(&heads),
                        nums(&flat),
                        nums(values)
                    ));
                }
                PanelData::Points {
                    positions,
                    values,
                    bounds,
                    boxed,
                } => {
                    let flat: Vec<f64> = positions.iter().flatten().copied().collect();
                    out.push_str(&format!(
                        "\"kind\":\"points\",\"boxed\":{boxed},\"b\":{},\"p\":{},\"v\":{}",
                        nums(bounds),
                        nums(&flat),
                        nums(values)
                    ));
                }
            }
            out.push('}');
        }
        out.push_str("],\"r\":[");
        for (ri, r) in f.readings.iter().enumerate() {
            if ri > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"d\":{},\"l\":{},\"u\":{},\"v\":{}}}",
                quote(&r.domain),
                quote(&r.label),
                quote(r.unit),
                num(r.value)
            ));
        }
        out.push_str("]}");
    }
    out.push_str("]}");
    out
}

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if c.is_control() => out.push(' '),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

const STYLE: &str = r##"
:root{--bg:#0d1116;--card:#151b24;--rule:#242d3a;--ink:#d7dde6;--dim:#7f8b9a;--key:#5fd0c0;
 --warm:#ffc247;--plot:#0a0e13;--grid:#1c242f}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);line-height:1.5;
 font-family:ui-sans-serif,system-ui,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;-webkit-font-smoothing:antialiased}
.wrap{max-width:1180px;margin:0 auto;padding:30px 18px 60px;display:flex;flex-direction:column;gap:16px}
.eyebrow{font-family:ui-monospace,Consolas,monospace;font-size:11px;letter-spacing:.16em;
 text-transform:uppercase;color:var(--dim);margin:0}
h1{font-size:clamp(21px,3.2vw,29px);font-weight:620;letter-spacing:-.02em;margin:4px 0 2px;text-wrap:balance}
.sub{margin:0;color:var(--dim);font-size:14px}
.bar{position:sticky;top:0;z-index:5;display:flex;align-items:center;gap:10px;padding:10px 12px;
 background:var(--card);border:1px solid var(--rule);border-radius:3px;flex-wrap:wrap}
.bar button{appearance:none;background:var(--bg);color:var(--ink);border:1px solid var(--rule);
 border-radius:3px;font:inherit;font-size:13px;padding:6px 12px;cursor:pointer}
#play{min-width:72px}
.bar button:hover{border-color:var(--key)}
.bar button:focus-visible,select:focus-visible{outline:2px solid var(--warm);outline-offset:1px}
.spd{font-size:12px;color:var(--dim);display:flex;align-items:center;gap:5px}
select{background:var(--bg);color:var(--ink);border:1px solid var(--rule);border-radius:3px;
 font:inherit;font-size:12px;padding:4px 6px}
input[type=range]{flex:1;min-width:120px;accent-color:var(--key)}
.tick{font-family:ui-monospace,Consolas,monospace;font-size:12px;color:var(--dim);
 font-variant-numeric:tabular-nums;min-width:186px;text-align:right}
.card{background:var(--card);border:1px solid var(--rule);border-radius:3px;overflow:hidden}
.head{display:flex;align-items:baseline;gap:12px;padding:12px 14px 10px;border-bottom:1px solid var(--rule)}
h2{margin:0;font-size:15px;font-weight:620;letter-spacing:-.01em}
.kind{font-family:ui-monospace,Consolas,monospace;font-size:11px;color:var(--dim);letter-spacing:.04em}
.read{margin-left:auto;font-family:ui-monospace,Consolas,monospace;font-size:11.5px;
 color:var(--key);font-variant-numeric:tabular-nums;white-space:nowrap}
.png{appearance:none;background:var(--bg);color:var(--dim);border:1px solid var(--rule);
 border-radius:3px;font:inherit;font-family:ui-monospace,Consolas,monospace;font-size:10.5px;
 letter-spacing:.06em;padding:3px 8px;cursor:pointer;flex:none}
.png:hover{color:var(--ink);border-color:var(--key)}
.png:focus-visible{outline:2px solid var(--warm);outline-offset:1px}
@media print{.png{display:none}}
canvas.view{display:block;width:100%;height:auto;background:var(--plot);touch-action:none}
canvas[data-kind=scene],canvas[data-kind=volume],canvas[data-kind=layout]{cursor:grab}
canvas[data-kind=scene]:active,canvas[data-kind=volume]:active,canvas[data-kind=layout]:active{cursor:grabbing}
.legend{display:flex;flex-wrap:wrap;gap:4px 14px;padding:10px 14px 2px}
.legend button{appearance:none;background:none;border:0;padding:2px 0;cursor:pointer;font:inherit;
 font-family:ui-monospace,Consolas,monospace;font-size:11.5px;color:var(--ink);
 display:flex;align-items:center;gap:6px;font-variant-numeric:tabular-nums}
.legend button[aria-pressed=false]{color:var(--dim);text-decoration:line-through}
.legend i{width:11px;height:11px;border-radius:2px;display:inline-block;flex:none}
.legend button[aria-pressed=false] i{opacity:.3}
.cap{margin:0;padding:9px 14px 12px;font-family:ui-monospace,Consolas,monospace;
 font-size:11.5px;color:var(--dim)}
@media (prefers-reduced-motion:reduce){#play{outline:1px dashed var(--dim)}}
@media print{
 body{background:#fff;color:#000}
 .bar{display:none}
 .card{break-inside:avoid;border-color:#bbb}
 .wrap{max-width:none}
}
"##;

const VIEWER: &str = r##"
(function(){
"use strict";
var RUN = JSON.parse(document.getElementById("run").textContent);
var MAPS = JSON.parse(document.getElementById("maps").textContent);
var F = RUN.frames, N = F.length;
var frame = 0, playing = true, last = 0, speed = 1;
var cam = {az:0.7, el:0.4, dist:2.5}, drag = null, dragging = null;

/* ---- colour ---------------------------------------------------------------------------------
   The two tables come from `pantometry-view::ramp`, which is where their properties are pinned:
   lightness that never falls across the sequential scale, no step that stands out, and two arms
   of the diverging scale that mirror. Nothing here invents a colour. */
function unhex(s){
  var t = new Uint8Array(s.length/2);
  for (var i=0;i<t.length;i++) t[i] = parseInt(s.substr(2*i,2),16);
  return t;
}
var SEQ = unhex(MAPS.seq), DIV = unhex(MAPS.div), LEVELS = SEQ.length/3;

/* ---- one scale per panel, across the whole run -----------------------------------------------
   A frame that rescales makes a quantity look constant while it changes by orders of magnitude.

   A range that straddles zero also gets the diverging scale *centred on zero* rather than on the
   middle of the range: for -100 to +300 the neutral colour belongs at the value 0, which is a
   quarter of the way up, and putting it halfway would colour +100 as if it were the datum. */
var range = {};
F[0].panels.forEach(function(p, i){
  var lo = Infinity, hi = -Infinity;
  F.forEach(function(f){ f.panels[i].v.forEach(function(x){
    if(x===null||!isFinite(x))return; if(x<lo)lo=x; if(x>hi)hi=x; }); });
  if(!isFinite(lo)){ lo = 0; hi = 1; }
  var signed = lo < 0 && hi > 0;
  var reach = Math.max(Math.abs(lo), Math.abs(hi)) || 1;
  range[p.name] = {lo:lo, hi:hi, signed:signed, reach:reach,
                   table: signed ? DIV : SEQ,
                   /* where a value sits on its own table, in [0,1] */
                   at: signed ? function(v){ return 0.5 + 0.5*v/reach; }
                              : function(v){ return (v-lo)/((hi-lo)||1); },
                   /* and the inverse, for labelling the colour bar */
                   of: signed ? function(u){ return (2*u-1)*reach; }
                              : function(u){ return lo + u*((hi-lo)||1); }};
});
function css(R, u){
  var i = 3*Math.max(0, Math.min(LEVELS-1, (u*(LEVELS-1))|0));
  return "rgb("+R.table[i]+","+R.table[i+1]+","+R.table[i+2]+")";
}

/* ---- readings, grouped by unit ---------------------------------------------------------------
   Every series used to be normalised to its own range, so six of them filled the frame and two
   temperatures forty kelvin apart drew the same line. Series that share a unit now share an axis
   and can be compared by eye; series that do not are different axes and are labelled as such. */
var series = [], groups = [];
if (F[0].r.length) {
  F[0].r.forEach(function(r, i){
    /* A reading that is missing from a later frame is a hole, not a zero. Substituting 0 was
       both a wrong point and a wrong axis: one absent joule dragged the scale to the origin. */
    var vals = F.map(function(f){ return f.r[i] && f.r[i].l === r.l ? f.r[i].v : NaN; });
    var lo = Infinity, hi = -Infinity, missing = 0;
    vals.forEach(function(v){ if(v!==v){missing++;return;} if(v<lo)lo=v; if(v>hi)hi=v; });
    if(!isFinite(lo)){ lo = 0; hi = 1; }
    series.push({key:r.d+"."+r.l, vals:vals, unit:r.u, missing:missing, on:true, lo:lo, hi:hi});
  });
  series.forEach(function(s){
    var g = groups.filter(function(q){ return q.unit === s.unit; })[0];
    if(!g){ g = {unit:s.unit, members:[]}; groups.push(g); }
    g.members.push(s);
  });
}
function regroup(){
  groups.forEach(function(g){
    var lo = Infinity, hi = -Infinity;
    g.members.forEach(function(s){ if(!s.on) return;
      if(s.lo<lo)lo=s.lo; if(s.hi>hi)hi=s.hi; });
    if(!isFinite(lo)){ lo = 0; hi = 1; }
    if(hi === lo){ hi = lo + (Math.abs(lo) || 1) * 1e-3; }
    g.lo = lo; g.hi = hi;
  });
}
regroup();

/* A field value that is not there arrives as `null`, and `null` is 0 in arithmetic — so
   every reader goes through this and every one of them checks. */
function pv(p, at){ var x = p.v[at]; return x === null ? NaN : x; }

/* ---- numbers on the page ---------------------------------------------------------------------- */
function nice(x){
  if (x === 0) return "0";
  if (!isFinite(x)) return "-";
  var a = Math.abs(x);
  if (a >= 1e5 || a < 1e-3) return x.toExponential(2);
  return x.toFixed(a < 1 ? 4 : a < 100 ? 2 : 1);
}
/* Lengths.

   **One unit per axis, chosen from the span**, not per label from its own magnitude. Choosing per
   label put "200 mm, 400 mm, 600 mm, 800 mm, 1 m, 1.2 m" on one axis of a room, which is a scale
   a reader has to convert in their head halfway up. This workspace draws espresso pucks at 6 cm
   and orbits at 3e11 m, so the unit cannot be fixed either.

   And trailing zeros are only trailing **after a decimal point**. `.replace(/\.?0+$/,"")` turned
   500 into 5, so a 500 µm channel was captioned as a 5 µm one and three slice depths in four were
   wrong. It is the kind of wrong that looks like a plausible number. */
function trim(s){ return s.indexOf(".") < 0 ? s : s.replace(/0+$/,"").replace(/\.$/,""); }
function unitFor(span){
  var a = Math.abs(span);
  if (!isFinite(a) || a === 0) return {k:1, u:"m"};
  if (a >= 1e4) return {k:1, u:"m", exp:true};     /* an orbit is not kilometres */
  if (a >= 1) return {k:1, u:"m"};
  if (a >= 1e-3) return {k:1e3, u:"mm"};
  if (a >= 1e-6) return {k:1e6, u:"µm"};
  return {k:1e9, u:"nm"};
}
function inUnit(v, uf){
  if (v === 0) return "0";
  if (uf.exp) return (v*uf.k).toExponential(2).replace("e+","e");
  var x = v*uf.k, a = Math.abs(x);
  return trim(a >= 100 ? x.toFixed(0) : a >= 10 ? x.toFixed(1)
            : a >= 1 ? x.toFixed(2) : x.toFixed(3));
}
/* A single length with its unit attached, for a caption. */
function len(m){ var uf = unitFor(m); return inUnit(m, uf) + " " + uf.u; }
/* Several lengths that belong together, on the unit the largest of them wants. */
function lens(){
  var vs = [].slice.call(arguments), uf = unitFor(Math.max.apply(null, vs.map(Math.abs)));
  return vs.map(function(v){ return inUnit(v, uf); }).join(" x ") + " " + uf.u;
}
/* Tick positions on a 1-2-5 ladder. Chosen rather than "the range in five equal parts" because a
   tick at 0.3714 is a number nobody reads; the ladder puts them where a person would. */
function ticks(lo, hi, want){
  if(!(hi > lo)) return [lo];
  var raw = (hi-lo)/Math.max(1, want), mag = Math.pow(10, Math.floor(Math.log(raw)/Math.LN10));
  var n = raw/mag, step = mag * (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10);
  /* `k*step`, not a running `t += step`: accumulating a float 40 times puts a tick at
     0.30000000000000004, and the only thing standing between that and the page is that every
     formatter here happens to round. Multiply instead and it cannot arise. */
  var out = [], k0 = Math.ceil(lo/step - 1e-9);
  for(var k = k0; out.length < 40; k++){
    var t = k*step;
    if(t > hi + step*1e-6) break;
    out.push(Math.abs(t) < step*1e-6 ? 0 : t);
  }
  return out;
}

/* ---- canvases ---------------------------------------------------------------------------------
   Sized from the element's own width and the display's pixel ratio, then the context is scaled so
   every draw below works in CSS pixels. The backing store used to be a fixed 1400 wide however
   wide the element was: on a 2x display that is half the resolution the screen can show, and the
   13 px labels landed at about 10. */
var DPR = Math.min(3, window.devicePixelRatio || 1);
var views = [].slice.call(document.querySelectorAll("canvas.view")).map(function(c){
  return {c:c, ctx:c.getContext("2d"), kind:c.dataset.kind, panel:c.dataset.panel,
          aspect: Number(c.dataset.aspect) || 2.2,
          slot:c.dataset.slot || c.dataset.panel || "series", hit:null, w:0, h:0};
});
function resize(v){
  var cw = v.c.clientWidth || v.c.width || 1200;
  var w = Math.max(320, Math.round(cw)), h = Math.round(w / v.aspect);
  if (v.w === w && v.h === h) return false;
  v.w = w; v.h = h;
  v.c.width = Math.round(w * DPR); v.c.height = Math.round(h * DPR);
  v.c.style.height = h + "px";
  return true;
}
function begin(v){
  var x = v.ctx;
  x.setTransform(DPR, 0, 0, DPR, 0, 0);
  x.fillStyle = "#0a0e13"; x.fillRect(0, 0, v.w, v.h);
  x.textBaseline = "alphabetic";
  return x;
}
var MONO = "12px ui-monospace,Consolas,monospace";
var SMALL = "10px ui-monospace,Consolas,monospace";

function panelOf(f, name){
  for (var i=0;i<f.panels.length;i++) if (f.panels[i].name===name) return f.panels[i];
  return null;
}

/* ---- the colour bar, on every view that carries a colour ------------------------------------
   Including the bodies view, which had none: a scene of four spheres said "colour is m/s" in its
   caption and offered nothing to read a speed off. */
function bar(x, v, R, unit){
  var bw = Math.min(260, v.w*0.32), bx = v.w - 16 - bw, by = v.h - 26;
  for(var k=0;k<bw;k++){
    x.fillStyle = css(R, k/(bw-1));
    x.fillRect(bx+k, by, 1, 9);
  }
  x.strokeStyle = "#2c3644"; x.lineWidth = 1;
  x.strokeRect(bx+0.5, by+0.5, bw-1, 9);
  x.fillStyle = "#8b96a5"; x.font = SMALL; x.textAlign = "center";
  var lo = R.of(0), hi = R.of(1);
  ticks(lo, hi, 4).forEach(function(t){
    var u = (t-lo)/((hi-lo)||1);
    if(u < -1e-9 || u > 1+1e-9) return;
    x.fillRect(bx+u*(bw-1), by+9, 1, 3);
    x.fillText(nice(t), bx+u*(bw-1), by+22);
  });
  x.textAlign = "right"; x.fillStyle = "#7f8b9a"; x.font = MONO;
  x.fillText(unit || "(no unit)", bx-10, by+9);
}

/* ---- 1D: a profile, over a ghost of every frame --------------------------------------------- */
function drawProfile(v, f){
  var x = begin(v), w = v.w, h = v.h, p = panelOf(f, v.panel), R = range[v.panel];
  var pad = {l:74, r:20, t:16, b:38};
  var iw = w-pad.l-pad.r, ih = h-pad.t-pad.b;
  var e = p.e, x0 = e[0], x1 = e[3];
  /* An axis in metres, not in cells. A profile across a 400 mm bar was labelled "cell 0" to
     "cell 60", which says how finely it was sampled and nothing about the bar. */
  var flat = !(x1 > x0);
  var sy = function(val){ return pad.t + ih*(1-(val-R.lo)/((R.hi-R.lo)||1)); };
  var sx = function(i){ return pad.l + iw*(p.v.length>1 ? i/(p.v.length-1) : 0.5); };

  x.strokeStyle="#1c242f"; x.lineWidth=1; x.fillStyle="#7f8b9a"; x.font=SMALL;
  ticks(R.lo, R.hi, 5).forEach(function(t){
    var yy = sy(t);
    if(yy < pad.t-1 || yy > pad.t+ih+1) return;
    x.beginPath(); x.moveTo(pad.l,yy); x.lineTo(w-pad.r,yy); x.stroke();
    x.textAlign="right"; x.fillText(nice(t), pad.l-8, yy+3.5);
  });
  var uf = unitFor(x1-x0);
  if(!flat) ticks(x0, x1, 6).forEach(function(t){
    var u = (t-x0)/(x1-x0), xx = pad.l+iw*u;
    x.strokeStyle="#161d27"; x.beginPath(); x.moveTo(xx,pad.t); x.lineTo(xx,pad.t+ih); x.stroke();
    x.textAlign="center"; x.fillText(inUnit(t, uf), xx, h-16);
  });

  /* every frame, faint: the envelope the run covers */
  x.strokeStyle="rgba(125,145,175,0.16)"; x.lineWidth=1.4;
  F.forEach(function(ff){
    var q = panelOf(ff, v.panel); if(!q) return;
    x.beginPath();
    q.v.forEach(function(val,i){ i?x.lineTo(sx(i),sy(val)):x.moveTo(sx(i),sy(val)); });
    x.stroke();
  });
  x.strokeStyle = css(R, 0.86); x.lineWidth=2.6; x.lineJoin="round"; x.beginPath();
  p.v.forEach(function(val,i){ i?x.lineTo(sx(i),sy(val)):x.moveTo(sx(i),sy(val)); });
  x.stroke();

  x.fillStyle="#7f8b9a"; x.font=MONO; x.textAlign="left";
  x.fillText(p.unit, pad.l, pad.t-4);
  x.textAlign="right"; x.fillText(flat ? "cell index" : "x / "+uf.u, w-pad.r, h-16);

  v.hit = function(px, py){
    if(px < pad.l-6 || px > w-pad.r+6 || py < pad.t-6 || py > pad.t+ih+6) return null;
    var i = Math.round((px-pad.l)/iw*(p.v.length-1));
    if(i < 0 || i >= p.v.length) return null;
    var val = pv(p, i);
    var at = flat ? ("cell "+i)
                  : inUnit(x0 + (x1-x0)*(p.v.length>1?i/(p.v.length-1):0.5), uf)+" "+uf.u;
    return {text: at+"   "+(val===val ? nice(val)+" "+p.unit : "no sample"),
            mark: {x:sx(i), y:val===val ? sy(val) : null}};
  };
  cap(v, p.v.length+" samples over "+(flat ? "a collapsed axis" : len(x1-x0))
       +" · scale fixed across the run · faint lines are every frame");
}

/* ---- 2D: a heatmap ------------------------------------------------------------------------ */
function drawHeat(v, f){
  var x = begin(v), w = v.w, h = v.h, p = panelOf(f, v.panel), R = range[v.panel];
  var pad = {l:62, r:16, t:14, b:52};
  var e = p.e, ex = e[3]-e[0], ey = e[4]-e[1];
  var aw = w-pad.l-pad.r, ah = h-pad.t-pad.b;
  /* Square cells unless the extent says otherwise: a 6 m by 3 m room is drawn twice as wide as
     it is tall, because it is. Fitting the grid to the panel instead made every room square. */
  var wantAspect = (ex > 0 && ey > 0) ? (ex/ey) : (p.nx/p.ny);
  var pw = Math.min(aw, ah*wantAspect), ph = pw/wantAspect;
  if(ph > ah){ ph = ah; pw = ph*wantAspect; }
  var ox = pad.l + (aw-pw)/2, oy = pad.t + (ah-ph)/2;
  var cw = pw/p.nx, ch = ph/p.ny;

  for (var j=0;j<p.ny;j++) for (var i=0;i<p.nx;i++){
    var val = pv(p, j*p.nx+i);
    if(!isFinite(val)) continue;
    x.fillStyle = css(R, R.at(val));
    x.fillRect(Math.floor(ox+i*cw), Math.floor(oy+(p.ny-1-j)*ch),
               Math.ceil(cw)+1, Math.ceil(ch)+1);
  }
  x.strokeStyle="#2c3644"; x.lineWidth=1; x.strokeRect(ox+0.5, oy+0.5, pw-1, ph-1);

  /* One unit for both axes, from the larger span: a room 4.4 m by 3.1 m is metres both ways,
     and a die 4 mm by 0.8 mm is millimetres both ways. */
  var uf = unitFor(Math.max(ex, ey));
  x.fillStyle="#7f8b9a"; x.font=SMALL;
  if(ex > 0) ticks(e[0], e[3], 6).forEach(function(t){
    var xx = ox + pw*(t-e[0])/ex;
    x.fillRect(xx, oy+ph, 1, 4);
    x.textAlign="center"; x.fillText(inUnit(t, uf), xx, oy+ph+16);
  });
  if(ey > 0) ticks(e[1], e[4], 4).forEach(function(t){
    var yy = oy + ph*(1-(t-e[1])/ey);
    x.fillRect(ox-4, yy, 4, 1);
    x.textAlign="right"; x.fillText(inUnit(t, uf), ox-8, yy+3.5);
  });
  x.font=MONO; x.textAlign="left"; x.fillStyle="#7f8b9a";
  x.fillText("x, y / "+uf.u, pad.l, pad.t-1);
  bar(x, v, R, p.unit);

  v.hit = function(px, py){
    var i = Math.floor((px-ox)/cw), j = p.ny-1-Math.floor((py-oy)/ch);
    if(i < 0 || j < 0 || i >= p.nx || j >= p.ny) return null;
    var val = pv(p, j*p.nx+i);
    var at = (ex > 0 ? inUnit(e[0]+ex*(p.nx>1?i/(p.nx-1):0.5), uf) : "-")
           + ", " + (ey > 0 ? inUnit(e[1]+ey*(p.ny>1?j/(p.ny-1):0.5), uf) : "-")
           + " " + uf.u;
    return {text:"["+i+","+j+"]  "+at+"   "+(val===val ? nice(val)+" "+p.unit : "no sample")};
  };
  cap(v, p.nx+" x "+p.ny+" samples over "+lens(ex, ey)
       +" · one colour scale across every frame"
       +(R.signed ? " · diverging, neutral at zero" : ""));
}

/* ---- 3D field: every slice, laid out as a montage ------------------------------------------ */
function drawSlices(v, f){
  var x = begin(v), w = v.w, h = v.h, p = panelOf(f, v.panel), R = range[v.panel];
  var pad = 14, gap = 7, aw = w-pad*2, ah = h-pad*2-34;
  var e = p.e, ez = e[5]-e[2], zf = unitFor(ez);
  /* Enough columns that the tiles are as large as possible while all of them fit. Searched
     rather than derived: the tile size depends on the aspect of both the grid and the canvas,
     and nz is small enough that trying every column count is free. */
  var best={s:0,cols:1};
  for(var c=1;c<=p.nz;c++){
    var rows=Math.ceil(p.nz/c);
    var s=Math.min((aw-(c-1)*gap)/(c*p.nx), (ah-(rows-1)*gap-11*rows)/(rows*p.ny));
    if(s>best.s) best={s:s,cols:c};
  }
  var s=best.s, cols=best.cols, rows=Math.ceil(p.nz/cols);
  var tw=cols*p.nx*s+(cols-1)*gap, th=rows*(p.ny*s+11)+(rows-1)*gap;
  var ox=(w-tw)/2, oy=pad+(ah-th)/2+11;
  var tiles = [];
  for(var k=0;k<p.nz;k++){
    var cx=ox+(k%cols)*(p.nx*s+gap), cy=oy+Math.floor(k/cols)*(p.ny*s+11+gap);
    tiles.push({k:k, x:cx, y:cy});
    for(var j=0;j<p.ny;j++) for(var i=0;i<p.nx;i++){
      var val=pv(p, k*p.nx*p.ny+j*p.nx+i);
      if(!isFinite(val)) continue;
      x.fillStyle=css(R, R.at(val));
      x.fillRect(Math.floor(cx+i*s), Math.floor(cy+(p.ny-1-j)*s), Math.ceil(s)+1, Math.ceil(s)+1);
    }
    /* The slice's depth, not just its number: "z5" says which of nine and "22 mm" says where. */
    x.fillStyle="#7f8b9a"; x.font=SMALL; x.textAlign="left";
    x.fillText("z"+k+(ez>0 ? "  "+inUnit(e[2]+ez*(p.nz>1?k/(p.nz-1):0.5), zf)+" "+zf.u : ""),
               cx, cy-3);
  }
  bar(x, v, R, p.unit);

  v.hit = function(px, py){
    for(var t=0;t<tiles.length;t++){
      var i = Math.floor((px-tiles[t].x)/s), j = p.ny-1-Math.floor((py-tiles[t].y)/s);
      if(i < 0 || j < 0 || i >= p.nx || j >= p.ny) continue;
      var val = pv(p, tiles[t].k*p.nx*p.ny+j*p.nx+i);
      return {text:"["+i+","+j+","+tiles[t].k+"]   "
                   +(val===val ? nice(val)+" "+p.unit : "no sample")};
    }
    return null;
  };
  cap(v, p.nx+" x "+p.ny+" x "+p.nz+" over "+lens(e[3]-e[0], e[4]-e[1], ez)
       +" · all "+p.nz+" slices, z increasing · one colour scale across every frame");
}

/* ---- 3D field: a raycast ------------------------------------------------------------------- */
/* Rendered into a small buffer and scaled up. The softness that costs is honest — this view is
   for shape, and the montage beside it carries the numbers.

   Two buffer sizes. Measured in node with a stub canvas, in the main realm: the espresso scene's
   volume takes 30 ms a frame at 220x150x48, and a drag fires a pointermove far more often than
   16 ms apart. So a drag renders at a quarter of the pixels and settles to full quality when the
   pointer stops. */
var VOL = {w:260, h:170}, VOL_DRAFT = {w:120, h:80};

/* Opacity from value, and the choice is not cosmetic.

   A signed field — pressure, zero mean — must be transparent in the middle and opaque at both
   extremes, or a standing wave renders as a solid block. A one-sided field — kelvin, 293 upward —
   must be transparent at the low end, or a block at ambient renders as a solid block for the
   opposite reason. So the transfer function is chosen from the run's own range rather than fixed,
   and `signed` is decided once from whether that range straddles zero. */
function opacity(t, signed){
  var a = signed ? Math.abs(2*t - 1) : t;
  return a*a;                     /* squared, so the quiet bulk clears out of the way */
}

function drawVolume(v, f){
  var x = begin(v), w = v.w, h = v.h, p = panelOf(f, v.panel), R = range[v.panel];
  var buf = (dragging === v) ? VOL_DRAFT : VOL, VW = buf.w, VH = buf.h;
  /* Two samples per cell along the longest axis, capped. Forty-eight steps through a nine-cell
     grid is five times the sampling the data can support and costs exactly that much. */
  var STEPS = Math.max(12, Math.min(64, 2*Math.max(p.nx, p.ny, p.nz)));

  /* The box in world units, longest axis normalised to 1, and **the box's own proportions**:
     a slab 40 mm by 40 mm by 4 mm is drawn as a slab. Normalising the sample counts instead drew
     it by how finely it was sampled, so asking for more slices made the block taller. */
  var e = p.e, sx = e[3]-e[0], sy = e[4]-e[1], sz = e[5]-e[2];
  if(!(sx > 0 && sy > 0 && sz > 0)){ sx = p.nx; sy = p.ny; sz = p.nz; }
  var m = Math.max(sx, sy, sz), ex = [sx/m, sy/m, sz/m];
  var ca=Math.cos(cam.az), sa=Math.sin(cam.az), ce=Math.cos(cam.el), se=Math.sin(cam.el);
  /* Camera basis: right, up, forward. The same angles the bodies view uses, so dragging one
     rotates the other and a reader is never looking at two different orientations. */
  var fwd=[-sa*ce, -se, -ca*ce], right=[ca, 0, -sa], up=[-sa*se, ce, -ca*se];
  var eye=[-fwd[0]*cam.dist, -fwd[1]*cam.dist, -fwd[2]*cam.dist];

  var img=x.createImageData(VW, VH), d=img.data, aspect=VW/VH, k=0, lit=0;
  var tbl = R.table, nxy = p.nx*p.ny;
  for(var py=0; py<VH; py++){
    var sv=(1 - 2*(py+0.5)/VH)*0.75;
    for(var px=0; px<VW; px++){
      var su=(2*(px+0.5)/VW - 1)*0.75*aspect;
      var dir=[fwd[0]+right[0]*su+up[0]*sv, fwd[1]+right[1]*su+up[1]*sv, fwd[2]+right[2]*su+up[2]*sv];
      var ln=Math.hypot(dir[0],dir[1],dir[2]); dir=[dir[0]/ln,dir[1]/ln,dir[2]/ln];

      /* Slab test against the box centred on the origin. */
      var t0=-1e9, t1=1e9, hit=true;
      for(var ax=0; ax<3; ax++){
        var half=ex[ax]/2, o=eye[ax], dd=dir[ax];
        if(Math.abs(dd)<1e-9){ if(o<-half||o>half){hit=false;break;} continue; }
        var a1=(-half-o)/dd, b1=(half-o)/dd, blo=Math.min(a1,b1), bhi=Math.max(a1,b1);
        if(blo>t0)t0=blo; if(bhi<t1)t1=bhi;
      }
      var r=10,g=14,b=19, alpha=0;
      if(hit && t1>t0 && t1>0){
        if(t0<0)t0=0;
        var dt=(t1-t0)/STEPS;
        for(var st=0; st<STEPS && alpha<0.985; st++){
          var tt=t0+dt*(st+0.5);
          /* World point -> grid index, with the box centred on the origin. */
          var gx=((eye[0]+dir[0]*tt)/ex[0]+0.5)*(p.nx-1);
          var gy=((eye[1]+dir[1]*tt)/ex[1]+0.5)*(p.ny-1);
          var gz=((eye[2]+dir[2]*tt)/ex[2]+0.5)*(p.nz-1);
          if(gx<0||gy<0||gz<0||gx>p.nx-1||gy>p.ny-1||gz>p.nz-1) continue;
          var i0=Math.floor(gx), j0=Math.floor(gy), k0=Math.floor(gz);
          var i1=Math.min(i0+1,p.nx-1), j1=Math.min(j0+1,p.ny-1), k1=Math.min(k0+1,p.nz-1);
          var fx=gx-i0, fy=gy-j0, fz=gz-k0;
          var c000=pv(p,k0*nxy+j0*p.nx+i0), c100=pv(p,k0*nxy+j0*p.nx+i1);
          var c010=pv(p,k0*nxy+j1*p.nx+i0), c110=pv(p,k0*nxy+j1*p.nx+i1);
          var c001=pv(p,k1*nxy+j0*p.nx+i0), c101=pv(p,k1*nxy+j0*p.nx+i1);
          var c011=pv(p,k1*nxy+j1*p.nx+i0), c111=pv(p,k1*nxy+j1*p.nx+i1);
          var z0=(c000*(1-fx)+c100*fx)*(1-fy)+(c010*(1-fx)+c110*fx)*fy;
          var z1=(c001*(1-fx)+c101*fx)*(1-fy)+(c011*(1-fx)+c111*fx)*fy;
          var val=z0*(1-fz)+z1*fz;

          if(!isFinite(val)) continue;
          var norm=R.at(val);
          var a=opacity(norm, R.signed)*0.16;
          if(a<=0.0008) continue;
          var ci=3*Math.max(0, Math.min(LEVELS-1, (norm*(LEVELS-1))|0));
          var contrib=a*(1-alpha);
          r+= (tbl[ci]-10)*contrib;
          g+= (tbl[ci+1]-14)*contrib;
          b+= (tbl[ci+2]-19)*contrib;
          alpha+=contrib;
        }
      }
      if(alpha>0.05) lit++;
      d[k++]=r; d[k++]=g; d[k++]=b; d[k++]=255;
    }
  }

  var tmp=document.createElement("canvas"); tmp.width=VW; tmp.height=VH;
  tmp.getContext("2d").putImageData(img,0,0);
  x.imageSmoothingEnabled=true;
  var scale=Math.min(w/VW,(h-34)/VH);
  x.drawImage(tmp,(w-VW*scale)/2,(h-34-VH*scale)/2,VW*scale,VH*scale);
  bar(x, v, R, p.unit);
  v.hit = null;

  /* **Say when the picture is nearly empty.** A localised feature -- one hot cell in a block
     of 729 -- is a small bright dot with everything else transparent, which is correct and
     reads exactly like a broken renderer. Making it look bigger would be making the picture
     lie, so the caption reports how much of the frame carries anything instead.

     Measured rather than guessed: the exponent in `opacity` was tried at 1, 1.5 and 2 against
     three real scenes and moved the occupied fraction from 0.2% to 0.1% on the hot spot. The
     transfer function is not what makes a point source small. */
  var frac = lit/(VW*VH);
  var note = frac < 0.03
    ? " · " + (100*frac).toFixed(1) + "% of the frame is occupied: the feature is that "
      + "small against the whole volume, and the montage below shows where"
    : "";
  cap(v, lens(sx, sy, sz)
       +" · composited along each ray, so this shows shape and not values"
       +" · the montage below carries the numbers · drag to rotate, scroll to zoom"
       + note);
}

/* ---- 3D: paths, depth sorted --------------------------------------------------------------- */
function drawLayout(v, f){
  var x=begin(v), w=v.w, h=v.h, p=panelOf(f,v.panel), R=range[v.panel];
  var b=p.b, s={c:[(b[0]+b[3])/2,(b[1]+b[4])/2,(b[2]+b[5])/2],
                span:Math.max(b[3]-b[0],b[4]-b[1],b[5]-b[2])||1};

  /* Sorted by mean depth and drawn back to front, so a ray in front covers one behind rather
     than whichever happened to be last in the array. */
  var n=p.s.length, order=[];
  for(var k=0;k<n;k++){
    var from=p.s[k], to=(k+1<n?p.s[k+1]:p.p.length/3), dsum=0;
    for(var i=from;i<to;i++) dsum+=project([p.p[3*i],p.p[3*i+1],p.p[3*i+2]],s,w,h).d;
    order.push({k:k, from:from, to:to, d:dsum/Math.max(1,to-from)});
  }
  order.sort(function(a,b2){ return b2.d-a.d; });

  order.forEach(function(o){
    x.beginPath();
    for(var i=o.from;i<o.to;i++){
      var q=project([p.p[3*i],p.p[3*i+1],p.p[3*i+2]],s,w,h);
      i===o.from ? x.moveTo(q.x,q.y) : x.lineTo(q.x,q.y);
    }
    x.strokeStyle=css(R, R.at(p.v[o.k]));
    /* Nearer paths a little heavier, which is the only depth cue a line has. */
    x.lineWidth=Math.max(0.6, 2.6/Math.max(0.4,o.d));
    x.globalAlpha=0.9;
    x.stroke();
  });
  x.globalAlpha=1;
  bar(x, v, R, p.unit);
  v.hit = null;
  cap(v, p.s.length+" paths across "+len(Math.max(b[3]-b[0],b[4]-b[1],b[5]-b[2]))
       +" · colour is "+p.unit+" · drag to rotate, scroll to zoom");
}

/* ---- 3D: bodies, depth sorted -------------------------------------------------------------- */
function project(pt, s, w, h){
  var X=(pt[0]-s.c[0])/s.span, Y=(pt[1]-s.c[1])/s.span, Z=(pt[2]-s.c[2])/s.span;
  var ca=Math.cos(cam.az), sa=Math.sin(cam.az);
  var x1=X*ca-Z*sa, z1=X*sa+Z*ca;
  var ce=Math.cos(cam.el), se=Math.sin(cam.el);
  var y1=Y*ce-z1*se, z2=Y*se+z1*ce;
  var d=z2+cam.dist; if(d<0.05)d=0.05;
  var fo=(Math.min(w,h)*0.60)/d;
  return {x:w/2+x1*fo, y:h/2-y1*fo, d:d, s:fo/Math.min(w,h)};
}
function drawScene(v, f){
  var x=begin(v), w=v.w, h=v.h, p=panelOf(f,v.panel), R=range[v.panel];
  var b=p.b, s={c:[(b[0]+b[3])/2,(b[1]+b[4])/2,(b[2]+b[5])/2],
                span:Math.max(b[3]-b[0],b[4]-b[1],b[5]-b[2])||1};
  if(p.boxed){
    var cs=[];
    for(var i=0;i<8;i++) cs.push(project([i&1?b[3]:b[0], i&2?b[4]:b[1], i&4?b[5]:b[2]], s, w, h));
    var E=[[0,1],[0,2],[0,4],[1,3],[1,5],[2,3],[2,6],[3,7],[4,5],[4,6],[5,7],[6,7]];
    x.strokeStyle="rgba(120,145,180,0.32)"; x.lineWidth=1.3;
    E.forEach(function(e){ x.beginPath(); x.moveTo(cs[e[0]].x,cs[e[0]].y); x.lineTo(cs[e[1]].x,cs[e[1]].y); x.stroke(); });
  }
  var pts=[];
  for(var k=0;k<p.v.length;k++){
    var q=project([p.p[3*k],p.p[3*k+1],p.p[3*k+2]], s, w, h);
    q.val=p.v[k]; q.i=k; pts.push(q);
  }
  pts.sort(function(a,b2){ return b2.d-a.d; });
  var base = p.v.length > 40 ? 6 : 14;
  pts.forEach(function(q){
    q.r=Math.max(1.5, base*q.s*44);
    x.beginPath(); x.arc(q.x,q.y,q.r,0,6.2832);
    x.fillStyle=css(R, R.at(q.val)); x.fill();
  });
  bar(x, v, R, p.unit);

  v.hit = function(px, py){
    /* Front to back, so the body a reader is pointing at is the one they can see. */
    for(var i=pts.length-1;i>=0;i--){
      var q=pts[i];
      if(Math.hypot(px-q.x, py-q.y) <= q.r+3)
        return {text:"body "+q.i+"   "+nice(q.val)+" "+p.unit};
    }
    return null;
  };
  cap(v, p.v.length+" bodies across "+len(s.span)+" · colour is "+p.unit
       +" · drag to rotate, scroll to zoom");
}

/* ---- scalars over time --------------------------------------------------------------------- */
function drawSeries(v){
  var x=begin(v), w=v.w, h=v.h;
  if(!series.length) return;
  var pad={l:78,r:78,t:20,b:36}, iw=w-pad.l-pad.r, ih=h-pad.t-pad.b;
  var live = groups.filter(function(g){ return g.members.some(function(s){ return s.on; }); });

  /* The first two unit-groups get a labelled axis, left and right. Beyond two there is nowhere
     to put a third without it overlapping the plot, so the caption says which ones are unlabelled
     rather than the chart implying they share the axis they are drawn against. */
  x.strokeStyle="#1c242f"; x.lineWidth=1;
  live.slice(0,2).forEach(function(g, gi){
    var tk = ticks(g.lo, g.hi, 5);
    x.fillStyle="#7f8b9a"; x.font=SMALL;
    tk.forEach(function(t){
      var yy = pad.t + ih*(1-(t-g.lo)/((g.hi-g.lo)||1));
      if(yy < pad.t-1 || yy > pad.t+ih+1) return;
      if(gi===0){ x.strokeStyle="#1c242f"; x.beginPath(); x.moveTo(pad.l,yy); x.lineTo(pad.l+iw,yy); x.stroke(); }
      x.textAlign = gi ? "left" : "right";
      x.fillText(nice(t), gi ? pad.l+iw+8 : pad.l-8, yy+3.5);
    });
    x.font=MONO; x.textAlign = gi ? "left" : "right";
    x.fillText(g.unit || "-", gi ? pad.l+iw+8 : pad.l-8, pad.t-6);
  });

  var t0 = F[0].t, t1 = F[N-1].t;
  x.fillStyle="#7f8b9a"; x.font=SMALL;
  ticks(t0, t1, 6).forEach(function(t){
    var xx = pad.l + iw*(t1>t0 ? (t-t0)/(t1-t0) : 0.5);
    x.fillRect(xx, pad.t+ih, 1, 4);
    x.textAlign="center"; x.fillText(nice(t), xx, pad.t+ih+17);
  });
  x.textAlign="right"; x.font=MONO; x.fillText("t / s", pad.l+iw, h-6);

  series.forEach(function(s, si){
    if(!s.on) return;
    var g = groups.filter(function(q){ return q.unit === s.unit; })[0];
    var span=(g.hi-g.lo)||1;
    x.strokeStyle=s.colour; x.lineWidth=2; x.lineJoin="round"; x.beginPath();
    var pen = false;
    s.vals.forEach(function(val,i){
      if(val !== val){ pen = false; return; }        /* a hole is a break, not a line to zero */
      var xx=pad.l+iw*(N>1?i/(N-1):0.5), yy=pad.t+ih*(1-(val-g.lo)/span);
      pen ? x.lineTo(xx,yy) : x.moveTo(xx,yy);
      pen = true;
    });
    x.stroke();
  });

  /* where we are */
  var cx=pad.l+iw*(N>1?frame/(N-1):0.5);
  x.strokeStyle="rgba(95,208,192,0.75)"; x.lineWidth=1.5;
  x.beginPath(); x.moveTo(cx,pad.t); x.lineTo(cx,pad.t+ih); x.stroke();

  v.hit = function(px, py){
    if(px < pad.l-4 || px > pad.l+iw+4) return null;
    var i = Math.round((px-pad.l)/iw*(N-1));
    if(i < 0 || i >= N) return null;
    var best = null, bd = 1e9;
    series.forEach(function(s){
      if(!s.on || s.vals[i] !== s.vals[i]) return;
      var g = groups.filter(function(q){ return q.unit === s.unit; })[0];
      var yy = pad.t+ih*(1-(s.vals[i]-g.lo)/((g.hi-g.lo)||1));
      var d = Math.abs(yy-py);
      if(d < bd){ bd = d; best = s; }
    });
    if(!best) return null;
    return {text:"t = "+nice(F[i].t)+" s   "+best.key+" = "+nice(best.vals[i])
                 +(best.unit ? " "+best.unit : "")};
  };

  var many = function(n, one, more){ return n + " " + (n === 1 ? one : (more || one + "s")); };
  var unlabelled = live.length > 2
    ? " · " + many(live.length-2, "further unit") + " drawn without an axis" : "";
  var holes = series.reduce(function(a,s){ return a + s.missing; }, 0);
  cap(v, many(series.length, "scalar")+" in "+many(groups.length, "unit group")
       + (groups.length > 1 ? " · series sharing a unit share an axis" : "")
       + " · the line marks the current frame"
       + unlabelled
       + (holes ? " · " + many(holes, "reading") + " absent and drawn as breaks" : ""));
}

function cap(v, text){
  /* keyed by slot, not by panel: a volume has two views of the same panel and each writes its
     own caption. Keying by panel made the second overwrite the first. */
  var el = document.getElementById("cap-" + v.slot);
  if (el) el.textContent = text;
}

/* ---- the legend, which is also the switch --------------------------------------------------- */
function buildLegend(){
  var box = document.getElementById("legend");
  if(!box || !series.length) return;
  series.forEach(function(s, i){
    /* Spread along the sequential scale rather than given arbitrary hues: the scale is
       already ordered and already distinguishable end to end, and a legend that reuses it
       means the same colour means the same thing everywhere on the page. */
    s.colour = css({table:SEQ}, series.length > 1 ? 0.12 + 0.76*i/(series.length-1) : 0.5);
    var b = document.createElement("button");
    b.setAttribute("aria-pressed", "true");
    b.innerHTML = "<i></i><span></span>";
    b.children[0].style.background = s.colour;
    s.label = b.children[1];
    b.onclick = function(){
      s.on = !s.on;
      b.setAttribute("aria-pressed", s.on ? "true" : "false");
      regroup(); mark("series"); render();
    };
    box.appendChild(b);
    s.button = b;
  });
}
function refreshLegend(){
  series.forEach(function(s){
    if(!s.label) return;
    var v = s.vals[frame];
    s.label.textContent = s.key + "  "
      + (v === v ? nice(v) + (s.unit ? " " + s.unit : "") : "absent");
  });
}

/* ---- what to redraw, and when ----------------------------------------------------------------
   Every view used to be redrawn on every pointermove of a drag, so rotating a cube also re-ran
   a 30 ms raycast for a volume nobody was touching. A view is redrawn when the frame moves, when
   the camera moves and it is a camera view, or when it is resized. */
var dirty = {};
function mark(which){
  views.forEach(function(v){
    if(which === "all"
       || (which === "frame")
       || (which === "camera" && (v.kind==="scene"||v.kind==="volume"||v.kind==="layout"))
       || which === v.slot) dirty[v.slot] = true;
  });
}
function drawOne(v){
  var f = F[frame];
  if(v.kind==="profile") drawProfile(v,f);
  else if(v.kind==="heatmap") drawHeat(v,f);
  else if(v.kind==="slices") drawSlices(v,f);
  else if(v.kind==="volume") drawVolume(v,f);
  else if(v.kind==="layout") drawLayout(v,f);
  else if(v.kind==="scene") drawScene(v,f);
  else drawSeries(v);
}
function render(){
  views.forEach(function(v){
    if(resize(v)) dirty[v.slot] = true;
    if(!dirty[v.slot]) return;
    dirty[v.slot] = false;
    drawOne(v);
  });
  refreshLegend();
  document.getElementById("tick").textContent =
    "frame "+(frame+1)+" / "+N+"   t = "+nice(F[frame].t)+" s";
  document.getElementById("scrub").value=String(frame);
}

/* ---- transport -------------------------------------------------------------------------------- */
var scrub=document.getElementById("scrub"), play=document.getElementById("play");
scrub.max=String(N-1);
function goto_(i){
  frame = ((i % N) + N) % N;
  mark("frame"); render();
}
function pause(){ playing=false; play.textContent="Play"; }
scrub.oninput=function(){ pause(); goto_(Number(scrub.value)); };
play.onclick=function(){ playing=!playing; play.textContent=playing?"Pause":"Play"; };
document.getElementById("back").onclick=function(){ pause(); goto_(frame-1); };
document.getElementById("fwd").onclick=function(){ pause(); goto_(frame+1); };
document.getElementById("speed").onchange=function(e){ speed = Number(e.target.value) || 1; };

/* Keys, because comparing two frames means pressing the same thing twice and a 12 px button is
   the wrong target for that. Ignored while a form control has focus, so the speed select still
   works with the keyboard. */
window.addEventListener("keydown", function(e){
  var tag = (e.target && e.target.tagName) || "";
  if(tag === "SELECT" || tag === "INPUT") return;
  if(e.key === " "){ e.preventDefault(); play.onclick(); }
  else if(e.key === "ArrowLeft"){ e.preventDefault(); pause(); goto_(frame-1); }
  else if(e.key === "ArrowRight"){ e.preventDefault(); pause(); goto_(frame+1); }
  else if(e.key === "Home"){ e.preventDefault(); pause(); goto_(0); }
  else if(e.key === "End"){ e.preventDefault(); pause(); goto_(N-1); }
});

/* ---- pointing at things ------------------------------------------------------------------------ */
/* **Save this view.** A report ends up in a document, and until now the only way to get a figure
   out of one was a screenshot at whatever the window happened to be. The canvas is already the
   size the display asked for, so `toDataURL` gives the full-resolution frame. Nothing leaves the
   machine: it is a data URL and an anchor, and the page has never had a network call in it. */
views.forEach(function(v){
  var b = document.getElementById("png-" + v.slot);
  if(!b) return;
  b.onclick = function(){
    var a = document.createElement("a");
    a.href = v.c.toDataURL("image/png");
    a.download = v.slot.replace(/[^\w.-]+/g, "_") + "-frame" + (frame+1) + ".png";
    a.click();
  };
});

views.forEach(function(v){
  var out = document.getElementById("read-" + v.slot);
  function at(e){
    var r = v.c.getBoundingClientRect();
    return [(e.clientX - r.left) * (v.w / (r.width || v.w)),
            (e.clientY - r.top) * (v.h / (r.height || v.h))];
  }
  v.c.addEventListener("pointermove", function(e){
    if(drag) return;
    if(!out || !v.hit) return;
    var q = at(e), got = v.hit(q[0], q[1]);
    out.textContent = got ? got.text : "";
  });
  v.c.addEventListener("pointerleave", function(){ if(out) out.textContent = ""; });
});

views.filter(function(v){
  return v.kind==="scene"||v.kind==="volume"||v.kind==="layout";
}).forEach(function(v){
  v.c.addEventListener("pointerdown",function(e){
    drag={x:e.clientX,y:e.clientY}; dragging=v; v.c.setPointerCapture(e.pointerId);
  });
  v.c.addEventListener("pointermove",function(e){
    if(!drag)return;
    cam.az+=(e.clientX-drag.x)*0.008;
    cam.el=Math.max(-1.5,Math.min(1.5,cam.el+(e.clientY-drag.y)*0.006));
    drag={x:e.clientX,y:e.clientY}; mark("camera"); render();
  });
  function stop(){
    if(!drag) return;
    drag=null;
    /* Settle to full quality once the pointer stops, so the draft resolution is never what a
       reader is left looking at. */
    dragging=null; mark("camera"); render();
  }
  v.c.addEventListener("pointerup", stop);
  v.c.addEventListener("pointercancel", stop);
  v.c.addEventListener("wheel",function(e){
    e.preventDefault();
    cam.dist=Math.max(1.2,Math.min(9,cam.dist*(1+e.deltaY*0.0011)));
    mark("camera"); render();
  },{passive:false});
});

if(window.ResizeObserver){
  var ro = new ResizeObserver(function(){ render(); });
  views.forEach(function(v){ ro.observe(v.c); });
} else {
  window.addEventListener("resize", function(){ mark("all"); render(); });
}

buildLegend();
if (window.matchMedia("(prefers-reduced-motion: reduce)").matches){ pause(); }
/* 16 frames a second at 1x. A run is not a film and a reader is comparing frames, not watching
   one; faster than this and a 12-frame run is over before it has been seen. */
function loop(now){
  if(playing && now-last > 62/speed){ last=now; frame=(frame+1)%N; mark("frame"); render(); }
  requestAnimationFrame(loop);
}
mark("all");
render();
requestAnimationFrame(loop);
})();
"##;
