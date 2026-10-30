/*!
 * DeGenPrime web app.
 * One file, no dependencies, no build step. The Rust engine runs as
 * WebAssembly in a Web Worker (website, WordPress plugin) or inside the local
 * `degenprime gui` server (desktop app); both answer the same JSON API.
 * Copyright (c) 2026 Richard Allen White III, RAW Lab, UNC Charlotte.
 * Free for academic and non-commercial use; see LICENSE.
 */
(function () {
  'use strict';
  if (!window.customElements || window.customElements.get('degenprime-app')) return;

  var CURRENT = document.currentScript;
  var SCRIPT_BASE = CURRENT && CURRENT.src ? CURRENT.src.replace(/[^/]*$/, '') : 'assets/';
  var VERSION = '2.3.0';
  var VIEWS = ['design', 'results', 'check', 'manual', 'about'];
  var STORE = 'degenprime.v2.';
  var REPO = 'https://github.com/raw-lab/degenprime';
  var MAFFT_DOI = 'https://doi.org/10.1093/molbev/mst010';
  // Above this size a loaded file is not copied into the text box: a
  // textarea holding megabytes of sequence makes the page crawl.
  var SHOW_LIMIT = 2 * 1024 * 1024;
  var ALIGN_DEFAULTS = { align: 'auto', align_mode: 'auto', adjust_direction: true };
  var ALIGN_MODES = [['auto', 'Automatic (mafft --auto)'], ['fftns2', 'FFT-NS-2 (fast)'], ['fftnsi', 'FFT-NS-i (iterative)'],
    ['linsi', 'L-INS-i (most accurate)'], ['ginsi', 'G-INS-i (global)'], ['einsi', 'E-INS-i (long insertions)']];
  var LAB = 'https://www.rawlab.org';
  // Licence notices of the built-in aligner (src/mafft/LICENSE-MAFFT.txt), which
  // must accompany the engine wherever it is distributed.
  var MAFFT_NOTICE = "rust-MAFFT: LICENSE-MIT\n-----------------------\n\nMIT License\n\nCopyright (c) 2026 Lucas Goiriz (luksgrin)\n\nPermission is hereby granted, free of charge, to any person obtaining a copy\nof this software and associated documentation files (the \"Software\"), to deal\nin the Software without restriction, including without limitation the rights\nto use, copy, modify, merge, publish, distribute, sublicense, and/or sell\ncopies of the Software, and to permit persons to whom the Software is\nfurnished to do so, subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\nOUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\nSOFTWARE.\n\n\nrust-MAFFT: LICENSE-BSD (MAFFT)\n-------------------------------\n\nBSD 3-Clause License\n\nThis license covers the algorithmic constructs ported from MAFFT\n(https://mafft.cbrc.jp/alignment/software/) and any verbatim or\nnear-verbatim translations of MAFFT's C source. Original Rust code,\ninfrastructure, and bindings authored for rust-MAFFT are MIT-licensed\n(see LICENSE-MIT).\n\nUpstream MAFFT:\n\n    MAFFT: multiple sequence alignment program\n    Copyright (c) 2009 Kazutaka Katoh\n\nRedistribution and use in source and binary forms, with or without\nmodification, are permitted provided that the following conditions are\nmet:\n\n1. Redistributions of source code must retain the above copyright\n   notice, this list of conditions and the following disclaimer.\n\n2. Redistributions in binary form must reproduce the above copyright\n   notice, this list of conditions and the following disclaimer in the\n   documentation and/or other materials provided with the\n   distribution.\n\n3. The name of the author may not be used to endorse or promote\n   products derived from this software without specific prior written\n   permission.\n\nTHIS SOFTWARE IS PROVIDED BY THE AUTHOR \"AS IS\" AND ANY EXPRESS OR\nIMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\nWARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\nDISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,\nINDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\n(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\nSERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)\nHOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,\nSTRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING\nIN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE\nPOSSIBILITY OF SUCH DAMAGE.\n";
  var DOI = 'https://doi.org/10.1093/bioadv/vbae044';

  var FALLBACK_DESIGN = {
    min_primer_len: 20, max_primer_len: 22, min_amplicon: 100, max_amplicon: null, begin: null, end: null,
    min_temp: 58, max_temp: 62, max_tm_diff: 5, min_gc: 40, max_gc: 60, max_gc_end: 80,
    max_degeneracy: 15, max_deletions: 15, max_poly_x: 4, max_self_any: 8, max_self_end: 3,
    max_pair_any: 8, max_pair_end: 3, primer_conc: 200, salt_conc: 50, max_primers: 5,
    min_3prime_distance: 3, rank_by: 'score'
  };
  var FALLBACK_CHECK = Object.assign({}, FALLBACK_DESIGN, { min_temp: 50, max_temp: 65, min_gc: 30, max_gc: 70 });
  var MASK = { A: 1, C: 2, G: 4, T: 8, U: 8, R: 5, Y: 10, S: 6, W: 9, K: 12, M: 3, B: 14, D: 13, H: 11, V: 7, N: 15 };
  var FORMAT_NAMES = { fasta: 'FASTA', clustal: 'CLUSTAL', raw: 'plain sequence' };
  var FILTER_NAMES = {
    degeneracy: 'Degeneracy', deletions: 'Gap columns', gc_content: 'GC content', gc_content_3prime: '3′ end GC',
    repeats: 'Repeats and runs', complementary_ends: 'Complementary ends', temperature: 'Melting temperature',
    self_complementarity: 'Self-complementarity'
  };
  function logoMark(assets) {
    return '<img class="mark" src="' + esc(assets) + 'mark.png" width="40" height="40" alt="" decoding="async">';
  }

  // [name, label, unit, min, max, step, optional]
  var PARAM_GROUPS = [
    { legend: 'Primer length', fields: [['min_primer_len', 'Shortest', 'nt', 18, 25, 1], ['max_primer_len', 'Longest', 'nt', 18, 25, 1]] },
    { legend: 'Melting temperature', fields: [['min_temp', 'Min', '°C', 50, 65, 0.5], ['max_temp', 'Max', '°C', 50, 65, 0.5], ['max_tm_diff', 'Max ΔTm', '°C', 0, 20, 0.5]] },
    { legend: 'GC content', fields: [['min_gc', 'Min', '%', 0, 100, 1], ['max_gc', 'Max', '%', 0, 100, 1], ['max_gc_end', 'Last 5 bases', '%', 0, 100, 20]] },
    { legend: 'Product size', fields: [['min_amplicon', 'Min', 'bp', 1, 1000000, 10], ['max_amplicon', 'Max', 'bp', 1, 1000000, 10, 'optional']] },
    { legend: 'Design region', note: 'Alignment columns, 1-based. Leave empty to search the whole alignment.', fields: [['region_from', 'From column', '', 1, 100000000, 1, 'optional'], ['region_to', 'To column', '', 1, 100000000, 1, 'optional']] },
    { legend: 'Reaction', fields: [['primer_conc', 'Primer', 'nM', 50, 100000, 10], ['salt_conc', 'Na⁺', 'mM', 1, 5000, 1]] }
  ];
  var ADVANCED = [
    ['max_degeneracy', 'Degenerate positions', '%', 0, 100, 1], ['max_deletions', 'Gap columns', '%', 0, 100, 1],
    ['max_poly_x', 'Longest single-base run', 'nt', 1, 12, 1], ['min_3prime_distance', 'Min 3′ distance, pairs', 'nt', 0, 50, 1],
    ['max_self_any', 'Self-compl., anywhere', '', 0, 40, 1], ['max_self_end', 'Self-compl., 3′ end', '', 0, 40, 1],
    ['max_pair_any', 'Pair compl., anywhere', '', 0, 40, 1], ['max_pair_end', 'Pair compl., 3′ end', '', 0, 40, 1]
  ];
  var OUTPUT_FIELDS = [['max_primers', 'Pairs to return', '', 1, 50, 1]];
  var CHECK_FIELDS = [['primer_conc', 'Primer', 'nM', 50, 100000, 10], ['salt_conc', 'Na⁺', 'mM', 1, 5000, 1], ['min_temp', 'Tm min', '°C', 50, 65, 0.5], ['max_temp', 'Tm max', '°C', 50, 65, 0.5], ['min_gc', 'GC min', '%', 0, 100, 1], ['max_gc', 'GC max', '%', 0, 100, 1]];
  var INT_FIELDS = { min_primer_len: 1, max_primer_len: 1, min_amplicon: 1, max_amplicon: 1, region_from: 1, region_to: 1, max_poly_x: 1, min_3prime_distance: 1, max_primers: 1 };
  var FIELD_META = {};
  PARAM_GROUPS.forEach(function (g) { g.fields.forEach(function (f) { FIELD_META[f[0]] = { f: f, legend: g.legend }; }); });
  ADVANCED.forEach(function (f) { FIELD_META[f[0]] = { f: f, legend: 'Advanced' }; });
  OUTPUT_FIELDS.forEach(function (f) { FIELD_META[f[0]] = { f: f, legend: 'Output' }; });

  // ------------------------------------------------------------------ helpers
  function esc(s) {
    return String(s == null ? '' : s).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }
  function fmt(v, d) { return (v == null || !isFinite(v)) ? '–' : Number(v).toFixed(d == null ? 1 : d); }
  function fmtInt(v) { return (v == null || !isFinite(v)) ? '–' : Math.round(Number(v)).toLocaleString('en-US'); }
  function plural(n, word) { return fmtInt(n) + ' ' + word + (n === 1 ? '' : 's'); }
  function save(key, value) { try { localStorage.setItem(STORE + key, JSON.stringify(value)); } catch (e) { /* storage off */ } }
  function restore(key) { try { var v = localStorage.getItem(STORE + key); return v ? JSON.parse(v) : null; } catch (e) { return null; } }
  function baseCls(c) { return (c === 'A' || c === 'C' || c === 'G' || c === 'T') ? 'b' + c : (c === '-' ? 'bgap' : 'bx'); }
  function seqHTML(s) {
    var o = '';
    String(s || '').split('').forEach(function (c) { o += '<span class="' + baseCls(c) + '">' + esc(c) + '</span>'; });
    return o;
  }
  function prop(k, v) { return '<div><dt>' + esc(k) + '</dt><dd>' + esc(v) + '</dd></div>'; }

  // ------------------------------------------------------------------ engines
  // The WebAssembly core. Self-contained so it can also run inside a Worker
  // created from its source text.
  function makeCore() {
    var inst = null, loading = null, lastPanic = '';
    var enc = new TextEncoder(), dec = new TextDecoder();
    function ensure(url) {
      if (inst) return Promise.resolve();
      if (!loading) {
        var imports = { env: { dgp_panic: function (p, l) {
          try { lastPanic = dec.decode(new Uint8Array(inst.exports.memory.buffer, p, l).slice()); } catch (e) { lastPanic = 'unknown error'; }
        } } };
        loading = (function () {
          var streaming = WebAssembly.instantiateStreaming
            ? WebAssembly.instantiateStreaming(fetch(url), imports)
            : Promise.reject(new Error('no streaming'));
          return streaming.catch(function () {
            return fetch(url).then(function (r) {
              if (!r.ok) throw new Error('could not download the engine (' + r.status + ')');
              return r.arrayBuffer();
            }).then(function (buf) { return WebAssembly.instantiate(buf, imports); });
          });
        })().then(function (res) { inst = res.instance; }, function (err) { loading = null; throw err; });
      }
      return loading;
    }
    function put(bytes) {
      var p = inst.exports.dgp_alloc(bytes.length);
      new Uint8Array(inst.exports.memory.buffer, p, bytes.length).set(bytes);
      return p;
    }
    function call(ep, body) {
      var e = enc.encode(ep), b = enc.encode(body);
      var ep_p = put(e), b_p = put(b), out;
      try { out = inst.exports.dgp_call(ep_p, e.length, b_p, b.length); } catch (err) {
        inst = null; loading = null;
        throw new Error('The engine stopped unexpectedly' + (lastPanic ? ' (' + lastPanic + ')' : '') + '. Please report this input to the developers.');
      }
      var mem = inst.exports.memory.buffer;
      var len = new DataView(mem).getUint32(out, true);
      var text = dec.decode(new Uint8Array(mem, out + 4, len).slice());
      inst.exports.dgp_free(out, len + 4);
      inst.exports.dgp_free(ep_p, e.length);
      inst.exports.dgp_free(b_p, b.length);
      return text;
    }
    return { ensure: ensure, call: call };
  }

  function WasmEngine(url) {
    this.kind = 'wasm';
    this.url = new URL(url, location.href).href;
    this.worker = null; this.core = null; this.failed = false;
    this.pending = new Map(); this.seq = 0;
  }
  WasmEngine.prototype.startWorker = function () {
    if (this.worker || this.failed) return;
    var self = this;
    try {
      var src = 'var core=(' + makeCore.toString() + ')();self.onmessage=function(e){var d=e.data;core.ensure(d.url).then(function(){' +
        'self.postMessage({id:d.id,ok:true,out:core.call(d.ep,d.body)});}).catch(function(err){self.postMessage({id:d.id,ok:false,err:String(err&&err.message||err)});});};';
      this.worker = new Worker(URL.createObjectURL(new Blob([src], { type: 'text/javascript' })));
      this.worker.onmessage = function (e) {
        var p = self.pending.get(e.data.id);
        if (!p) return;
        self.pending.delete(e.data.id);
        if (e.data.ok) p.resolve(e.data.out); else p.reject(new Error(e.data.err));
      };
      this.worker.onerror = function (e) { if (e && e.preventDefault) e.preventDefault(); self.fallback(); };
    } catch (e) {
      this.failed = true; this.worker = null;
    }
  };
  // Workers can be blocked (strict CSP, file://). Run on the main thread instead.
  WasmEngine.prototype.fallback = function () {
    var self = this;
    this.failed = true;
    if (this.worker) { try { this.worker.terminate(); } catch (e) { /* ignore */ } this.worker = null; }
    var waiting = Array.from(this.pending.values());
    this.pending.clear();
    waiting.forEach(function (p) { self.direct(p.ep, p.body).then(p.resolve, p.reject); });
  };
  WasmEngine.prototype.direct = function (ep, body) {
    if (!this.core) this.core = makeCore();
    var core = this.core;
    return core.ensure(this.url).then(function () { return core.call(ep, body); });
  };
  // Stop whatever the worker is doing (a long alignment, say): terminate it
  // and start a fresh one on the next call.
  WasmEngine.prototype.cancel = function () {
    var waiting = Array.from(this.pending.values());
    this.pending.clear();
    if (this.worker) { try { this.worker.terminate(); } catch (e) { /* ignore */ } this.worker = null; }
    waiting.forEach(function (p) { p.reject(new Error('cancelled')); });
    return waiting.length > 0;
  };
  WasmEngine.prototype.call = function (ep, payload) {
    var self = this, body = JSON.stringify(payload || {});
    this.startWorker();
    var text = this.worker
      ? new Promise(function (resolve, reject) {
          var id = ++self.seq;
          self.pending.set(id, { resolve: resolve, reject: reject, ep: ep, body: body });
          self.worker.postMessage({ id: id, url: self.url, ep: ep, body: body });
        })
      : this.direct(ep, body);
    return text.then(JSON.parse);
  };

  function NativeEngine(token) { this.kind = 'native'; this.token = token; this.inflight = new Set(); }
  NativeEngine.prototype.call = function (ep, payload) {
    var self = this, ctl = window.AbortController ? new AbortController() : null;
    if (ctl) this.inflight.add(ctl);
    var done = function () { if (ctl) self.inflight.delete(ctl); };
    return fetch('/api/' + ep, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', 'X-DGP-Token': this.token || '' },
      body: JSON.stringify(payload || {}),
      signal: ctl ? ctl.signal : undefined
    }).then(function (r) {
      done();
      if (r.status === 403) throw new Error('the desktop app was restarted; reload this page');
      if (!r.ok) throw new Error('the desktop app answered ' + r.status + '; restart degenprime gui and reload this page');
      return r.json();
    }, function (err) {
      done();
      if (err && err.name === 'AbortError') throw new Error('cancelled');
      throw new Error('the desktop app is not running; start it again with degenprime gui');
    });
  };
  // The desktop app finishes the job in the background; the page stops waiting.
  NativeEngine.prototype.cancel = function () {
    var n = this.inflight.size;
    this.inflight.forEach(function (c) { c.abort(); });
    this.inflight.clear();
    return n > 0;
  };

  function ensureFontFace(base) {
    if (document.getElementById('dgp-font-face')) return;
    var st = document.createElement('style');
    st.id = 'dgp-font-face';
    st.textContent = '@font-face{font-family:"DGP Recursive";src:url("' + base + 'fonts/recursive-latin-mono-normal.woff2") format("woff2");font-weight:300 1000;font-style:normal;font-display:swap}';
    (document.head || document.documentElement).appendChild(st);
  }

  // ------------------------------------------------------------------ forms
  function fieldHTML(f, value) {
    var v = value == null ? '' : value, optional = f[6] === 'optional';
    return '<label class="num"><span class="lab">' + esc(f[1]) + '</span><span class="box"><input type="number" name="' + f[0] +
      '" value="' + esc(v) + '" min="' + f[3] + '" max="' + f[4] + '" step="' + f[5] + '" inputmode="decimal"' +
      (optional ? ' placeholder="none"' : '') + '>' + (f[2] ? '<span class="unit">' + esc(f[2]) + '</span>' : '') + '</span></label>';
  }

  function readDesignForm(form, base) {
    var out = Object.assign({}, base), errors = [], ui = {};
    Object.keys(FIELD_META).forEach(function (name) {
      var el = form.elements[name];
      if (!el) return;
      var meta = FIELD_META[name], f = meta.f, raw = String(el.value).trim();
      if (raw === '') {
        if (f[6] === 'optional') ui[name] = null;
        else errors.push(meta.legend + ': ' + f[1].toLowerCase() + ' is required.');
        return;
      }
      var n = Number(raw);
      if (!isFinite(n)) { errors.push(meta.legend + ': ' + f[1].toLowerCase() + ' must be a number.'); return; }
      if (INT_FIELDS[name] && n < f[3]) { errors.push(meta.legend + ': ' + f[1].toLowerCase() + ' must be at least ' + f[3] + '.'); return; }
      if (n < 0) { errors.push(meta.legend + ': ' + f[1].toLowerCase() + ' cannot be negative.'); return; }
      ui[name] = INT_FIELDS[name] ? Math.round(n) : n;
    });
    Object.keys(ui).forEach(function (k) { if (k !== 'region_from' && k !== 'region_to') out[k] = ui[k]; });
    if ('region_from' in ui) out.begin = ui.region_from == null ? null : Math.max(0, ui.region_from - 1);
    if ('region_to' in ui) out.end = ui.region_to == null ? null : ui.region_to;
    if (form.elements.rank_by) out.rank_by = form.elements.rank_by.value === 'tm_diff' ? 'tm_diff' : 'score';
    if (out.begin != null && out.end != null && out.end <= out.begin) errors.push('Design region: the last column must come after the first.');
    if (out.max_amplicon != null && out.max_amplicon < out.min_amplicon) errors.push('Product size: the maximum is smaller than the minimum.');
    return { settings: out, errors: errors };
  }

  // ------------------------------------------------------------------ charts
  function niceStep(span, target) {
    var raw = span / target, pow = Math.pow(10, Math.floor(Math.log10(Math.max(raw, 1e-9))));
    var steps = [1, 2, 5, 10];
    for (var i = 0; i < steps.length; i++) if (raw <= steps[i] * pow) return steps[i] * pow;
    return 10 * pow;
  }

  function arrow(a, b, y, dir) {
    var c = (a + b) / 2, w = Math.max(8, b - a), l = c - w / 2, r = c + w / 2, head = Math.min(6, w * 0.5), h = 4;
    var pts = dir > 0
      ? [[l, y - h], [r - head, y - h], [r, y], [r - head, y + h], [l, y + h]]
      : [[r, y - h], [l + head, y - h], [l, y], [l + head, y + h], [r, y + h]];
    return '<polygon points="' + pts.map(function (p) { return p[0].toFixed(1) + ',' + p[1].toFixed(1); }).join(' ') + '"/>';
  }

  // Conservation along the alignment, with primer lanes underneath.
  function skylineSVG(cons, pairs, sel, compact) {
    var cv = cons && cons.conservation;
    if (!cv || !cv.length) return '';
    var L = cv.length, codes = cons.codes || '';
    var W = 1000, padL = compact ? 1 : 44, padR = compact ? 1 : 10, plotW = W - padL - padR;
    var barTop = compact ? 2 : 12, barH = compact ? 46 : 112, base = barTop + barH;
    var laneH = 18, lanesTop = base + 30, lanes = compact ? 0 : pairs.length;
    var H = compact ? base + 2 : lanesTop + lanes * laneH + 8;
    var x = function (c) { return padL + (c / L) * plotW; };
    var bins = Math.max(1, Math.min(L, Math.floor(plotW / (compact ? 1.6 : 2.2))));
    var g = '', bands = '';
    if (!compact) {
      var s = pairs.filter(function (p) { return p.rank === sel; })[0];
      if (s) [[s.forward_from - 1, s.forward_to], [s.reverse_from - 1, s.reverse_to]].forEach(function (ab) {
        bands += '<rect class="footband" x="' + x(ab[0]).toFixed(1) + '" y="' + (barTop - 6) + '" width="' + Math.max(3, x(ab[1]) - x(ab[0])).toFixed(1) + '" height="' + (barH + 6) + '"/>';
      });
    }
    for (var b = 0; b < bins; b++) {
      var c0 = Math.floor((b * L) / bins), c1 = Math.max(c0 + 1, Math.floor(((b + 1) * L) / bins));
      var sum = 0, n = 0;
      for (var c = c0; c < c1; c++) if (codes[c] !== '-') { sum += cv[c]; n++; }
      var x0 = x(c0), w = Math.max(0.5, x(c1) - x0 - (bins > 300 ? 0.25 : 0.6));
      if (!n) { g += '<rect class="gapcol" x="' + x0.toFixed(1) + '" y="' + (base - 2) + '" width="' + w.toFixed(2) + '" height="2"/>'; continue; }
      var v = sum / n / 100, h = Math.max(1.5, v * barH);
      g += '<rect class="cbar' + (v >= 0.999 ? ' full' : '') + '" x="' + x0.toFixed(1) + '" y="' + (base - h).toFixed(1) + '" width="' + w.toFixed(2) +
        '" height="' + h.toFixed(1) + '" fill-opacity="' + (0.22 + 0.78 * v * v).toFixed(2) + '"/>';
    }
    g += bands;
    g += '<line class="axis" x1="' + padL + '" x2="' + (W - padR) + '" y1="' + (base + 0.5) + '" y2="' + (base + 0.5) + '"/>';
    if (!compact) {
      g += '<text class="ylab" x="' + (padL - 8) + '" y="' + (barTop + 9) + '" text-anchor="end">100%</text>';
      g += '<text class="ylab" x="' + (padL - 8) + '" y="' + base + '" text-anchor="end">0</text>';
      var step = niceStep(L, 10);
      for (var t = step; t <= L; t += step) {
        var xx = x(t - 0.5).toFixed(1);
        g += '<line class="tick" x1="' + xx + '" x2="' + xx + '" y1="' + base + '" y2="' + (base + 5) + '"/><text class="xlab" x="' + xx + '" y="' + (base + 18) + '" text-anchor="middle">' + fmtInt(t) + '</text>';
      }
      pairs.forEach(function (p, i) {
        var y = lanesTop + i * laneH + laneH / 2;
        var f0 = x(p.forward_from - 1), f1 = x(p.forward_to), r0 = x(p.reverse_from - 1), r1 = x(p.reverse_to);
        var cls = p.rank === sel ? 'lane sel' : (sel ? 'lane dim' : 'lane');
        g += '<g class="' + cls + '" data-action="select" data-rank="' + p.rank + '">';
        g += '<rect class="hit" x="0" y="' + (y - laneH / 2) + '" width="' + W + '" height="' + laneH + '"/>';
        g += '<text class="lanelab" x="' + (padL - 8) + '" y="' + (y + 4) + '" text-anchor="end">' + p.rank + '</text>';
        g += '<line class="prod" x1="' + f1.toFixed(1) + '" x2="' + r0.toFixed(1) + '" y1="' + y + '" y2="' + y + '"/>';
        g += arrow(f0, f1, y, 1) + arrow(r0, r1, y, -1);
        g += '<title>Pair ' + p.rank + ': ' + p.amplicon_size + ' bp; forward ' + p.forward_from + '–' + p.forward_to + ', reverse ' + p.reverse_from + '–' + p.reverse_to + '</title></g>';
      });
    }
    return '<svg class="skyline' + (compact ? ' compact' : '') + '" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Conservation along the alignment' +
      (compact ? '' : ', with primer positions') + '">' + g + '</svg>';
  }

  // A simulated agarose gel: ladder plus one lane per pair (log-size migration).
  function gelSVG(pairs, sel) {
    var sizes = pairs.map(function (p) { return p.amplicon_size; });
    var maxP = Math.max.apply(null, [100].concat(sizes));
    var LADDER = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000];
    var hi = LADDER.filter(function (v) { return v >= maxP * 1.3; })[0] || Math.ceil(maxP * 1.3);
    var lo = 75, ladder = LADDER.filter(function (v) { return v <= hi; });
    var laneW = 46, gap = 12, left = 50, top = 22, gelH = 250, lanes = pairs.length + 1;
    var W = left + lanes * (laneW + gap) + 4, H = top + gelH + 30;
    var y = function (bp) {
      bp = Math.min(hi, Math.max(lo, bp));
      return top + 20 + ((Math.log(hi) - Math.log(bp)) / (Math.log(hi) - Math.log(lo))) * (gelH - 34);
    };
    var lx = function (i) { return left + i * (laneW + gap); };
    var s = '<svg class="gel" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Predicted products on a simulated agarose gel">' +
      '<defs><filter id="dgpglow" x="-20%" y="-300%" width="140%" height="700%"><feGaussianBlur stdDeviation="1.7" result="b"/><feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter></defs>' +
      '<rect class="gel-bg" x="0" y="0" width="' + W + '" height="' + H + '" rx="8"/>';
    for (var i = 0; i < lanes; i++) s += '<rect class="well" x="' + (lx(i) + 4) + '" y="' + top + '" width="' + (laneW - 8) + '" height="7" rx="1.5"/>';
    var labelled = { 100: 1, 200: 1, 300: 1, 500: 1, 1000: 1, 1500: 1, 2000: 1, 3000: 1, 5000: 1, 10000: 1 };
    ladder.forEach(function (bp) {
      s += '<rect class="band ladder" x="' + (lx(0) + 5) + '" y="' + (y(bp) - 1.5).toFixed(1) + '" width="' + (laneW - 10) + '" height="3" rx="1.5"/>';
      if (labelled[bp]) s += '<text class="gel-lab" x="' + (left - 8) + '" y="' + (y(bp) + 3.5).toFixed(1) + '" text-anchor="end">' + (bp >= 1000 ? (bp / 1000) + ' kb' : bp) + '</text>';
    });
    pairs.forEach(function (p, k) {
      var xx = lx(k + 1), yy = y(p.amplicon_size), on = p.rank === sel ? ' sel' : '';
      s += '<g class="gel-lane' + on + '" data-action="select" data-rank="' + p.rank + '"><rect class="hit" x="' + xx + '" y="' + top + '" width="' + laneW + '" height="' + gelH + '"/>' +
        '<rect class="band prod" x="' + (xx + 5) + '" y="' + (yy - 2.5).toFixed(1) + '" width="' + (laneW - 10) + '" height="5" rx="2" filter="url(#dgpglow)"/>' +
        '<title>Pair ' + p.rank + ': ' + p.amplicon_size + ' bp</title></g>';
    });
    s += '<text class="gel-lane-lab" x="' + (lx(0) + laneW / 2) + '" y="' + (H - 10) + '" text-anchor="middle">L</text>';
    pairs.forEach(function (p, k) {
      s += '<text class="gel-lane-lab' + (p.rank === sel ? ' sel' : '') + '" x="' + (lx(k + 1) + laneW / 2) + '" y="' + (H - 10) + '" text-anchor="middle">' + p.rank + '</text>';
    });
    return s + '</svg>';
  }

  // Forward versus reverse Tm, with the Tm window as a box.
  function tmSVG(pairs, set, sel) {
    var tms = [];
    pairs.forEach(function (p) { tms.push(p.forward_tm, p.reverse_tm); });
    var lo = Math.floor(Math.min.apply(null, [set.min_tm].concat(tms)) - 1);
    var hi = Math.ceil(Math.max.apply(null, [set.max_tm].concat(tms)) + 1);
    var S = 220, L = 44, T = 10, W = L + S + 14, H = T + S + 40;
    var sx = function (v) { return L + ((v - lo) / (hi - lo)) * S; };
    var sy = function (v) { return T + S - ((v - lo) / (hi - lo)) * S; };
    var s = '<svg class="tmplot" viewBox="0 0 ' + W + ' ' + H + '" role="img" aria-label="Forward versus reverse primer Tm">' +
      '<rect class="plot-bg" x="' + L + '" y="' + T + '" width="' + S + '" height="' + S + '"/>' +
      '<rect class="tm-window" x="' + sx(set.min_tm).toFixed(1) + '" y="' + sy(set.max_tm).toFixed(1) + '" width="' + (sx(set.max_tm) - sx(set.min_tm)).toFixed(1) +
      '" height="' + (sy(set.min_tm) - sy(set.max_tm)).toFixed(1) + '"/>' +
      '<line class="diag" x1="' + sx(lo) + '" y1="' + sy(lo) + '" x2="' + sx(hi) + '" y2="' + sy(hi) + '"/>';
    var step = (hi - lo) > 8 ? 2 : 1;
    for (var v = Math.ceil(lo / step) * step; v <= hi; v += step) {
      s += '<text class="axlab" x="' + sx(v).toFixed(1) + '" y="' + (T + S + 15) + '" text-anchor="middle">' + v + '</text>' +
        '<text class="axlab" x="' + (L - 6) + '" y="' + (sy(v) + 4).toFixed(1) + '" text-anchor="end">' + v + '</text>';
    }
    s += '<text class="axtitle" x="' + (L + S / 2) + '" y="' + (H - 4) + '" text-anchor="middle">Forward primer Tm, °C</text>' +
      '<text class="axtitle" transform="translate(11 ' + (T + S / 2) + ') rotate(-90)" text-anchor="middle">Reverse primer Tm, °C</text>';
    pairs.forEach(function (p) {
      var on = p.rank === sel, cx = sx(p.forward_tm), cy = sy(p.reverse_tm);
      s += '<g class="pt' + (on ? ' sel' : '') + '" data-action="select" data-rank="' + p.rank + '"><circle cx="' + cx.toFixed(1) + '" cy="' + cy.toFixed(1) + '" r="' + (on ? 6 : 4.5) + '"/>' +
        '<text x="' + (cx + 8).toFixed(1) + '" y="' + (cy - 5 - ((p.rank - 1) % 3) * 9).toFixed(1) + '">' + p.rank + '</text>' +
        '<title>Pair ' + p.rank + ': forward ' + fmt(p.forward_tm, 1) + ' °C, reverse ' + fmt(p.reverse_tm, 1) + ' °C</title></g>';
    });
    return s + '</svg>';
  }

  // ------------------------------------------------------------------ result pieces
  function siteRowHTML(primer, consensus, seq) {
    var out = '';
    for (var i = 0; i < seq.length; i++) {
      var c = seq[i], p = primer[i] || 'N', k = consensus[i] || 'N';
      if (p === '-') { out += c === '-' ? '<span class="gap">-</span>' : '<span class="ins" title="Base in a column the primer skips">' + esc(c) + '</span>'; continue; }
      if (c === k) { out += '<span class="same">.</span>'; continue; }
      if (c === 'N') { out += '<span class="gap">n</span>'; continue; }
      if (c === '-') { out += '<span class="miss" title="Gap under the primer">-</span>'; continue; }
      var cm = MASK[c] || 15, pm = MASK[p] || 15;
      out += (cm & pm) === cm
        ? '<span class="var ' + baseCls(c) + '" title="Covered by the degenerate base ' + esc(p) + '">' + esc(c) + '</span>'
        : '<span class="miss" title="Not covered by the primer">' + esc(c) + '</span>';
    }
    return out;
  }

  function siteBlockHTML(title, b) {
    var h = '<div class="site"><h4>' + esc(title) + '<span class="cols">columns ' + fmtInt(b.from) + '–' + fmtInt(b.to) + '</span></h4><div class="site-grid">';
    h += '<span class="sid primer">Primer</span><span class="sq primer seq">' + seqHTML(b.primer) + '</span>';
    b.rows.forEach(function (row) {
      h += '<span class="sid" title="' + esc(row.id) + '">' + esc(row.id) + '</span><span class="sq">' + siteRowHTML(b.primer, b.consensus, row.seq) + '</span>';
    });
    if (b.hidden) h += '<span class="sid">+' + fmtInt(b.hidden) + ' more</span><span></span>';
    return h + '</div></div>';
  }

  function sitesHTML(site) {
    return '<div class="sites">' + siteBlockHTML('Forward primer', site.forward) + siteBlockHTML('Reverse primer, read 5′→3′', site.reverse) +
      '<p class="legend"><span><span class="k same">.</span>same as the consensus</span><span><span class="k bG">G</span>variant covered by the degenerate code</span>' +
      '<span><span class="k miss">A</span>not covered</span><span><span class="k">-</span>gap</span></p></div>';
  }

  function explainHTML(e) {
    if (!e) return '';
    var rows = [['Candidates considered', 'considered'], ['Too degenerate', 'degeneracy'], ['Too many gap columns', 'deletions'], ['GC content out of range', 'gc_content'],
      ['3′ end too GC-rich', 'gc_3prime'], ['Repeats or long runs', 'repeats'], ['Complementary ends', 'complementary_ends'], ['Tm too low', 'tm_low'],
      ['Tm too high', 'tm_high'], ['Self-complementary', 'self_complementarity'], ['Passed every filter', 'ok']];
    var h = '<div class="explain-grid"><table class="stats"><thead><tr><th scope="col">Primers</th><th scope="col" class="n">Forward</th><th scope="col" class="n">Reverse</th></tr></thead><tbody>';
    rows.forEach(function (r) {
      h += '<tr' + (r[1] === 'ok' || r[1] === 'considered' ? ' class="strong"' : '') + '><th scope="row">' + r[0] + '</th><td class="n">' + fmtInt(e.forward[r[1]]) + '</td><td class="n">' + fmtInt(e.reverse[r[1]]) + '</td></tr>';
    });
    var p = e.pairs;
    var prow = [['Forward × reverse combinations', p.combinations], ['Reverse primer upstream or overlapping', p.overlap], ['Product too short', p.product_too_short], ['Product too long', p.product_too_long],
      ['Examined in detail', p.examined], ['ΔTm too large', p.tm_diff], ['Primers could form a dimer', p.pair_complementarity], ['Near-duplicate of a better pair', p.too_similar], ['Returned', p.ok]];
    h += '</tbody></table><table class="stats"><thead><tr><th scope="col">Pairs</th><th scope="col" class="n">Count</th></tr></thead><tbody>';
    prow.forEach(function (r, i) {
      h += '<tr' + (i === 0 || i === 4 || i === prow.length - 1 ? ' class="strong"' : '') + '><th scope="row">' + r[0] + '</th><td class="n">' + fmtInt(r[1]) + '</td></tr>';
    });
    return h + '</tbody></table></div><p class="hint">Each candidate is counted under the first filter it fails. The first four pair counts cover every combination; the search then examines only pairs that could still beat the ones already found, and the remaining counts refer to those.</p>';
  }

  function oligoCard(title, a) {
    var pass = a.filters.every(function (f) { return f.passed; });
    var range = a.fold > 1 && a.tm_min != null ? '<span class="range">' + fmtInt(a.fold) + ' variants, ' + fmt(a.tm_min, 1) + '–' + fmt(a.tm_max, 1) + ' °C</span>' : '';
    return '<article class="oligo"><h2>' + esc(title) + '</h2><p class="oligo-seq seq">5′-' + seqHTML(a.sequence) + '-3′</p>' +
      '<div class="readout"><span class="lab">Tm</span><span class="val">' + fmt(a.tm, 1) + '</span><span class="unit">°C</span>' + range + '</div><dl class="props">' +
      prop('Length', a.actual_size + ' nt') + prop('GC content', fmt(a.gc_content, 1) + ' %') +
      prop('Degenerate positions', a.degeneracy_count + (a.fold > 1 ? ' (' + fmtInt(a.fold) + '-fold)' : '')) +
      prop('Self-complementarity', fmt(a.self_any, 0) + ' anywhere, ' + fmt(a.self_end, 0) + ' at 3′') +
      prop('3′ end stability', fmt(a.end_stability, 1) + ' kcal/mol') + prop('Hairpin heuristic', a.hairpin_score + ' of 3') +
      '</dl><h3>Filters <span class="verdict ' + (pass ? 'ok' : 'bad') + '">' + (pass ? 'all passed' : 'some failed') + '</span></h3><ul class="checks">' +
      a.filters.map(function (f) {
        return '<li class="' + (f.passed ? 'ok' : 'bad') + '"><span class="mark" aria-hidden="true">' + (f.passed ? '✓' : '✗') + '</span><span>' +
          esc(FILTER_NAMES[f.name] || f.name) + '</span><span class="why">' + esc(f.reason || 'passed') + '</span></li>';
      }).join('') + '</ul>' +
      (a.self_dimer ? '<h3>Strongest 3′ self-dimer</h3><pre class="duplex">' + esc(a.self_dimer) + '</pre>' : '') + '</article>';
  }

  function pairCard(p) {
    return '<article class="pair-check"><h2>The pair</h2><dl class="props">' + prop('ΔTm', fmt(p.tm_diff, 1) + ' °C') +
      prop('Pair complementarity', fmt(p.pair_any, 0) + ' anywhere, ' + fmt(p.pair_end, 0) + ' at 3′') + prop('3′–3′ dimer', p.dimer_score + ' bp') + '</dl>' +
      (p.diagram ? '<h3>Strongest 3′ cross-dimer</h3><pre class="duplex">' + esc(p.diagram) + '</pre>' : '<p class="hint">No 3′ cross-dimer found.</p>') +
      '<p class="hint">When designing, DeGenPrime drops pairs that score above 8 anywhere or 3 at the 3′ end, the Primer3 defaults.</p></article>';
  }

  // ------------------------------------------------------------------ static pages
  var MANUAL = '<div class="manual"><nav class="toc" aria-label="Manual sections">' +
    [['m-quick', 'Quick start'], ['m-input', 'Input'], ['m-how', 'How primers are chosen'], ['m-params', 'Parameters'], ['m-results', 'Reading the results'],
      ['m-tm', 'Melting temperature'], ['m-cli', 'Command line and desktop app'], ['m-privacy', 'Privacy']].map(function (t) {
      return '<a href="#/manual" data-action="jump" data-target="' + t[0] + '">' + t[1] + '</a>';
    }).join('') + '</nav><article class="prose">' +
    '<h1>Manual</h1><p class="lede">DeGenPrime turns an alignment of related DNA sequences into ranked, degenerate PCR primer pairs. This page explains what it does and how to read the results.</p>' +
    '<h2 id="m-quick">Quick start</h2><ol>' +
    '<li>On <a href="#/design">Design</a>, paste aligned sequences in FASTA or CLUSTAL format, choose a file, or load the example.</li>' +
    '<li>Keep the default parameters or adjust them. They suit standard Taq PCR with 20–22 nt primers.</li>' +
    '<li>Select <b>Design primers</b>. Results lists the best pairs, shows where they bind and how every sequence matches them, and offers CSV, report, FASTA and JSON downloads.</li></ol>' +
    '<h2 id="m-input">Input</h2>' +
    '<p>DeGenPrime reads FASTA (gaps as <code>-</code> or <code>.</code>), CLUSTAL as written by Clustal Omega, MAFFT and MUSCLE, or a single plain sequence. Spaces, digits and <code>*</code> are ignored, <code>U</code> is read as <code>T</code>, and protein sequences are refused with a message.</p>' +
    '<p>Aligned or not: when the sequences differ in length, DeGenPrime aligns them first with MAFFT, built in as a pure-Rust port (rust-MAFFT) that gives the same alignments as C MAFFT. Nothing needs to be installed, and on the website the alignment also runs in your browser. <b>Automatic</b> picks the strategy the way <code>mafft --auto</code> does: L-INS-i for small sets, FFT-NS-i for medium ones, FFT-NS-2 for large ones. Sequences that lie on the opposite strand are reverse-complemented first (as <code>mafft --adjustdirection</code>), and the Results page says which. Download the alignment from Results to check it. Choose <b>Never</b> to pad unaligned sequences with gaps instead, as version 2.2 did.</p>' +
    '<p>Ambiguity codes in your sequences are respected: an <code>R</code> counts half as A and half as G, and widens the primer’s code to cover both. <code>N</code> is treated as missing data.</p>' +
    '<h2 id="m-how">How primers are chosen</h2><ol>' +
    '<li><b>Consensus.</b> Each alignment column becomes one position: the IUPAC code covering every base seen there, the most common base, and its conservation, the share of sequences carrying it. Gaps count against conservation, and a column where most sequences have a gap is left out of the primer.</li>' +
    '<li><b>Candidates.</b> Every window of the allowed primer lengths, on both strands, inside the design region.</li>' +
    '<li><b>Filters</b>, in this order: degeneracy, gap columns, GC content, GC in the last five bases, single-base runs and dinucleotide repeats, complementary ends, melting temperature, and self-complementarity. Each candidate is counted under the first filter it fails, in <b>Filter statistics</b> on the Results page.</li>' +
    '<li><b>Pairing.</b> Forward and reverse primers must not overlap, must give a product within the size limits, and must differ in Tm by no more than the maximum ΔTm. Pairs are ranked by' +
    '<div class="formula">score = P<sub>f</sub> + P<sub>r</sub> + |Tm<sub>f</sub> − Tm<sub>opt</sub>| + |Tm<sub>r</sub> − Tm<sub>opt</sub>| + |ΔTm|</div>' +
    'where P = Σ(1 − conservation) × 100 over the primer and Tm<sub>opt</sub> is the middle of the Tm window. One column at 80% conservation costs 20 points, so the most conserved primers come first and Tm breaks ties. Lower is better.</li>' +
    '<li><b>Final checks.</b> Pairs whose primers could anneal to each other are dropped, and a pair is skipped when both of its 3′ ends lie within 3 nt of a better pair’s.</li></ol>' +
    '<h2 id="m-params">Parameters</h2><div class="table-scroll"><table><thead><tr><th>Setting</th><th>Default</th><th>Command line</th><th>Meaning</th></tr></thead><tbody>' +
    [['Primer length', '20–22 nt', '--min_primer_len --max_primer_len', 'Allowed range 18–25 nt.'],
      ['Melting temperature', '58–62 °C', '--min_temp --max_temp', 'Tm window, allowed range 50–65 °C. Its middle is the optimum used for ranking.'],
      ['Max ΔTm', '5 °C', '--max_tm_diff', 'Largest Tm difference within a pair.'],
      ['GC content', '40–60 %', '--min_gc --max_gc', 'Over the whole primer.'],
      ['GC in the last 5 bases', '80 %', '--max_gc_end', 'At most four G or C in the 3′ pentamer.'],
      ['Product size', 'at least 100 bp', '--amplicon --max_amplicon', 'In alignment columns. No maximum unless you set one.'],
      ['Design region', 'whole alignment', '--begin --end', 'Keeps both primers, and so the product, inside these columns.'],
      ['Primer and Na⁺', '200 nM, 50 mM', '--primer_conc --salt_conc', 'Reaction conditions for the Tm model.'],
      ['Pairs to return', '5', '--max_primers', 'From 1 to 50.'],
      ['Rank pairs by', 'conservation, then Tm', '--rank score|tm_diff', 'tm_diff restores the version 2.1 ordering, smallest ΔTm first.'],
      ['Align', 'when needed', '--align auto|always|never', 'When to align with the built-in MAFFT: when lengths differ, always (realign), or never (pad).'],
      ['MAFFT strategy', 'automatic', '--align_mode auto|fftns2|fftnsi|linsi|ginsi|einsi', 'Automatic follows mafft --auto.'],
      ['Opposite strand', 'reverse-complement', '--keep_strand', 'Sequences on the other strand are flipped before aligning unless you keep them.'],
      ['Degenerate positions', '15 %', '--max_degeneracy', 'Share of ambiguous positions allowed in a primer.'],
      ['Gap columns', '15 %', '--max_deletions', 'Share of mostly-gap columns allowed under a primer.'],
      ['Longest single-base run', '4 nt', '--max_poly_x', 'Runs such as AAAAA are rejected.'],
      ['Complementarity limits', '8, 3, 8, 3', '--max_self_any --max_self_end --max_pair_any --max_pair_end', 'Self- and pair complementarity, anywhere and at the 3′ end. These are the Primer3 defaults.'],
      ['Min 3′ distance', '3 nt', '--min_3prime_distance', 'Skips near-duplicate pairs; 0 turns this off.']].map(function (r) {
      return '<tr><td>' + r[0] + '</td><td>' + r[1] + '</td><td>' + r[2].split(' ').map(function (c) { return '<code>' + esc(c) + '</code>'; }).join(' ') + '</td><td>' + r[3] + '</td></tr>';
    }).join('') + '</tbody></table></div>' +
    '<h2 id="m-results">Reading the results</h2><ul>' +
    '<li><b>Location</b> is the primer’s footprint in alignment columns, 1-based and inclusive, on the forward strand. Reverse primers are written 5′→3′, as you would order them.</li>' +
    '<li><b>Degenerate</b> counts ambiguous positions; the figure in brackets is the fold, the number of distinct oligos in the degenerate mix.</li>' +
    '<li><b>Self compl.</b> and <b>pair compl.</b> are Primer3-style alignment scores, anywhere and at the 3′ end: one point per base pair, minus one per mismatch.</li>' +
    '<li><b>Product size</b> is counted in alignment columns. Insertions and deletions in a particular sequence change its real product size.</li>' +
    '<li><b>Binding sites</b> show every sequence under each primer: a dot matches the consensus, a coloured letter is a variant the degenerate code covers, and a red letter is not covered.</li></ul>' +
    '<h2 id="m-tm">Melting temperature</h2>' +
    '<p>Tm uses the SantaLucia (1998) unified nearest-neighbour parameters, the salt correction ΔS + 0.368 × (N − 1) × ln[Na⁺], and a quarter of the primer concentration. No initiation terms are added, which matches the published DeGenPrime, so values stay comparable with the original tool; expect differences of 1–3 °C from Primer3. A degenerate base uses the average of the stacks it stands for, and <a href="#/check">Check a primer</a> also reports the range across every variant.</p>' +
    '<h2 id="m-cli">Command line and desktop app</h2>' +
    '<pre><code>degenprime design -i sequences.fasta -o results --detailed --fasta\ndegenprime align -i sequences.fasta -o aligned.fasta   # MAFFT, pure Rust\ndegenprime test -p GCTACCTGAACACCATCCCG --partner CCCTGGATGTTGTASTCGCC\ndegenprime gui      # this interface, running on your computer\ndegenprime doctor   # checks your installation</code></pre>' +
    '<p>The command line writes the same CSV and JSON you can download here. Its <code>--begin</code> and <code>--end</code> are 0-based, with the end excluded, as in earlier versions; this page uses 1-based columns.</p>' +
    '<h2 id="m-privacy">Privacy</h2>' +
    '<p>On the website the engine, alignment included, runs in your browser as WebAssembly, so sequences are never uploaded or stored. Your parameter choices are remembered in this browser only. The desktop app runs entirely on your computer.</p>' +
    '</article></div>';

  function aboutHTML(native, assets) {
    var cite = 'Fulghum B, Tanker SH, White RA III (2024). DeGenPrime provides robust primer design and optimization unlocking the biosphere. Bioinformatics Advances 4(1): vbae044. doi:10.1093/bioadv/vbae044';
    return '<div class="about-grid"><article class="prose"><img class="about-logo" src="' + esc(assets) + 'logo.webp" width="280" height="270" alt="DeGenPrime: an orange robot holding a DNA double helix" decoding="async">' +
      '<h1>About DeGenPrime</h1>' +
      '<p class="lede">DeGenPrime designs PCR primers that work across a group of related sequences: a gene family, a viral lineage, or the variants of a marker gene in an environmental survey. Where the sequences differ, the primer carries an IUPAC ambiguity code, so one degenerate primer matches them all.</p>' +
      '<p>DeGenPrime is written in Rust. The same engine powers the command line, the desktop app and this website, where it runs in your browser as WebAssembly. Sequences that are not aligned yet are aligned with MAFFT, built in as a pure-Rust port.</p>' +
      '<h2>How to cite</h2><p>If DeGenPrime contributes to published work, please cite:</p>' +
      '<div class="cite">Fulghum B, Tanker SH, White RA III (2024). DeGenPrime provides robust primer design and optimization unlocking the biosphere. <i>Bioinformatics Advances</i> 4(1): vbae044. <a href="' + DOI + '" target="_blank" rel="noopener">doi:10.1093/bioadv/vbae044</a></div>' +
      '<div class="cite">White RA III (2026). DeGenPrime-2.0: a robust, streamlined, and accurate primer design via functional programming in Rust. <i>bioRxiv</i>, preprint.</div>' +
      '<p><button type="button" class="btn ghost" data-action="copy" data-label="the citation" data-text="' + esc(cite) + '">Copy citation</button></p>' +
      '<p>If DeGenPrime aligned your sequences, please also cite MAFFT:</p>' +
      '<div class="cite">Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. <i>Molecular Biology and Evolution</i> 30(4): 772–780. <a href="' + MAFFT_DOI + '" target="_blank" rel="noopener">doi:10.1093/molbev/mst010</a></div>' +
      '<h2>Licence</h2><p>Free for academic, research and other non-commercial use. Commercial use needs written permission from the author; see the LICENSE file in the source code. The built-in aligner is rust-MAFFT by Lucas Goiriz (MIT), a port of MAFFT by Kazutaka Katoh (BSD-3-Clause).</p>' +
      '<details class="notices"><summary>Licence notices of the built-in aligner</summary><pre>' + esc(MAFFT_NOTICE) + '</pre></details></article>' +
      '<aside class="facts-card"><dl>' +
      '<div><dt>Version</dt><dd>' + VERSION + '</dd></div>' +
      '<div><dt>Engine</dt><dd>' + (native ? 'Native Rust, on this computer' : 'Rust compiled to WebAssembly, in your browser') + '</dd></div>' +
      '<div><dt>Aligner</dt><dd>MAFFT 7.526 algorithms, pure Rust (rust-MAFFT 0.2.0)</dd></div>' +
      '<div><dt>Author</dt><dd>Richard Allen White III, RAW Lab, University of North Carolina at Charlotte</dd></div>' +
      '<div><dt>Source code</dt><dd><a href="' + REPO + '" target="_blank" rel="noopener">github.com/raw-lab/degenprime</a></dd></div>' +
      '<div><dt>Lab</dt><dd><a href="' + LAB + '" target="_blank" rel="noopener">www.rawlab.org</a></dd></div>' +
      '<div><dt>Contact</dt><dd><a href="mailto:rwhit101@charlotte.edu">rwhit101@charlotte.edu</a></dd></div>' +
      '</dl></aside></div>';
  }

  // ------------------------------------------------------------------ the element
  class DegenPrimeApp extends HTMLElement {
    connectedCallback() {
      if (this._ready) return;
      this._ready = true;
      var cfg = window.DGP_CONFIG || {};
      this.assets = (this.getAttribute('data-assets') || SCRIPT_BASE).replace(/\/?$/, '/');
      this.host = this.getAttribute('data-host') || 'site';
      this.native = cfg.engine === 'native';
      this.engine = this.native ? new NativeEngine(cfg.token) : new WasmEngine(this.assets + 'degenprime.wasm');
      var v = this.getAttribute('data-view');
      this.initialView = VIEWS.indexOf(v) >= 0 ? v : 'design';
      this.defaults = { design: FALLBACK_DESIGN, check: FALLBACK_CHECK };
      this.state = {
        view: this.initialView, input: '', filename: '', bigInput: false, inspect: null, result: null, selected: 0, busy: false,
        design: Object.assign({}, FALLBACK_DESIGN, restore('design') || {}),
        align: Object.assign({}, ALIGN_DEFAULTS, restore('align') || {}),
        check: { primer: '', partner: '', settings: Object.assign({}, FALLBACK_CHECK, restore('check') || {}), result: null, error: '' }
      };
      this.inspectSeq = 0;
      ensureFontFace(this.assets);
      this.root = this.attachShadow({ mode: 'open' });
      this.root.innerHTML = '<style>:host{display:block}.dgp{visibility:hidden}</style><link rel="stylesheet" href="' + esc(this.assets) + 'degenprime.css">' + this.shell();
      var self = this;
      setTimeout(function () { var d = self.$('.dgp'); if (d) d.style.visibility = 'visible'; }, 2500);
      this.bind();
      this._onHash = function () { self.route(true); };
      window.addEventListener('hashchange', this._onHash);
      this.route(false);
      this.loadDefaults();
    }

    disconnectedCallback() { if (this._onHash) window.removeEventListener('hashchange', this._onHash); }
    $(sel) { return this.root.querySelector(sel); }
    $$(sel) { return Array.prototype.slice.call(this.root.querySelectorAll(sel)); }

    shell() {
      var tabs = [['design', 'Design'], ['results', 'Results'], ['check', 'Check a primer'], ['manual', 'Manual'], ['about', 'About']];
      return '<div class="dgp" data-host="' + esc(this.host) + '"><header class="top">' +
        '<a class="brand" href="#/design" aria-label="DeGenPrime: design primers">' + logoMark(this.assets) +
        '<span class="wordmark"><b>DeGenPrime</b><span class="basebar" aria-hidden="true"><i></i><i></i><i></i><i></i></span></span></a>' +
        '<nav class="tabs" aria-label="DeGenPrime">' + tabs.map(function (t) { return '<a href="#/' + t[0] + '" data-view="' + t[0] + '">' + t[1] + '</a>'; }).join('') + '</nav>' +
        '<span class="engine" data-engine title="' + (this.native ? 'Calculations run in the degenprime program on this computer.' : 'Calculations run in this browser tab with WebAssembly. Sequences are not uploaded.') + '">' +
        (this.native ? 'Desktop app, on this computer' : 'Runs in your browser') + '</span></header>' +
        '<div class="progress" aria-hidden="true"></div><main class="view" tabindex="-1"></main>' +
        '<footer class="foot"><span>DeGenPrime ' + VERSION + '. Sequences are analysed ' + (this.native ? 'on this computer' : 'in your browser') + ' and never uploaded.</span>' +
        '<nav aria-label="Links"><a href="' + REPO + '" target="_blank" rel="noopener">Source code</a><a href="' + LAB + '" target="_blank" rel="noopener">RAW Lab</a><a href="#/about">How to cite</a></nav></footer>' +
        '<div class="toast" role="status" aria-live="polite"></div></div>';
    }

    route(fromHash) {
      var h = (location.hash || '').replace(/^#\/?/, '').split(/[/?&]/)[0];
      if (VIEWS.indexOf(h) >= 0) this.state.view = h;
      else if (!fromHash) this.state.view = this.initialView;
      else return;
      this.render();
      if (fromHash) {
        if (this.getBoundingClientRect().top < 0) this.scrollIntoView({ block: 'start' });
        var m = this.$('main');
        if (m) m.focus({ preventScroll: true });
      }
    }

    go(view) {
      if (location.hash === '#/' + view) this.route(true);
      else location.hash = '#/' + view;
    }

    render() {
      var v = this.state.view;
      this.$$('.tabs a').forEach(function (a) {
        if (a.getAttribute('data-view') === v) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current');
      });
      var html = v === 'results' ? this.viewResults() : v === 'check' ? this.viewCheck() : v === 'manual' ? MANUAL : v === 'about' ? aboutHTML(this.native, this.assets) : this.viewDesign();
      var main = this.$('main');
      main.innerHTML = (this.engineError ? '<p class="notice bad">' + esc(this.engineError) + '</p>' : '') + html;
      main.setAttribute('data-view', v);
      if (v === 'design') this.afterDesign();
    }

    loadDefaults() {
      var self = this;
      return this.engine.call('defaults').then(function (d) {
        if (d && d.design) self.defaults = { design: d.design, check: d.check };
      }).catch(function (e) {
        self.engineError = 'The primer engine could not start: ' + e.message + '.' +
          (self.native ? '' : ' If this site sends a Content-Security-Policy header, it must allow \'wasm-unsafe-eval\'.');
        var st = self.$('[data-engine]');
        if (st) { st.classList.add('err'); st.textContent = 'Engine unavailable'; }
        self.render();
      });
    }

    // ------------------------------------------------------------ design
    viewDesign() {
      var s = this.state;
      return '<section class="design"><h1>Design degenerate primers from related sequences</h1>' +
        '<p class="lede">Paste related DNA sequences, aligned or not: unaligned sequences are aligned with the built-in MAFFT. DeGenPrime finds conserved stretches on both strands, writes the variable positions as IUPAC codes and ranks primer pairs by conservation, then by melting temperature.</p>' +
        '<form class="design-form" data-form="design" novalidate><div class="seq-col">' +
        '<div class="field-head"><label for="dgp-input">Sequences</label><span class="hint">FASTA or CLUSTAL, aligned or not. A single sequence works too.</span></div>' +
        '<div class="drop" data-drop><textarea id="dgp-input" class="seqs" name="input" wrap="off" spellcheck="false" autocomplete="off" autocapitalize="off" ' +
        'placeholder="&gt;sequence_1&#10;ATGCCGCACCACGAGTTCGAGTGCAGCAAGGTG&#10;&gt;sequence_2&#10;ATGCCGTACCACACCTTCGAGTGCAGCGAGTGC"></textarea></div>' +
        '<div class="seq-actions"><label class="btn ghost">Choose file<input type="file" class="vh" data-file accept=".fa,.fasta,.fas,.fna,.ffn,.frn,.aln,.clust,.clustal,.txt,text/plain"></label>' +
        '<button type="button" class="btn ghost" data-action="example">Load example</button><button type="button" class="btn ghost" data-action="clear">Clear</button>' +
        '<span class="filename" data-filename>' + esc(s.filename) + '</span></div>' +
        '<div class="inspect" data-inspect aria-live="polite">' + this.inspectHTML() + '</div></div>' +
        '<div class="param-col">' + this.paramsHTML(s.design) + '</div></form></section>';
    }

    paramsHTML(set) {
      var ui = Object.assign({}, set, { region_from: set.begin != null ? set.begin + 1 : null, region_to: set.end != null ? set.end : null });
      var h = PARAM_GROUPS.map(function (g) {
        return '<fieldset><legend>' + g.legend + '</legend><div class="row">' + g.fields.map(function (f) { return fieldHTML(f, ui[f[0]]); }).join('') + '</div>' +
          (g.note ? '<p class="note">' + g.note + '</p>' : '') + '</fieldset>';
      }).join('');
      h += '<fieldset><legend>Output</legend><div class="row">' + fieldHTML(OUTPUT_FIELDS[0], set.max_primers) +
        '<label class="selectfield"><span class="lab">Rank pairs by</span><select name="rank_by">' +
        '<option value="score"' + (set.rank_by !== 'tm_diff' ? ' selected' : '') + '>Conservation, then Tm</option>' +
        '<option value="tm_diff"' + (set.rank_by === 'tm_diff' ? ' selected' : '') + '>Smallest ΔTm (2.1)</option></select></label></div></fieldset>';
      h += '<details class="advanced"><summary>Advanced filters</summary><div class="grid2">' + ADVANCED.map(function (f) { return fieldHTML(f, set[f[0]]); }).join('') +
        '</div><p class="note">Complementarity limits are Primer3-style alignment scores, one point per base pair. Primer3’s defaults are 8 anywhere and 3 at the 3′ end.</p></details>';
      var al = this.state.align;
      h += '<fieldset class="align-set"><legend>Alignment</legend><div class="row">' +
        '<label class="selectfield"><span class="lab">Align with MAFFT</span><select name="align">' +
        [['auto', 'When lengths differ'], ['always', 'Always (realign)'], ['never', 'Never (pad with gaps)']].map(function (o) {
          return '<option value="' + o[0] + '"' + (al.align === o[0] ? ' selected' : '') + '>' + o[1] + '</option>';
        }).join('') + '</select></label>' +
        '<label class="selectfield"><span class="lab">Strategy</span><select name="align_mode">' +
        ALIGN_MODES.map(function (o) { return '<option value="' + o[0] + '"' + (al.align_mode === o[0] ? ' selected' : '') + '>' + o[1] + '</option>'; }).join('') +
        '</select></label></div><label class="checkline"><input type="checkbox" name="adjust_direction"' + (al.adjust_direction ? ' checked' : '') + '> Reverse-complement sequences on the opposite strand</label>' +
        '<p class="note">MAFFT is built in (pure Rust), so nothing needs installing' + (this.native ? '' : ', and the alignment runs in your browser') + '.</p></fieldset>';
      return h + '<div class="run"><button type="submit" class="btn primary">Design primers</button><button type="button" class="btn ghost" data-action="reset">Reset</button>' +
        '<button type="button" class="btn ghost cancel" data-action="cancel" hidden>Cancel</button></div>' +
        '<p class="form-error" data-error role="alert"></p>';
    }

    afterDesign() {
      this.fillTextarea();
      if (this.state.input.trim() && !this.state.inspect) this.runInspect();
      if (this.state.busy) this.setBusy(true);
    }

    // Large files stay out of the text box (a textarea with megabytes in
    // it makes the page crawl); the box says what is loaded instead.
    fillTextarea() {
      var ta = this.$('#dgp-input'), s = this.state;
      if (!ta) return;
      if (s.bigInput) {
        ta.value = '';
        ta.readOnly = true;
        var size = ' (' + fmt(s.input.length / 1048576, 1) + ' MB)';
        ta.placeholder = (s.filename ? s.filename + size + ' is loaded' : 'The sequences' + size + ' are loaded') + '\nbut not shown, to keep the page responsive.\nSelect Clear to paste other sequences.';
      } else {
        ta.readOnly = false;
        ta.placeholder = '>sequence_1\nATGCCGCACCACGAGTTCGAGTGCAGCAAGGTG\n>sequence_2\nATGCCGTACCACACCTTCGAGTGCAGCGAGTGC';
        if (ta.value !== s.input) ta.value = s.input;
      }
    }

    inspectHTML() {
      var s = this.state;
      if (!s.input.trim()) return '<p class="hint">Paste sequences, drop a file onto the box, or load the example: five <i>anfD</i> nitrogenase genes aligned with MAFFT.</p>';
      var r = s.inspect;
      if (!r) return '<p class="hint">Reading the sequences…</p>';
      if (!r.ok) return '<p class="notice bad">' + esc(r.error || 'These sequences could not be read.') + '</p>';
      if (!r.input) return '';
      var c = r.consensus, inp = r.input, called = 0, full = 0;
      var notes = (r.notes || []).map(function (w) { return '<p class="notice info">' + esc(w) + '</p>'; }).join('');
      if (!c) {
        return '<p class="summary"><b>' + plural(inp.sequences, 'sequence') + '</b>, up to ' + fmtInt(inp.columns) + ' nt, ' + (FORMAT_NAMES[inp.format] || inp.format) + ', not aligned.</p>' + notes +
          (this.state.align.align === 'never' ? '<p class="notice warn">Alignment is set to Never, so these sequences will be padded with gaps instead. Change it under Alignment.</p>' : '');
      }
      for (var i = 0; i < c.codes.length; i++) if (c.codes[i] !== '-') { called++; if (c.conservation[i] >= 100) full++; }
      var pct = called ? Math.round((100 * full) / called) : 0;
      return '<p class="summary"><b>' + plural(inp.sequences, 'sequence') + '</b>, ' + fmtInt(inp.columns) + ' alignment columns, ' + (FORMAT_NAMES[inp.format] || inp.format) + '. ' +
        (inp.sequences > 1 ? pct + '% of columns are identical in every sequence.' : 'Primers from a single sequence will not be degenerate.') + '</p>' +
        skylineSVG(c, [], 0, true) + '<div class="strip-axis"><span>1</span><span>conservation along the alignment</span><span>' + fmtInt(inp.columns) + '</span></div>' +
        notes + (r.warnings || []).map(function (w) { return '<p class="notice warn">' + esc(w) + '</p>'; }).join('');
    }

    paintInspect() { var el = this.$('[data-inspect]'); if (el) el.innerHTML = this.inspectHTML(); }

    scheduleInspect() {
      var self = this;
      clearTimeout(this._inspectTimer);
      this._inspectTimer = setTimeout(function () { self.runInspect(); }, 400);
    }

    runInspect() {
      var self = this, text = this.state.input, token = ++this.inspectSeq;
      if (!text.trim()) { this.state.inspect = null; this.paintInspect(); return Promise.resolve(); }
      return this.engine.call('inspect', { input: text, filename: this.state.filename || null }).then(function (r) {
        if (token === self.inspectSeq) { self.state.inspect = r; self.paintInspect(); }
      }).catch(function (e) {
        if (token !== self.inspectSeq) return;
        if (e.message === 'cancelled') { self.state.inspect = null; return; }
        self.state.inspect = { ok: false, error: e.message };
        self.paintInspect();
      });
    }

    setInput(text, filename) {
      this.state.input = text;
      this.state.filename = filename || '';
      this.state.bigInput = text.length > SHOW_LIMIT;
      this.state.inspect = null;
      this.fillTextarea();
      var fn = this.$('[data-filename]');
      if (fn) fn.textContent = this.state.filename;
      this.paintInspect();
      return this.runInspect();
    }

    readFile(file) {
      var self = this;
      var max = this.native ? 256 : 64;
      if (file.size > max * 1024 * 1024) { this.toast('That file is larger than ' + max + ' MB' + (this.native ? '; use the command line for it.' : '; use the desktop app or the command line for it.')); return; }
      var rd = new FileReader();
      rd.onload = function () { self.setInput(String(rd.result || ''), file.name); };
      rd.onerror = function () { self.toast('Could not read ' + file.name + '.'); };
      rd.readAsText(file);
    }

    loadExample() {
      var self = this;
      return fetch(this.assets + 'example_anfD.clust').then(function (r) {
        if (!r.ok) throw new Error('HTTP ' + r.status);
        return r.text();
      }).then(function (t) { return self.setInput(t, 'example_anfD.clust'); }).catch(function (e) {
        self.toast('The example could not be loaded (' + e.message + ').');
      });
    }

    setBusy(on, label) {
      this.state.busy = on;
      var d = this.$('.dgp');
      if (d) d.classList.toggle('busy', on);
      var b = this.$('form[data-form="design"] button[type="submit"]');
      if (b) { b.disabled = on; b.textContent = on ? (label || 'Designing…') : 'Design primers'; }
      var c = this.$('[data-action="cancel"]');
      if (c) c.hidden = !on;
    }

    readAlign(form) {
      var a = Object.assign({}, this.state.align);
      if (form.elements.align) a.align = form.elements.align.value;
      if (form.elements.align_mode) a.align_mode = form.elements.align_mode.value;
      if (form.elements.adjust_direction) a.adjust_direction = !!form.elements.adjust_direction.checked;
      this.state.align = a;
      save('align', a);
      return a;
    }

    cancelDesign() {
      if (!this.state.busy) return;
      this.engine.cancel();
      this.runToken = (this.runToken || 0) + 1;
      this.setBusy(false);
      var e = this.$('[data-error]');
      if (e) e.textContent = 'Cancelled.';
      this.state.inspect = null;
      if (this.state.input.trim()) this.runInspect();
    }

    runDesign(form) {
      var self = this, err = this.$('[data-error]');
      var read = readDesignForm(form, this.state.design);
      if (read.errors.length) { err.textContent = read.errors.join(' '); return; }
      if (!this.state.input.trim()) { err.textContent = 'Paste sequences or load the example first.'; this.$('#dgp-input').focus(); return; }
      err.textContent = '';
      this.state.design = read.settings;
      save('design', read.settings);
      var al = this.readAlign(form);
      var insp = this.state.inspect, unaligned = insp && insp.input && insp.input.aligned_input === false && insp.input.sequences > 1;
      var aligning = al.align === 'always' || (al.align === 'auto' && unaligned);
      var t0 = performance.now(), token = this.runToken = (this.runToken || 0) + 1;
      this.setBusy(true, aligning ? 'Aligning and designing…' : 'Designing…');
      return this.engine.call('design', {
        input: this.state.input, settings: read.settings, filename: this.state.filename || null,
        generated_at: new Date().toISOString(), align: al.align, align_mode: al.align_mode, adjust_direction: al.adjust_direction
      }).then(function (r) {
        if (token !== self.runToken) return;
        self.setBusy(false);
        r.elapsed = performance.now() - t0;
        if (!r.ok) { var e2 = self.$('[data-error]'); if (e2) e2.textContent = r.error || 'The design failed.'; return; }
        self.state.result = r;
        self.state.selected = r.pairs.length ? r.pairs[0].rank : 0;
        self.go('results');
      }).catch(function (e) {
        if (token !== self.runToken || e.message === 'cancelled') return;
        self.setBusy(false);
        var e2 = self.$('[data-error]');
        if (e2) e2.textContent = e.message;
      });
    }

    // ------------------------------------------------------------ results
    viewResults() {
      var r = this.state.result;
      if (!r) {
        return '<section class="empty"><h1>No results yet</h1><p class="lede">Design primers from an alignment and the ranked pairs appear here, with where they bind and how each sequence matches them.</p>' +
          '<p><a class="btn primary" href="#/design">Design primers</a></p></section>';
      }
      var inp = r.input || {}, pairs = r.pairs || [], sel = this.state.selected;
      var h = '<header class="res-head"><div><h1>' + (pairs.length ? plural(pairs.length, 'primer pair') : 'No primer pairs passed') + '</h1>' +
        '<p class="sub">' + (inp.filename ? esc(inp.filename) : 'Pasted sequences') + ': ' + plural(inp.sequences, 'sequence') + ', ' + fmtInt(inp.columns) + ' alignment columns' +
        (inp.aligned_with ? ', aligned with ' + esc(inp.aligned_with) : '') + '. Designed in ' + fmt((r.elapsed || 0) / 1000, 2) + ' s ' + (this.native ? 'on this computer' : 'in your browser') + '.</p></div>';
      var dl = pairs.length ? [['csv', 'Download CSV'], ['txt', 'Report'], ['fasta', 'FASTA'], ['json', 'JSON']] : [];
      if (r.alignment) dl.push(['alignment', 'Alignment']);
      if (dl.length) {
        h += '<div class="downloads">' + dl.map(function (d) {
          return '<button type="button" class="btn ghost" data-action="download" data-kind="' + d[0] + '">' + d[1] + '</button>';
        }).join('') + '</div>';
      }
      h += '</header>' + (r.notes || []).map(function (w) { return '<p class="notice info">' + esc(w) + '</p>'; }).join('') +
        (r.warnings || []).map(function (w) { return '<p class="notice warn">' + esc(w) + '</p>'; }).join('');
      if (pairs.length) {
        h += '<section><h2>Where the primers bind</h2><p class="hint">Bar height is conservation across the alignment. Arrows are primers and dashed lines the products between them; select a pair to highlight it.</p>' +
          '<div class="map scroll-x" data-map>' + skylineSVG(r.consensus, pairs, sel, false) + '</div></section>';
        h += '<section><h2>Primer pairs</h2><div class="table-wrap">' + this.pairsTable(r) + '</div></section>';
        h += '<section class="charts" data-charts>' + this.chartsHTML(r) + '</section>';
      }
      return h + '<details class="explain"' + (pairs.length ? '' : ' open') + '><summary>Filter statistics</summary>' + explainHTML(r.explain) + '</details>';
    }

    pairsTable(r) {
      var self = this, sel = this.state.selected;
      var h = '<table class="pairs"><thead><tr><th scope="col">Primer</th><th scope="col">Sequence, 5′→3′</th><th scope="col">Location</th><th scope="col" class="n">Length</th>' +
        '<th scope="col" class="n">Tm °C</th><th scope="col" class="n">GC %</th><th scope="col" class="n">Degenerate</th>' +
        '<th scope="col" class="n" title="Primer3-style self-complementarity, anywhere / at the 3′ end">Self compl.</th><th scope="col"><span class="vh">Copy</span></th></tr></thead>';
      r.pairs.forEach(function (p) {
        var site = (r.sites || []).filter(function (s) { return s.rank === p.rank; })[0];
        h += '<tbody class="pair' + (p.rank === sel ? ' sel' : '') + '" data-rank="' + p.rank + '"><tr class="pair-head"><th colspan="9" scope="rowgroup">' +
          '<button type="button" class="pair-select" data-action="select" data-rank="' + p.rank + '">Pair ' + p.rank + '</button><span class="facts">' +
          '<span><b>' + fmtInt(p.amplicon_size) + '</b> bp product</span><span>ΔTm <b>' + fmt(p.tm_difference, 1) + '</b> °C</span>' +
          '<span>pair compl. <b>' + fmt(p.pair_any, 0) + '</b> / <b>' + fmt(p.pair_end, 0) + '</b></span><span>score <b>' + fmt(p.score, 1) + '</b></span></span></th></tr>' +
          self.primerRow('Forward', p.forward_primer, p.forward_from, p.forward_to, p.forward_length, p.forward_tm, p.forward_gc, p.forward_degeneracy, p.forward_fold, p.forward_self_any, p.forward_self_end, 'pair ' + p.rank + ' forward primer') +
          self.primerRow('Reverse', p.reverse_primer, p.reverse_from, p.reverse_to, p.reverse_length, p.reverse_tm, p.reverse_gc, p.reverse_degeneracy, p.reverse_fold, p.reverse_self_any, p.reverse_self_end, 'pair ' + p.rank + ' reverse primer');
        if (site) {
          h += '<tr class="sites-row"><td colspan="9"><details data-sites="' + p.rank + '"><summary>Binding sites in ' + plural(site.forward.rows.length + site.forward.hidden, 'sequence') + '</summary>' + sitesHTML(site) + '</details></td></tr>';
        }
        h += '</tbody>';
      });
      return h + '</table>';
    }

    primerRow(label, seq, from, to, len, tm, gc, deg, fold, sany, send, copyLabel) {
      return '<tr><th scope="row">' + label + '</th><td class="seq">' + seqHTML(seq) + '</td><td class="loc">' + fmtInt(from) + '–' + fmtInt(to) + '</td>' +
        '<td class="n">' + len + '</td><td class="n">' + fmt(tm, 1) + '</td><td class="n">' + fmt(gc, 0) + '</td>' +
        '<td class="n">' + (deg ? deg + ' <span class="fold">(' + fmtInt(fold) + '×)</span>' : '0') + '</td><td class="n">' + fmt(sany, 0) + ' / ' + fmt(send, 0) + '</td>' +
        '<td><button type="button" class="btn tiny" data-action="copy" data-text="' + esc(seq) + '" data-label="' + esc(copyLabel) + '">Copy</button></td></tr>';
    }

    chartsHTML(r) {
      var sel = this.state.selected, set = r.settings || { min_tm: 58, max_tm: 62 };
      return '<figure>' + gelSVG(r.pairs, sel) + '<figcaption><b>Expected products.</b> A simulated agarose gel: lane L is a DNA ladder, then one lane per pair. Real products vary with insertions and deletions in each sequence.</figcaption></figure>' +
        '<figure>' + tmSVG(r.pairs, set, sel) + '<figcaption><b>Tm balance.</b> Forward against reverse primer Tm. The dashed box is your Tm window; points near the diagonal anneal at the same temperature.</figcaption></figure>';
    }

    updateSelection() {
      var r = this.state.result, sel = this.state.selected;
      if (!r) return;
      var map = this.$('[data-map]');
      if (map) map.innerHTML = skylineSVG(r.consensus, r.pairs, sel, false);
      var ch = this.$('[data-charts]');
      if (ch) ch.innerHTML = this.chartsHTML(r);
      this.$$('tbody.pair').forEach(function (tb) { tb.classList.toggle('sel', Number(tb.getAttribute('data-rank')) === sel); });
    }

    download(kind) {
      var r = this.state.result;
      if (!r) return;
      var stem = (r.input && r.input.filename ? r.input.filename.replace(/\.[^.]*$/, '') : 'degenprime') + (kind === 'alignment' ? '' : '_primers');
      var files = { csv: [r.csv, 'text/csv', '.csv'], txt: [r.report, 'text/plain', '_report.txt'], fasta: [r.fasta, 'text/plain', '.fasta'],
        json: [JSON.stringify(r, null, 2), 'application/json', '.json'], alignment: [r.alignment || '', 'text/plain', '_aligned.fasta'] };
      var f = files[kind];
      if (!f) return;
      var url = URL.createObjectURL(new Blob([f[0]], { type: f[1] + ';charset=utf-8' }));
      var a = document.createElement('a');
      a.href = url; a.download = stem + f[2]; a.style.display = 'none';
      document.body.appendChild(a); a.click(); a.remove();
      setTimeout(function () { URL.revokeObjectURL(url); }, 4000);
      this.toast('Downloading ' + stem + f[2]);
    }

    // ------------------------------------------------------------ check
    viewCheck() {
      var c = this.state.check, s = c.settings;
      return '<section class="check"><h1>Check a primer</h1><p class="lede">Melting temperature, GC content, degeneracy and dimer risk for primers you already have. IUPAC codes are welcome, and a partner primer adds the pair checks.</p>' +
        '<form class="check-form" data-form="check" novalidate>' +
        '<label class="seqfield"><span class="lab">Primer, 5′→3′</span><input name="primer" value="' + esc(c.primer) + '" autocomplete="off" spellcheck="false" autocapitalize="characters" placeholder="GCTACCTGAACACCATCCCG"></label>' +
        '<label class="seqfield"><span class="lab">Partner primer, 5′→3′ <span class="opt">(optional)</span></span><input name="partner" value="' + esc(c.partner) + '" autocomplete="off" spellcheck="false" autocapitalize="characters" placeholder="CCCTGGATGTTGTASTCGCC"></label>' +
        '<div class="cond">' + CHECK_FIELDS.map(function (f) { return fieldHTML(f, s[f[0]]); }).join('') + '</div>' +
        '<div class="run"><button type="submit" class="btn primary">Check primer</button><button type="button" class="btn ghost" data-action="check-example">Try an example pair</button></div>' +
        '<p class="form-error" data-error role="alert">' + esc(c.error) + '</p></form>' +
        '<div class="check-out" data-check-out>' + (c.result ? this.checkResultHTML(c.result) : '') + '</div></section>';
    }

    checkResultHTML(r) {
      return (r.warnings || []).map(function (w) { return '<p class="notice warn">' + esc(w) + '</p>'; }).join('') +
        oligoCard(r.partner ? 'Primer' : 'Your primer', r.primer) + (r.partner ? oligoCard('Partner primer', r.partner) : '') + (r.pair ? pairCard(r.pair) : '');
    }

    runCheck(form) {
      var self = this, c = this.state.check, err = this.$('[data-error]');
      c.primer = form.elements.primer.value;
      c.partner = form.elements.partner.value;
      var s = Object.assign({}, c.settings);
      CHECK_FIELDS.forEach(function (f) {
        var raw = String(form.elements[f[0]].value).trim(), n = Number(raw);
        if (raw !== '' && isFinite(n)) s[f[0]] = n;
      });
      c.settings = s;
      save('check', s);
      if (!c.primer.trim()) { err.textContent = 'Enter a primer sequence, 5′→3′.'; return Promise.resolve(); }
      err.textContent = '';
      return this.engine.call('check', { primer: c.primer, partner: c.partner.trim() ? c.partner : null, settings: s }).then(function (r) {
        var out = self.$('[data-check-out]'), e2 = self.$('[data-error]');
        if (!r.ok) { c.result = null; c.error = r.error || 'The check failed.'; if (e2) e2.textContent = c.error; if (out) out.innerHTML = ''; return; }
        c.result = r; c.error = '';
        if (out) out.innerHTML = self.checkResultHTML(r);
      }).catch(function (e) { var e2 = self.$('[data-error]'); if (e2) e2.textContent = e.message; });
    }

    // ------------------------------------------------------------ shared UI
    copy(text, label) {
      var self = this;
      var fallback = function () {
        var ta = document.createElement('textarea'), ok = false;
        ta.value = text; ta.setAttribute('readonly', ''); ta.style.position = 'fixed'; ta.style.opacity = '0';
        document.body.appendChild(ta); ta.select();
        try { ok = document.execCommand('copy'); } catch (e) { ok = false; }
        ta.remove();
        return ok;
      };
      var done = function (ok) { self.toast(ok ? 'Copied ' + label : 'Copy failed: select the text and copy it yourself.'); };
      if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(function () { done(true); }, function () { done(fallback()); });
      else done(fallback());
    }

    toast(msg) {
      var t = this.$('.toast');
      if (!t) return;
      t.textContent = msg;
      t.classList.add('show');
      clearTimeout(this._toastTimer);
      this._toastTimer = setTimeout(function () { t.classList.remove('show'); }, 2400);
    }

    bind() {
      var self = this, r = this.root;
      r.addEventListener('click', function (e) {
        var el = e.target.closest ? e.target.closest('[data-action]') : null;
        if (!el) return;
        var a = el.getAttribute('data-action');
        if (el.tagName === 'A') e.preventDefault();
        if (a === 'example') self.loadExample();
        else if (a === 'clear') self.setInput('', '');
        else if (a === 'cancel') self.cancelDesign();
        else if (a === 'reset') {
          self.state.design = Object.assign({}, self.defaults.design); save('design', self.state.design);
          self.state.align = Object.assign({}, ALIGN_DEFAULTS); save('align', self.state.align);
          self.render();
        }
        else if (a === 'copy') self.copy(el.getAttribute('data-text'), el.getAttribute('data-label') || 'sequence');
        else if (a === 'download') self.download(el.getAttribute('data-kind'));
        else if (a === 'select') { self.state.selected = Number(el.getAttribute('data-rank')); self.updateSelection(); }
        else if (a === 'jump') { var t = self.$('#' + el.getAttribute('data-target')); if (t) t.scrollIntoView({ block: 'start' }); }
        else if (a === 'check-example') {
          var f = self.$('form[data-form="check"]');
          if (!f) return;
          f.elements.primer.value = 'GCTACCTGAACACCATCCCG';
          f.elements.partner.value = 'CCCTGGATGTTGTASTCGCC';
          self.runCheck(f);
        }
      });
      r.addEventListener('submit', function (e) {
        var f = e.target;
        e.preventDefault();
        if (f.getAttribute('data-form') === 'design') self.runDesign(f);
        else if (f.getAttribute('data-form') === 'check') self.runCheck(f);
      });
      r.addEventListener('input', function (e) {
        var t = e.target;
        if (t.id === 'dgp-input') { if (t.readOnly) return; self.state.input = t.value; self.state.filename = ''; self.state.bigInput = false; var fn = self.$('[data-filename]'); if (fn) fn.textContent = ''; self.scheduleInspect(); }
        else if (t.form && t.form.getAttribute('data-form') === 'design' && t.name) { self.state.design = readDesignForm(t.form, self.state.design).settings; save('design', self.state.design); }
        else if (t.form && t.form.getAttribute('data-form') === 'check' && (t.name === 'primer' || t.name === 'partner')) self.state.check[t.name] = t.value;
      });
      r.addEventListener('change', function (e) {
        var t = e.target;
        if (t.matches && t.matches('input[data-file]') && t.files && t.files[0]) { self.readFile(t.files[0]); t.value = ''; }
        else if (t.name === 'rank_by' && t.form) { self.state.design = readDesignForm(t.form, self.state.design).settings; save('design', self.state.design); }
        else if ((t.name === 'align' || t.name === 'align_mode' || t.name === 'adjust_direction') && t.form) { self.readAlign(t.form); self.paintInspect(); }
      });
      r.addEventListener('dragover', function (e) { var d = e.target.closest && e.target.closest('[data-drop]'); if (d) { e.preventDefault(); d.classList.add('over'); } });
      r.addEventListener('dragleave', function (e) { var d = e.target.closest && e.target.closest('[data-drop]'); if (d) d.classList.remove('over'); });
      r.addEventListener('drop', function (e) {
        var d = e.target.closest && e.target.closest('[data-drop]');
        if (!d) return;
        e.preventDefault();
        d.classList.remove('over');
        var f = e.dataTransfer && e.dataTransfer.files && e.dataTransfer.files[0];
        if (f) self.readFile(f);
      });
    }
  }

  window.customElements.define('degenprime-app', DegenPrimeApp);
})();
