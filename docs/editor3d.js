/* Editor de contornos volumétrico: planos axial, coronal y sagital sobre un solo volumen de etiquetas.
   Herramientas: navegar, pincel, borrador, lazo (agregar o quitar), relleno por densidad (varita), máscara por
   rango de HU, interpolación entre cortes editados, deshacer. Las correcciones se guardan en el navegador y se
   descargan como cambios sobre la grilla del volumen ({formato: "contornos-editados-v2"}).
   EditorVolumen({ contenedor, base: 'datos/editor/', estructuras: [{e, n, c:[r,g,b], on}], fondo: 1, clave, pagina }) */
(function () {
  function h(tag, attrs, html) { var el = document.createElement(tag); for (var k in (attrs || {})) el.setAttribute(k, attrs[k]); if (html !== undefined) el.innerHTML = html; return el; }
  function gunzip(url) {
    return fetch(url).then(function (r) { if (!r.ok) throw new Error(url); return r.arrayBuffer(); }).then(function (buf) {
      if (typeof DecompressionStream === 'undefined') throw new Error('este navegador no descomprime gzip (DecompressionStream)');
      var ds = new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip')));
      return ds.arrayBuffer();
    }).then(function (b) { return new Uint8Array(b); });
  }
  // transformada de distancia euclídea 1D (Felzenszwalb y Huttenlocher), distancias al cuadrado
  function edt1d(f, n, d, v, z) {
    var k = 0; v[0] = 0; z[0] = -1e20; z[1] = 1e20;
    for (var q = 1; q < n; q++) {
      var s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]);
      while (s <= z[k]) { k--; s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2 * q - 2 * v[k]); }
      k++; v[k] = q; z[k] = s; z[k + 1] = 1e20;
    }
    k = 0;
    for (var q2 = 0; q2 < n; q2++) { while (z[k + 1] < q2) k++; d[q2] = (q2 - v[k]) * (q2 - v[k]) + f[v[k]]; }
  }
  function edt2d(mask, w, hh) {                       // distancia (vóxeles) al píxel más cercano con mask = 1
    var INF = 1e20, g = new Float64Array(w * hh), n = Math.max(w, hh), f = new Float64Array(n), d = new Float64Array(n), v = new Int32Array(n), z = new Float64Array(n + 1);
    for (var x = 0; x < w; x++) { for (var y = 0; y < hh; y++) f[y] = mask[y * w + x] ? 0 : INF; edt1d(f, hh, d, v, z); for (y = 0; y < hh; y++) g[y * w + x] = d[y]; }
    for (var y2 = 0; y2 < hh; y2++) { for (var x2 = 0; x2 < w; x2++) f[x2] = g[y2 * w + x2]; edt1d(f, w, d, v, z); for (x2 = 0; x2 < w; x2++) g[y2 * w + x2] = Math.sqrt(d[x2]); }
    return g;
  }

  function EditorVolumen(o) {
    var Z, Y, X, ct, lab, orig, vox = 2, nombres = [];
    var plano = 'ax', pos = { z: 0, y: 0, x: 0 }, herr = 'pincel', estr = o.estructuras[0].e, radio = 3, tol = 60;
    var mascara = false, mmin = -50, mmax = 150, win = { c: 40, w: 400 }, opacidad = 0.35;
    var pila = [], op = null, editados = { ax: {}, co: {}, sa: {} }, lazo = null, cursor = null, pintando = false;
    var raiz = h('div', { class: 'ev' });
    raiz.innerHTML =
      '<div class="fila ev-planos"><span class="lbl">Plano</span><button data-p="ax" class="primary">Axial</button><button data-p="co">Coronal</button><button data-p="sa">Sagital</button>' +
      '<span class="lbl">Corte</span><input type="range" class="ev-corte" min="0" max="1" value="0" aria-label="Corte"><span class="tag ev-pos"></span></div>' +
      '<div class="stage"><canvas class="ev-lienzo" width="900" height="900" tabindex="0" style="touch-action:none"></canvas><div class="hud ev-hud"></div></div>' +
      '<div class="fila"><span class="lbl">Herramienta</span>' +
      ['navegar:Navegar', 'pincel:Pincel', 'borrador:Borrador', 'lazo:Lazo (agregar)', 'lazo-quitar:Lazo (quitar)', 'relleno:Relleno'].map(function (q) { var p = q.split(':'); return '<button data-h="' + p[0] + '"' + (p[0] === 'pincel' ? ' class="primary"' : '') + '>' + p[1] + '</button>'; }).join('') + '</div>' +
      '<div class="fila"><span class="lbl">Estructura</span><select class="ev-estr" aria-label="Estructura"></select><span class="lbl">Tamaño</span><input type="range" class="ev-radio" min="1" max="15" value="3" style="flex:0 1 110px" aria-label="Tamaño del pincel">' +
      '<span class="lbl">Tolerancia</span><input type="range" class="ev-tol" min="10" max="300" value="60" style="flex:0 1 110px" aria-label="Tolerancia del relleno en HU"><span class="tag ev-tol-t">±60 HU</span>' +
      '<span class="lbl">Relleno</span><input type="range" class="ev-op" min="0" max="80" value="35" style="flex:0 1 90px" aria-label="Opacidad del relleno"></div>' +
      '<div class="fila"><label><input type="checkbox" class="ev-masc"> <b>Máscara por densidad</b></label><span class="lbl">de</span><input type="number" class="ev-mmin" value="-50" step="10" style="width:5.5em"><span class="lbl">a</span><input type="number" class="ev-mmax" value="150" step="10" style="width:5.5em"><span class="lbl">HU</span>' +
      '<select class="ev-pres" aria-label="Rangos predefinidos"><option value="">rangos…</option><option value="-50,150">tejido blando −50 a 150</option><option value="20,120">órgano sólido 20 a 120</option><option value="-10,25">líquido / orina −10 a 25</option><option value="200,1550">hueso &gt; 200</option><option value="-1000,-300">aire / pulmón</option><option value="-150,-30">grasa</option></select>' +
      '<span class="lbl">Ventana</span><select class="ev-win"><option value="40,400">partes blandas</option><option value="60,160">órganos</option><option value="400,1800">hueso</option><option value="-600,1500">pulmón</option></select></div>' +
      '<div class="fila"><button class="ev-interp">Interpolar entre cortes editados</button><button class="ev-deshacer">Deshacer</button><button class="ev-rest">Restaurar estructura en este corte</button>' +
      '<button class="ev-desc">Descargar correcciones</button><label><input type="file" class="ev-carg" accept=".json" hidden><span style="text-decoration:underline;cursor:pointer">Cargar correcciones</span></label><span class="tag ev-info">cargando volumen…</span></div>' +
      '<div class="fila ev-vis"></div>';
    o.contenedor.appendChild(raiz);
    var $ = function (s) { return raiz.querySelector(s); }, $$ = function (s) { return raiz.querySelectorAll(s); };
    var lienzo = $('.ev-lienzo'), ctx = lienzo.getContext('2d'), vista = { ox: 0, oy: 0, s: 1, w: 1, h: 1 };
    var visibles = {};
    o.estructuras.forEach(function (q) {
      visibles[q.e] = q.on !== false;
      var op_ = h('option', { value: q.e }, q.n); $('.ev-estr').appendChild(op_);
      var l = h('label', {}, '<input type="checkbox"' + (visibles[q.e] ? ' checked' : '') + '> <span style="display:inline-block;width:11px;height:11px;border-radius:2px;vertical-align:-1px;background:rgb(' + q.c.join(',') + ')"></span> ' + q.n);
      l.firstChild.addEventListener('change', function () { visibles[q.e] = this.checked; dibujar(); }); $('.ev-vis').appendChild(l);
    });
    var colores = {}; o.estructuras.forEach(function (q) { colores[q.e] = q.c; });

    // ---- geometría de planos ----
    function dims(p) { return p === 'ax' ? { w: X, h: Y, n: Z } : p === 'co' ? { w: X, h: Z, n: Y } : { w: Y, h: Z, n: X }; }
    function corte(p) { return p === 'ax' ? pos.z : p === 'co' ? pos.y : pos.x; }
    function indice(p, k, u, v) {                       // (u columna, v fila desde arriba) del corte k -> índice del volumen
      if (p === 'ax') return (k * Y + v) * X + u;
      if (p === 'co') return ((Z - 1 - v) * Y + k) * X + u;
      return ((Z - 1 - v) * Y + u) * X + k;
    }
    function hu(i) { return ct[i] * 10 - 1000; }
    function permitido(i) { if (orig[i] === 0) return false; if (mascara) { var v = hu(i); if (v < mmin || v > mmax) return false; } return true; }
    function fijar(i, val) {
      if (lab[i] === val) return;
      if (op && !op.vistos.has(i)) { op.vistos.add(i); op.idx.push(i); op.viejo.push(lab[i]); }
      lab[i] = val;
    }
    function quitar(i) { if (lab[i] === estr) fijar(i, orig[i] === estr ? o.fondo : orig[i]); }
    function empezarOp() { op = { idx: [], viejo: [], vistos: new Set() }; }
    function terminarOp(marcarCorte) {
      if (op && op.idx.length) { pila.push(op); if (pila.length > 60) pila.shift(); if (marcarCorte !== false) { (editados[plano][estr] = editados[plano][estr] || new Set()).add(corte(plano)); } guardar(); }
      op = null;
    }

    // ---- dibujo ----
    var lienzoCT = document.createElement('canvas'), lienzoEt = document.createElement('canvas');
    function dibujar() {
      if (!lab) return;
      var d = dims(plano), k = corte(plano), w = d.w, hh = d.h;
      var W = lienzo.width, H = lienzo.height, s = Math.min(W / w, H / hh);
      vista = { ox: (W - w * s) / 2, oy: (H - hh * s) / 2, s: s, w: w, h: hh };
      lienzoCT.width = w; lienzoCT.height = hh;
      var c1 = lienzoCT.getContext('2d'), im = c1.createImageData(w, hh), lo = win.c - win.w / 2;
      for (var v = 0; v < hh; v++) for (var u = 0; u < w; u++) {
        var i = indice(plano, k, u, v), g = (hu(i) - lo) / win.w * 255; g = g < 0 ? 0 : g > 255 ? 255 : g;
        var p = (v * w + u) * 4; im.data[p] = im.data[p + 1] = im.data[p + 2] = g; im.data[p + 3] = 255;
      }
      c1.putImageData(im, 0, 0);
      var E = 2; lienzoEt.width = w * E; lienzoEt.height = hh * E;
      var c2 = lienzoEt.getContext('2d'), im2 = c2.createImageData(w * E, hh * E), dd = im2.data;
      function et(uu, vv) { if (uu < 0 || vv < 0 || uu >= w || vv >= hh) return 0; return lab[indice(plano, k, uu, vv)]; }
      var a = Math.round(opacidad * 255);
      for (var vv = 0; vv < hh * E; vv++) for (var uu = 0; uu < w * E; uu++) {
        var U = uu >> 1, V = vv >> 1, e = et(U, V), col = colores[e]; if (!col || !visibles[e]) continue;
        var borde = (uu & 1) ? et(U + 1, V) !== e : et(U - 1, V) !== e; if (!borde) borde = (vv & 1) ? et(U, V + 1) !== e : et(U, V - 1) !== e;
        var q = (vv * w * E + uu) * 4;
        if (borde) { dd[q] = col[0]; dd[q + 1] = col[1]; dd[q + 2] = col[2]; dd[q + 3] = 255; }
        else if (e === estr && a > 0) { dd[q] = col[0]; dd[q + 1] = col[1]; dd[q + 2] = col[2]; dd[q + 3] = a; }
      }
      c2.putImageData(im2, 0, 0);
      ctx.fillStyle = '#000'; ctx.fillRect(0, 0, W, H);
      ctx.imageSmoothingEnabled = true; ctx.drawImage(lienzoCT, vista.ox, vista.oy, w * s, hh * s);
      ctx.imageSmoothingEnabled = false; ctx.drawImage(lienzoEt, vista.ox, vista.oy, w * s, hh * s); ctx.imageSmoothingEnabled = true;
      // líneas de los otros planos
      ctx.strokeStyle = 'rgba(255,255,255,0.25)'; ctx.lineWidth = 1; ctx.beginPath();
      var cu = plano === 'sa' ? pos.y : pos.x, cv = plano === 'ax' ? pos.y : (Z - 1 - pos.z);
      ctx.moveTo(vista.ox + (cu + 0.5) * s, vista.oy); ctx.lineTo(vista.ox + (cu + 0.5) * s, vista.oy + hh * s);
      ctx.moveTo(vista.ox, vista.oy + (cv + 0.5) * s); ctx.lineTo(vista.ox + w * s, vista.oy + (cv + 0.5) * s); ctx.stroke();
      if (lazo && lazo.length > 1) { ctx.strokeStyle = herr === 'lazo' ? '#ffffff' : '#ff6b6b'; ctx.lineWidth = 2; ctx.beginPath(); lazo.forEach(function (p, j) { var X_ = vista.ox + p[0] * s, Y_ = vista.oy + p[1] * s; if (j) ctx.lineTo(X_, Y_); else ctx.moveTo(X_, Y_); }); ctx.stroke(); }
      if (cursor && (herr === 'pincel' || herr === 'borrador')) { ctx.strokeStyle = herr === 'pincel' ? '#fff' : '#ff6b6b'; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(vista.ox + cursor[0] * s, vista.oy + cursor[1] * s, (radio - 0.5) * s, 0, 6.283); ctx.stroke(); }
      var nomP = { ax: 'axial', co: 'coronal', sa: 'sagital' }[plano], rot = plano === 'ax' ? 'D ← → I · anterior arriba' : plano === 'co' ? 'D ← → I · cabeza arriba' : 'anterior ← → posterior · cabeza arriba';
      $('.ev-hud').textContent = nomP + ' · corte ' + (k + 1) + ' de ' + d.n + ' · ' + rot;
      $('.ev-pos').textContent = (k * vox / 10).toFixed(1) + ' cm';
    }
    function aVox(ev) { var r = lienzo.getBoundingClientRect(), f = lienzo.width / r.width; return [((ev.clientX - r.left) * f - vista.ox) / vista.s, ((ev.clientY - r.top) * f - vista.oy) / vista.s]; }

    // ---- herramientas ----
    function pincel(pu, pv) {
      var d = dims(plano), k = corte(plano), rr = radio - 0.5;
      for (var v = Math.floor(pv - rr); v <= Math.ceil(pv + rr); v++) {
        if (v < 0 || v >= d.h) continue;
        for (var u = Math.floor(pu - rr); u <= Math.ceil(pu + rr); u++) {
          if (u < 0 || u >= d.w) continue;
          var bajo = (u === Math.floor(pu) && v === Math.floor(pv));
          if (!bajo && (u + 0.5 - pu) * (u + 0.5 - pu) + (v + 0.5 - pv) * (v + 0.5 - pv) > rr * rr) continue;
          var i = indice(plano, k, u, v);
          if (herr === 'pincel') { if (permitido(i)) fijar(i, estr); } else quitar(i);
        }
      }
    }
    function dentro(poly, x, y) { var c = false; for (var i = 0, j = poly.length - 1; i < poly.length; j = i++) { var xi = poly[i][0], yi = poly[i][1], xj = poly[j][0], yj = poly[j][1]; if (((yi > y) !== (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) c = !c; } return c; }
    function aplicarLazo() {
      if (!lazo || lazo.length < 3) return;
      var d = dims(plano), k = corte(plano), u0 = Infinity, u1 = -Infinity, v0 = Infinity, v1 = -Infinity;
      lazo.forEach(function (p) { u0 = Math.min(u0, p[0]); u1 = Math.max(u1, p[0]); v0 = Math.min(v0, p[1]); v1 = Math.max(v1, p[1]); });
      empezarOp();
      for (var v = Math.max(0, Math.floor(v0)); v <= Math.min(d.h - 1, Math.ceil(v1)); v++) for (var u = Math.max(0, Math.floor(u0)); u <= Math.min(d.w - 1, Math.ceil(u1)); u++) {
        if (!dentro(lazo, u + 0.5, v + 0.5)) continue;
        var i = indice(plano, k, u, v);
        if (herr === 'lazo') { if (permitido(i)) fijar(i, estr); } else quitar(i);
      }
      terminarOp();
    }
    function relleno(pu, pv) {
      var d = dims(plano), k = corte(plano), u = Math.floor(pu), v = Math.floor(pv);
      if (u < 0 || v < 0 || u >= d.w || v >= d.h) return;
      var semilla = hu(indice(plano, k, u, v)), visto = new Uint8Array(d.w * d.h), cola = [u, v], n = 0;
      visto[v * d.w + u] = 1; empezarOp();
      while (cola.length && n < 60000) {
        var vv = cola.pop(), uu = cola.pop(), i = indice(plano, k, uu, vv);
        if (!permitido(i) || Math.abs(hu(i) - semilla) > tol) continue;
        fijar(i, estr); n++;
        [[1, 0], [-1, 0], [0, 1], [0, -1]].forEach(function (q) { var a = uu + q[0], b = vv + q[1]; if (a >= 0 && b >= 0 && a < d.w && b < d.h && !visto[b * d.w + a]) { visto[b * d.w + a] = 1; cola.push(a, b); } });
      }
      terminarOp();
      info(n >= 60000 ? 'relleno cortado en 60 000 vóxeles: sube la máscara o baja la tolerancia' : null);
    }
    function interpolar() {
      var marcados = editados[plano][estr] ? Array.from(editados[plano][estr]).sort(function (a, b) { return a - b; }) : [];
      if (marcados.length < 2) { info('edita al menos dos cortes de esta estructura en este plano para interpolar'); return; }
      var d = dims(plano), n = d.w * d.h, hecho = 0;
      function mascaraCorte(k) { var m = new Uint8Array(n); for (var v = 0; v < d.h; v++) for (var u = 0; u < d.w; u++) m[v * d.w + u] = lab[indice(plano, k, u, v)] === estr ? 1 : 0; return m; }
      function sdf(m) { var fuera = edt2d(m, d.w, d.h), inv = new Uint8Array(n); for (var i = 0; i < n; i++) inv[i] = m[i] ? 0 : 1; var dentro_ = edt2d(inv, d.w, d.h), s = new Float64Array(n); for (i = 0; i < n; i++) s[i] = m[i] ? -dentro_[i] + 0.5 : fuera[i] - 0.5; return s; }
      empezarOp();
      for (var j = 0; j + 1 < marcados.length; j++) {
        var a = marcados[j], b = marcados[j + 1]; if (b - a < 2) continue;
        var sa = sdf(mascaraCorte(a)), sb = sdf(mascaraCorte(b));
        for (var k = a + 1; k < b; k++) {
          var t = (k - a) / (b - a);
          for (var v = 0; v < d.h; v++) for (var u = 0; u < d.w; u++) {
            var q = v * d.w + u, val = (1 - t) * sa[q] + t * sb[q], i = indice(plano, k, u, v);
            if (val <= 0) { if (permitido(i)) fijar(i, estr); } else if (lab[i] === estr) quitar(i);
          }
          hecho++;
        }
      }
      terminarOp(false);
      info(hecho ? hecho + ' cortes interpolados en el plano ' + { ax: 'axial', co: 'coronal', sa: 'sagital' }[plano] : 'los cortes editados son contiguos: no hay nada que interpolar');
      dibujar();
    }

    // ---- eventos ----
    lienzo.addEventListener('pointerdown', function (ev) {
      if (!lab) return; ev.preventDefault(); try { lienzo.focus({ preventScroll: true }); lienzo.setPointerCapture(ev.pointerId); } catch (e) {}
      var p = aVox(ev);
      if (herr === 'navegar') { navegar(p); return; }
      if (herr === 'relleno') { relleno(p[0], p[1]); dibujar(); return; }
      if (herr === 'lazo' || herr === 'lazo-quitar') { lazo = [p]; pintando = true; return; }
      empezarOp(); pintando = true; pincel(p[0], p[1]); dibujar();
    });
    lienzo.addEventListener('pointermove', function (ev) {
      if (!lab) return; var p = aVox(ev); cursor = p;
      if (pintando && lazo) lazo.push(p);
      else if (pintando && herr === 'navegar') navegar(p);
      else if (pintando) pincel(p[0], p[1]);
      dibujar();
    });
    function fin() { if (!pintando) return; pintando = false; if (lazo) { aplicarLazo(); lazo = null; } else terminarOp(); dibujar(); }
    lienzo.addEventListener('pointerup', fin); lienzo.addEventListener('pointercancel', fin);
    lienzo.addEventListener('pointerleave', function () { cursor = null; dibujar(); });
    function navegar(p) {
      var u = Math.max(0, Math.min(dims(plano).w - 1, Math.floor(p[0]))), v = Math.max(0, Math.min(dims(plano).h - 1, Math.floor(p[1])));
      if (plano === 'ax') { pos.x = u; pos.y = v; } else if (plano === 'co') { pos.x = u; pos.z = Z - 1 - v; } else { pos.y = u; pos.z = Z - 1 - v; }
      pintando = true; dibujar();
    }
    lienzo.addEventListener('keydown', function (ev) {
      if (ev.key === 'ArrowUp' || ev.key === 'ArrowDown' || ev.key === 'PageUp' || ev.key === 'PageDown') { ev.preventDefault(); mover(ev.key === 'ArrowUp' || ev.key === 'PageUp' ? 1 : -1); }
      if ((ev.ctrlKey || ev.metaKey) && ev.key.toLowerCase() === 'z') { ev.preventDefault(); deshacer(); }
    });
    lienzo.addEventListener('wheel', function (ev) { if (!lab) return; ev.preventDefault(); mover(ev.deltaY < 0 ? 1 : -1); }, { passive: false });
    function mover(dk) { var d = dims(plano), k = Math.max(0, Math.min(d.n - 1, corte(plano) + dk)); fijarCorte(k); }
    function fijarCorte(k) { if (plano === 'ax') pos.z = k; else if (plano === 'co') pos.y = k; else pos.x = k; $('.ev-corte').value = k; dibujar(); }
    $('.ev-corte').addEventListener('input', function () { fijarCorte(parseInt(this.value, 10)); });
    Array.prototype.forEach.call($$('.ev-planos button'), function (b) { b.addEventListener('click', function () {
      plano = b.getAttribute('data-p'); Array.prototype.forEach.call($$('.ev-planos button'), function (x) { x.className = x === b ? 'primary' : ''; });
      var d = dims(plano); $('.ev-corte').max = d.n - 1; $('.ev-corte').value = corte(plano); dibujar(); }); });
    Array.prototype.forEach.call(raiz.querySelectorAll('button[data-h]'), function (b) { b.addEventListener('click', function () {
      herr = b.getAttribute('data-h'); Array.prototype.forEach.call(raiz.querySelectorAll('button[data-h]'), function (x) { x.className = x === b ? 'primary' : ''; });
      lienzo.style.cursor = herr === 'navegar' ? 'move' : 'crosshair'; dibujar(); }); });
    $('.ev-estr').addEventListener('change', function () { estr = +this.value; dibujar(); });
    $('.ev-radio').addEventListener('input', function () { radio = +this.value; dibujar(); });
    $('.ev-tol').addEventListener('input', function () { tol = +this.value; $('.ev-tol-t').textContent = '±' + tol + ' HU'; });
    $('.ev-op').addEventListener('input', function () { opacidad = +this.value / 100; dibujar(); });
    $('.ev-masc').addEventListener('change', function () { mascara = this.checked; });
    $('.ev-mmin').addEventListener('change', function () { mmin = +this.value; });
    $('.ev-mmax').addEventListener('change', function () { mmax = +this.value; });
    $('.ev-pres').addEventListener('change', function () { if (!this.value) return; var p = this.value.split(','); mmin = +p[0]; mmax = +p[1]; $('.ev-mmin').value = mmin; $('.ev-mmax').value = mmax; mascara = true; $('.ev-masc').checked = true; this.value = ''; });
    $('.ev-win').addEventListener('change', function () { var p = this.value.split(','); win = { c: +p[0], w: +p[1] }; dibujar(); });
    function deshacer() { var u = pila.pop(); if (!u) return; for (var j = 0; j < u.idx.length; j++) lab[u.idx[j]] = u.viejo[j]; guardar(); dibujar(); }
    $('.ev-deshacer').addEventListener('click', deshacer);
    $('.ev-interp').addEventListener('click', interpolar);
    $('.ev-rest').addEventListener('click', function () {
      var d = dims(plano), k = corte(plano); empezarOp();
      for (var v = 0; v < d.h; v++) for (var u = 0; u < d.w; u++) { var i = indice(plano, k, u, v); if (lab[i] === estr || orig[i] === estr) fijar(i, orig[i]); }
      terminarOp(false); dibujar();
    });

    // ---- guardar, descargar, cargar ----
    function cambios() { var out = [], n = lab.length, i = 0; while (i < n) { if (lab[i] !== orig[i]) { var v = lab[i], j = i; while (j < n && lab[j] !== orig[j] && lab[j] === v) j++; out.push([i, j - i, v]); i = j; } else i++; } return out; }
    function exportar() { return { formato: 'contornos-editados-v2', fecha: new Date().toISOString(), meta: { pagina: o.pagina, forma_zyx: [Z, Y, X], voxel_mm: vox, nombres: nombres }, cambios: cambios() }; }
    function info(msg) {
      if (!lab) return; var c = cambios(), n = 0, vol = 0; c.forEach(function (r) { n += r[1]; });
      for (var i = 0; i < lab.length; i++) if (lab[i] === estr) vol++;
      $('.ev-info').textContent = (msg ? msg + ' · ' : '') + (n ? n + ' vóxeles corregidos' : 'sin cambios') + ' · ' + $('.ev-estr').selectedOptions[0].textContent + ': ' + (vol * vox * vox * vox / 1000).toFixed(1) + ' mL';
    }
    function guardar() { try { localStorage.setItem('contornos3d:' + o.clave, JSON.stringify(exportar())); } catch (e) {} info(); }
    function aplicar(dat) {
      if (dat.formato !== 'contornos-editados-v2' || !dat.meta || dat.meta.forma_zyx.join() !== [Z, Y, X].join()) throw new Error('forma');
      lab.set(orig); dat.cambios.forEach(function (r) { lab.fill(r[2], r[0], r[0] + r[1]); }); guardar(); dibujar();
    }
    $('.ev-desc').addEventListener('click', function () {
      var blob = new Blob([JSON.stringify(exportar())], { type: 'application/json' }), a = document.createElement('a');
      a.href = URL.createObjectURL(blob); a.download = 'contornos-' + o.clave + '.json'; document.body.appendChild(a); a.click();
      setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    });
    $('.ev-carg').addEventListener('change', function () {
      var f = this.files[0]; if (!f) return; var lr = new FileReader();
      lr.onload = function () { try { aplicar(JSON.parse(lr.result)); } catch (e) { $('.ev-info').textContent = 'el archivo no corresponde a este volumen'; } };
      lr.readAsText(f); this.value = '';
    });

    // ---- carga ----
    Promise.all([fetch(o.base + 'editor.json').then(function (r) { return r.json(); }), gunzip(o.base + 'ct.u8.gz'), gunzip(o.base + 'etiquetas.u8.gz')]).then(function (r) {
      Z = r[0].forma_zyx[0]; Y = r[0].forma_zyx[1]; X = r[0].forma_zyx[2]; vox = r[0].voxel_mm; nombres = r[0].nombres;
      ct = r[1]; lab = r[2]; orig = lab.slice();
      var k0 = o.corteInicial ? o.corteInicial(lab, Z, Y, X) : Math.floor(Z / 2);
      pos = { z: k0, y: Math.floor(Y / 2), x: Math.floor(X / 2) };
      $('.ev-corte').max = Z - 1; $('.ev-corte').value = pos.z;
      try { var g = localStorage.getItem('contornos3d:' + o.clave); if (g) aplicar(JSON.parse(g)); } catch (e) {}
      info(); dibujar();
    }).catch(function (e) { $('.ev-info').textContent = 'no se pudo cargar el volumen: ' + e.message; });
    return { _estado: function () { return { Z: Z, Y: Y, X: X, lab: lab, orig: orig, pos: pos, plano: plano }; } };
  }
  window.EditorVolumen = EditorVolumen;
})();
