/* Editor de contornos por corte (pincel, borrador, deshacer, copiar del corte vecino, descargar/cargar).
   Uso: var ed = EditorContornos({ ... }); la página le entrega:
     contenedor   elemento donde se arma la barra de herramientas
     lienzo       canvas sobre el que se dibuja (eventos de puntero)
     estructuras  [{e: etiqueta, n: nombre, c: [r,g,b]}]   (las que se pueden pintar)
     fondo        etiqueta que se escribe al borrar (tejido blando)
     corte()      índice del corte visible
     etiquetas(k) Promise<{w, h, e: Uint8Array}> del corte k (el editor modifica ese mismo arreglo)
     aEtiqueta(x, y)  coordenadas del lienzo (px CSS) -> [columna, fila] en la imagen de etiquetas (flotantes)
     escalaPx()   cuántos px CSS mide un píxel de etiqueta (para dibujar el cursor)
     redibujar()  vuelve a dibujar el corte
     clave        nombre para guardar en el navegador (página + caso)
     meta         objeto que se copia en el archivo descargado (grilla, caso, etc.)
   El archivo descargado guarda solo los píxeles cambiados: {formato, meta, cortes: {k: [[inicio, largo, etiqueta], ...]}} */
(function () {
  function EditorContornos(o) {
    var activo = false, estr = o.estructuras[0].e, radio = 3, modo = 'pincel', pintando = false;
    var originales = {}, deshacer = {}, cursor = null;
    var bar = document.createElement('div');
    bar.className = 'fila editor-barra';
    bar.innerHTML =
      '<label><input type="checkbox" class="ed-on"> <b>Editar contornos</b></label>' +
      '<span class="ed-herr" style="display:none">' +
      '<select class="ed-estr" aria-label="Estructura"></select>' +
      '<button class="ed-pincel primary" type="button">Pincel</button><button class="ed-borrar" type="button">Borrador</button>' +
      '<span class="lbl">Tamaño</span><input type="range" class="ed-radio" min="1" max="10" value="3" style="flex:0 1 110px" aria-label="Tamaño del pincel">' +
      '<button class="ed-deshacer" type="button">Deshacer</button><button class="ed-ant" type="button">Copiar del anterior</button><button class="ed-sig" type="button">Copiar del siguiente</button>' +
      '<button class="ed-rest" type="button">Restaurar corte</button>' +
      '<button class="ed-desc" type="button">Descargar correcciones</button><label class="ed-carg-l"><input type="file" class="ed-carg" accept=".json" hidden><span style="text-decoration:underline;cursor:pointer">Cargar correcciones</span></label>' +
      '<span class="tag ed-info"></span></span>';
    o.contenedor.appendChild(bar);
    var $ = function (s) { return bar.querySelector(s); };
    o.estructuras.forEach(function (q) { var op = document.createElement('option'); op.value = q.e; op.textContent = q.n; $('.ed-estr').appendChild(op); });
    function info() {
      var n = Object.keys(originales).filter(function (k) { return cambios(+k).length; }).length;
      $('.ed-info').textContent = n ? n + (n === 1 ? ' corte corregido' : ' cortes corregidos') : 'sin cambios';
    }
    function guardar() {
      try { localStorage.setItem('contornos:' + o.clave, JSON.stringify(exportar())); } catch (e) {}
      info();
    }
    function cambios(k) {
      var L = o._cache && o._cache[k], ori = originales[k]; if (!L || !ori) return [];
      var out = [], i = 0, n = ori.length;
      while (i < n) {
        if (L.e[i] !== ori[i]) { var v = L.e[i], j = i; while (j < n && L.e[j] !== ori[j] && L.e[j] === v) j++; out.push([i, j - i, v]); i = j; } else i++;
      }
      return out;
    }
    function exportar() {
      var cortes = {};
      Object.keys(originales).forEach(function (k) { var c = cambios(+k); if (c.length) cortes[k] = c; });
      return { formato: 'contornos-editados-v1', fecha: new Date().toISOString(), meta: o.meta, cortes: cortes };
    }
    o._cache = {};
    function conCorte(k) {
      return o.etiquetas(k).then(function (L) {
        o._cache[k] = L;
        if (!originales[k]) originales[k] = L.e.slice();
        return L;
      });
    }
    function aplicar(datos) {
      var ks = Object.keys(datos.cortes || {});
      return Promise.all(ks.map(function (k) {
        return conCorte(+k).then(function (L) { datos.cortes[k].forEach(function (r) { L.e.fill(r[2], r[0], r[0] + r[1]); }); });
      })).then(function () { info(); o.redibujar(); });
    }
    try { var g = localStorage.getItem('contornos:' + o.clave); if (g) aplicar(JSON.parse(g)); } catch (e) {}
    function pintar(ev) {
      var r = o.lienzo.getBoundingClientRect(), p = o.aEtiqueta(ev.clientX - r.left, ev.clientY - r.top), k = o.corte(), L = o._cache[k];
      cursor = [ev.clientX - r.left, ev.clientY - r.top];
      if (!L || !pintando) { o.redibujar(); return; }
      var cx = p[0], cy = p[1], rr = radio - 0.5;
      for (var y = Math.floor(cy - rr); y <= Math.ceil(cy + rr); y++) {
        if (y < 0 || y >= L.h) continue;
        for (var x = Math.floor(cx - rr); x <= Math.ceil(cx + rr); x++) {
          if (x < 0 || x >= L.w) continue;
          var bajo = (x === Math.floor(cx) && y === Math.floor(cy));   // el píxel bajo el puntero siempre entra
          if (!bajo && (x + 0.5 - cx) * (x + 0.5 - cx) + (y + 0.5 - cy) * (y + 0.5 - cy) > rr * rr) continue;
          var i = y * L.w + x, ori = originales[k][i];
          if (ori === 0) continue;                                   // nunca fuera del cuerpo
          if (modo === 'pincel') L.e[i] = estr;
          else if (L.e[i] === estr) L.e[i] = (ori === estr ? o.fondo : ori);   // borrar devuelve lo que había debajo
        }
      }
      o.redibujar();
    }
    o.lienzo.addEventListener('pointerdown', function (ev) {
      if (!activo) return; ev.preventDefault(); o.lienzo.setPointerCapture(ev.pointerId);
      var k = o.corte();
      function empezar(L) { (deshacer[k] = deshacer[k] || []).push(L.e.slice()); if (deshacer[k].length > 30) deshacer[k].shift(); pintando = true; pintar(ev); }
      if (o._cache[k] && originales[k]) empezar(o._cache[k]);          // corte ya cargado: sin demora
      else conCorte(k).then(empezar);
    });
    o.lienzo.addEventListener('pointermove', function (ev) { if (activo) pintar(ev); });
    ['pointerup', 'pointercancel'].forEach(function (t) { o.lienzo.addEventListener(t, function () { if (pintando) { pintando = false; guardar(); } }); });
    o.lienzo.addEventListener('pointerleave', function () { cursor = null; if (activo) o.redibujar(); });
    $('.ed-on').addEventListener('change', function () {
      activo = this.checked; $('.ed-herr').style.display = activo ? 'contents' : 'none'; o.lienzo.style.cursor = activo ? 'crosshair' : '';
      o.lienzo.style.touchAction = activo ? 'none' : ''; if (activo) conCorte(o.corte()); info(); o.redibujar();
    });
    $('.ed-estr').addEventListener('change', function () { estr = +this.value; });
    $('.ed-radio').addEventListener('input', function () { radio = +this.value; });
    $('.ed-pincel').addEventListener('click', function () { modo = 'pincel'; this.className = 'ed-pincel primary'; $('.ed-borrar').className = 'ed-borrar'; });
    $('.ed-borrar').addEventListener('click', function () { modo = 'borrar'; this.className = 'ed-borrar primary'; $('.ed-pincel').className = 'ed-pincel'; });
    $('.ed-deshacer').addEventListener('click', function () { var k = o.corte(), s = deshacer[k]; if (s && s.length) { o._cache[k].e.set(s.pop()); guardar(); o.redibujar(); } });
    $('.ed-rest').addEventListener('click', function () { var k = o.corte(); if (originales[k]) { (deshacer[k] = deshacer[k] || []).push(o._cache[k].e.slice()); o._cache[k].e.set(originales[k]); guardar(); o.redibujar(); } });
    function copiar(desde) {
      var k = o.corte(), j = k + desde;
      Promise.all([conCorte(k), conCorte(j)]).then(function (r) {
        var A = r[0], B = r[1]; if (A.w !== B.w || A.h !== B.h) return;
        (deshacer[k] = deshacer[k] || []).push(A.e.slice());
        for (var i = 0; i < A.e.length; i++) {
          if (originales[k][i] === 0) continue;
          var vb = B.e[i], va = A.e[i];
          if (vb === estr) A.e[i] = estr;                            // agrega la estructura del vecino
          else if (va === estr) A.e[i] = (originales[k][i] === estr ? o.fondo : originales[k][i]);   // quita lo que no está en el vecino
        }
        guardar(); o.redibujar();
      }).catch(function () {});
    }
    $('.ed-ant').addEventListener('click', function () { copiar(-1); });
    $('.ed-sig').addEventListener('click', function () { copiar(1); });
    $('.ed-desc').addEventListener('click', function () {
      var blob = new Blob([JSON.stringify(exportar())], { type: 'application/json' }), a = document.createElement('a');
      a.href = URL.createObjectURL(blob); a.download = 'contornos-' + o.clave.replace(/[^a-z0-9-]+/gi, '_') + '.json';
      document.body.appendChild(a); a.click(); setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
    });
    $('.ed-carg').addEventListener('change', function () {
      var f = this.files[0]; if (!f) return; var lr = new FileReader();
      lr.onload = function () { try { var d = JSON.parse(lr.result); if (d.formato !== 'contornos-editados-v1') throw 0; aplicar(d).then(guardar); } catch (e) { $('.ed-info').textContent = 'el archivo no es de correcciones de contornos'; } };
      lr.readAsText(f); this.value = '';
    });
    return {
      activo: function () { return activo; },
      estructura: function () { return estr; },
      cursor: function () { return cursor; },
      radio: function () { return radio; },
      // el editor necesita el arreglo vivo: la página debe pedir las etiquetas a través de esto
      etiquetas: function (k) { return o._cache[k] ? Promise.resolve(o._cache[k]) : conCorte(k); },
      dibujarCursor: function (ctx) {
        if (!activo || !cursor) return; var s = o.escalaPx(); ctx.save(); ctx.strokeStyle = modo === 'pincel' ? '#fff' : '#ff6b6b'; ctx.lineWidth = 1.5;
        var dpr = ctx.canvas.width / o.lienzo.getBoundingClientRect().width;
        ctx.beginPath(); ctx.arc(cursor[0] * dpr, cursor[1] * dpr, radio * s * dpr, 0, 6.283); ctx.stroke(); ctx.restore();
      }
    };
  }
  window.EditorContornos = EditorContornos;
})();
