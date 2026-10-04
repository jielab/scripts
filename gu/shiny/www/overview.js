/* Shared genome viewport. Coordinates sent to R are 0-based, half-open. */
window.guDensityBind = function (el, region) {
  // Wait for the second click before a Shiny redraw replaces the clicked bin.
  // onRender can run again on the same element; retain only our own handlers.
  if (el._guDensityHandlers && el.removeListener) {
    Object.entries(el._guDensityHandlers).forEach(([name, fn]) => el.removeListener(name, fn));
  }
  clearTimeout(el._guClickTimer);
  clearTimeout(el._guRangeTimer);
  // Plotly also emits relayout while initializing and resizing subplots.
  // Only a user gesture may change the shared genome viewport.
  if (el._guDensityGesture) {
    ['pointerdown', 'wheel', 'keydown'].forEach(name =>
      el.removeEventListener(name, el._guDensityGesture, true));
  }
  let userNavigation = false;
  el._guDensityGesture = function () { userNavigation = true; };
  ['pointerdown', 'wheel', 'keydown'].forEach(name =>
    el.addEventListener(name, el._guDensityGesture, true));
  let pending = null, openedAt = 0;
  function open() {
    if (!pending || Date.now() - openedAt < 450) return;
    clearTimeout(el._guClickTimer);
    Shiny.setInputValue('density_bin_open', region || pending.key, {priority: 'event'});
    pending = null; openedAt = Date.now();
  }
  const handlers = {plotly_click: function (e) {
    userNavigation = false;
    const p = e.points && e.points[0];
    if (!p || p.customdata == null || Date.now() - openedAt < 450) return;
    const key = p.customdata, now = Date.now();
    if (pending && pending.key === key && now - pending.time < 450) { open(); return; }
    clearTimeout(el._guClickTimer);
    pending = {key, time: now};
    el._guClickTimer = setTimeout(function () {
      pending = null;
      Shiny.setInputValue('density_bin_click', key, {priority: 'event'});
    }, 450);
  }, plotly_doubleclick: open, plotly_relayout: function (e) {
    if (!userNavigation) return;
    userNavigation = false;
    if (Object.keys(e).some(k => /^xaxis\d*\.autorange$/.test(k) && e[k])) {
      Shiny.setInputValue('density_range', {reset: true}, {priority: 'event'}); return;
    }
    const k = Object.keys(e).find(k => /^xaxis\d*\.range\[0\]$/.test(k));
    const a = Object.keys(e).find(k => /^xaxis\d*\.range$/.test(k));
    if (!k && !a) return;
    const range = a ? e[a] : [e[k], e[k.replace('[0]', '[1]')]];
    clearTimeout(el._guRangeTimer);
    el._guRangeTimer = setTimeout(() => Shiny.setInputValue('density_range', {start:range[0], end:range[1]}, {priority:'event'}), 200);
  }};
  Object.entries(handlers).forEach(([name, fn]) => el.on(name, fn));
  el._guDensityHandlers = handlers;
};
$(document).on('shiny:connected', function () {
  let browser, genome, applying = false, timer, pending, running = false, userNavigation = false;
  // Initial layout and server-driven IGV searches must not select a density region.
  ['pointerdown', 'wheel', 'keydown'].forEach(event => document.addEventListener(event, e => {
    if (!applying && e.target.closest && e.target.closest('#gu_igv')) userNavigation = true;
  }, true));
  function genomeConfig(id) {
    const db = id === 'chm13v2.0' ? 'hs1' : id;
    const base = 'https://hgdownload.soe.ucsc.edu/goldenPath/' + db;
    const local = 'gu_assets/genomes/hg19/';
    return {
      db,
      reference: {
        id: db, name: db, twoBitURL: base + '/bigZips/' + db + '.2bit',
        chromSizesURL: db === 'hg19' ? local + 'hg19.chrom.sizes' : base + '/bigZips/' + db + (db === 'hs1' ? '.chrom.sizes.txt' : '.chrom.sizes'),
        ...(db === 'hs1'
          ? {cytobandBbURL: 'https://hgdownload.soe.ucsc.edu/gbdb/hs1/cytoBandMapped/cytoBandMapped.bb'}
          : {cytobandURL: db === 'hg19' ? local + 'cytoBand.txt.gz' : base + '/database/cytoBandIdeo.txt.gz'})
      },
      genes: db === 'hg19'
        ? {format: 'refgene', url: local + 'ncbiRefSeq.txt.gz', indexed: false}
        : db === 'hs1'
          ? {format: 'biggenepred', url: 'https://hgdownload.soe.ucsc.edu/gbdb/hs1/ncbiRefSeq/ncbiRefSeqSelectCurated.bb'}
          : {format: 'refgene', url: base + '/database/ncbiRefSeq.txt.gz', indexed: false}
    };
  }
  async function loadGenes(owner, config) {
    // An annotation failure must not discard a working reference browser.
    try {
      let url = config.url;
      if (config.indexed === false) {
        // Fetch first so an unavailable file never creates a broken IGV track
        // or its blocking error dialog. IGV decompresses this named gzip File.
        const response = await fetch(url);
        if (!response.ok) throw new Error('Annotation HTTP ' + response.status);
        url = new File([await response.blob()], 'ncbiRefSeq.txt.gz');
      }
      if (owner !== browser) return;
      await owner.loadTrack({name: 'RefSeq genes', ...config, url, order: 1000});
    } catch (err) {
      console.warn('GU gene annotation:', err);
      if (owner !== browser) return;
      const status = document.getElementById('gu_igv_status');
      if (status) status.textContent = '基因注释暂不可用；仍可浏览参考序列和基因组位置。';
    }
  }
  async function update() {
    if (running) return;
    running = true;
    while (pending) {
      const r = pending; pending = null;
      const host = document.getElementById('gu_igv');
      if (!host || !window.igv) { running = false; pending = r; setTimeout(update, 300); return; }
      applying = true; userNavigation = false; clearTimeout(timer);
      host.style.pointerEvents = "none"; host.setAttribute("aria-busy", "true");
      try {
        if (!browser || genome !== r.genome) {
          if (browser) igv.removeBrowser(browser);
          browser = null;
          host.textContent = '';
          const config = genomeConfig(r.genome);
          if (r.reference && r.reference.fastaURL) {
            Object.assign(config.reference, r.reference);
            delete config.reference.twoBitURL;
            delete config.reference.chromSizesURL;
          }
          browser = await igv.createBrowser(host, {
            loadDefaultGenomes:false,
            reference:config.reference,
            locus:r.locus,showNavigation:true,showIdeogram:true,
            tracks:[]
          });
          genome = r.genome;
          const status = document.getElementById("gu_igv_status");
          if(status) status.textContent = "";
          browser.on('locuschange', function () {
            if (applying || !userNavigation) return;
            clearTimeout(timer);
            timer = setTimeout(function () {
              if (applying || !userNavigation) return;
              const loci = browser.currentLoci();
              const locus = Array.isArray(loci) ? loci[0] : loci;
              const m = /^([^:]+):([\d,]+)-([\d,]+)$/.exec(locus);
              if (m) {
                userNavigation = false;
                Shiny.setInputValue('igv_region', {chr:m[1],start:Number(m[2].replaceAll(',',''))-1,end:Number(m[3].replaceAll(',',''))}, {priority:'event'});
              }
            }, 250);
          });
          loadGenes(browser, config.genes);
        } else {
          const loci = browser.currentLoci();
          const current = (Array.isArray(loci) ? loci[0] : loci).replaceAll(',','');
          if (current !== r.locus) await browser.search(r.locus);
        }
      } catch (err) {
        console.error('GU IGV:', err);
        const status = document.getElementById('gu_igv_status'); if(status) status.textContent = ''; 
        if (browser) igv.removeBrowser(browser);
        host.textContent = 'IGV 参考序列加载失败；请使用 Open IGV / Open UCSC 链接。';
        browser = null;
      } finally { applying = false; host.style.pointerEvents = ""; host.setAttribute("aria-busy", "false"); }
    }
    running = false;
  }
  Shiny.addCustomMessageHandler('gu-igv-region', function (r) { pending = r; update(); });
  $(document).on('shown.bs.tab', function () {
    // Hidden panels have no measurable width. Refit headers and canvases once
    // the Overview / PhyML detail tab becomes visible again.
    setTimeout(function () {
      if ($.fn && $.fn.dataTable) $.fn.dataTable.tables({visible: true, api: true}).columns.adjust();
      const host = document.getElementById('gu_igv');
      if (host && host.offsetParent !== null && window.igv && igv.visibilityChange) igv.visibilityChange();
    }, 0);
  });
});
