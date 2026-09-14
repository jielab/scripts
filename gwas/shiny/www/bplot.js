/* Shared viewport events. Programmatic IGV searches must not echo into Shiny. */
(function () {
  "use strict";
  let browser = null, build = null, pending = null, working = false, applying = false;
  let blocked = false, view = null, igvTimer = null;
  const status = text => {
    const el = document.getElementById("igv-status");
    el.textContent = text;
    el.style.display = text ? "block" : "none";
  };
  const send = (name, value) => window.Shiny.setInputValue(name, value, {priority: "event"});
  const delay = ms => new Promise(resolve => setTimeout(resolve, ms));
  const esc = s => String(s).replace(/,/g, "");
  const lociOf = instance => {
    const loci = instance.currentLoci();
    return Array.isArray(loci) ? loci : loci ? [loci] : [];
  };
  function removeIgv(instance) {
    // IGV 3.0.0 leaves this observer alive during dispose. Disconnect before
    // removing tracks so it cannot inspect disposed tracks on a build switch.
    instance.roiManager?.observer?.disconnect();
    igv.removeBrowser(instance);
  }

  async function updateIgv() {
    if (working || !pending) return;
    working = true;
    while (pending) {
      const cfg = pending; pending = null;
      if (cfg.error) {
        status(cfg.error);
        if (browser) { removeIgv(browser); browser = null; build = null; }
        document.getElementById("igv").classList.remove("busy");
        continue;
      }
      applying = true;
      status("");
      try {
        if (!browser || build !== cfg.build) {
          if (browser) removeIgv(browser);
          browser = null; build = null;
          const instance = await igv.createBrowser(document.getElementById("igv"), {
            loadDefaultGenomes: false, reference: cfg.reference, locus: cfg.locus,
            showNavigation: true, showIdeogram: true, showCenterGuide: true,
            showCursorTrackingGuide: false, minimumBases: 40, tracks: []
          });
          browser = instance; build = cfg.build;
          instance.on("locuschange", function () {
            if (instance !== browser || applying || blocked) return;
            clearTimeout(igvTimer);
            igvTimer = setTimeout(function () {
              if (instance !== browser || applying || blocked) return;
              const loci = lociOf(instance);
              if (!loci || loci.length !== 1) return;
              if (/^all(?::|$)/i.test(loci[0])) { send("igv_event", {build, chr: "All"}); return; }
              const match = esc(loci[0]).match(/^(?:chr)?([0-9]+|X):([0-9]+(?:\.[0-9]+)?)-([0-9]+(?:\.[0-9]+)?)$/i);
              if (!match) return;
              const next = {build, chr: match[1].toUpperCase(), start: Math.round(Number(match[2])), end: Math.round(Number(match[3]))};
              if (view && view.build === build && view.chr === next.chr &&
                  Math.abs(view.start - next.start) <= 2 && Math.abs(view.end - next.end) <= 2) return;
              send("igv_event", next);
            }, 200);
          });
          if (cfg.genes) {
            try { await instance.loadTrack({name: "Genes", url: cfg.genes, format: "bed", type: "annotation",
              indexed: false, displayMode: "COLLAPSED", height: 55, color: "#356a91"}); }
            catch (error) { status("GRCh" + build + " · 基因注释载入失败"); console.warn(error); }
          }
        } else {
          const current = lociOf(browser);
          if (!current || current.length !== 1 || esc(current[0]).toLowerCase() !== esc(cfg.locus).toLowerCase()) {
            await browser.search(cfg.locus);
          }
        }
        if (!blocked) document.getElementById("igv").classList.remove("busy");
      } catch (error) {
        status("IGV 载入失败：" + error.message); console.error(error);
      } finally {
        await delay(80); applying = false;
      }
    }
    working = false;
  }

  window.bplot = {
    bind: function (el, context) {
      if (el._bplotHandlers) {
        const old = el._bplotHandlers;
        el.removeListener("plotly_relayout", old.range);
        el.removeListener("plotly_click", old.click);
        el.removeEventListener("mousemove", old.move);
        el.removeEventListener("mouseleave", old.leave);
        el.removeEventListener("dblclick", old.double, true);
        el.removeEventListener("mousedown", old.down, true);
        clearTimeout(old.timer); clearTimeout(old.clickTimer);
      }
      const handlers = {};
      const base = {build: context.build, chr: context.chr, page: context.page};
      const tip = document.getElementById("block-hover");
      const hit = event => {
        const layout = el._fullLayout;
        if (!layout || blocked) return null;
        const rect = el.getBoundingClientRect();
        const x = event.clientX - rect.left, y = event.clientY - rect.top;
        const axis = layout.xaxis;
        if (x < axis._offset || x > axis._offset + axis._length) return null;
        let pos = axis.p2d(x - axis._offset), chrom = context.chr;
        if (chrom === "All") {
          const entries = Object.entries(context.offsets).sort((a, b) => a[1] - b[1]);
          for (const [ch, offset] of entries) if (pos >= offset) chrom = ch;
          pos -= context.offsets[chrom];
        }
        const tracks = Array.isArray(context.tracks) ? context.tracks : [context.tracks];
        for (let i = 0; i < tracks.length; i++) {
          const ay = layout[i === 0 ? "yaxis" : "yaxis" + (i + 1)];
          if (!ay || y < ay._offset || y > ay._offset + ay._length) continue;
          const track = tracks[i];
          const block = track.blocks.find(b => String(b[0]) === chrom && pos > b[1] && pos <= b[2]);
          if (block) return {race: track.race, id: block[3], chrom, pos: Math.round(pos)};
        }
        return null;
      };
      handlers.leave = () => { if (tip) tip.style.display = "none"; };
      handlers.move = event => {
        const block = hit(event);
        if (!block || !tip) { handlers.leave(); return; }
        tip.textContent = "block: " + block.race + " " + block.id + " · 双击查看";
        tip.style.display = "block";
        tip.style.left = Math.min(event.clientX + 14, window.innerWidth - 240) + "px";
        tip.style.top = Math.max(8, event.clientY - 35) + "px";
      };
      handlers.down = event => {
        if (event.button !== 0) return;
        const now = performance.now(), last = handlers.lastDown;
        const block = hit(event);
        handlers.lastDown = {time: now, x: event.clientX, y: event.clientY};
        if (!block || !last || now - last.time > 500 || Math.hypot(event.clientX-last.x, event.clientY-last.y) > 6) return;
        handlers.lastDown = null;
        handlers.double(event);
      };
      handlers.double = event => {
        if (performance.now() - (handlers.opened || -1000) < 450) return;
        clearTimeout(handlers.clickTimer); clearTimeout(handlers.timer);
        const block = hit(event);
        if (!block) return;
        handlers.opened = performance.now();
        event.preventDefault(); event.stopPropagation(); handlers.leave();
        send("plot_event", {...base, type: "block", ...block});
      };
      handlers.range = function (event) {
        if (blocked) return;
        const axis = Object.keys(event).find(k => /^xaxis\d*\.(range|autorange)/.test(k));
        if (!axis) return;
        const name = axis.split(".")[0];
        const start = event[name + ".range[0]"] ?? event[name + ".range"]?.[0];
        const end = event[name + ".range[1]"] ?? event[name + ".range"]?.[1];
        const reset = event[name + ".autorange"] === true;
        if (!reset && (!Number.isFinite(start) || !Number.isFinite(end))) return;
        clearTimeout(handlers.timer);
        handlers.timer = setTimeout(() => send("plot_event", {...base, type: "range", start, end, reset}), 180);
      };
      handlers.click = function (event) {
        const point = event.points && event.points[0];
        if (blocked || !point || !point.customdata || performance.now() - (handlers.opened || -1000) < 450) return;
        clearTimeout(handlers.clickTimer);
        handlers.clickTimer = setTimeout(() => send("plot_event", {...base, type: "click", chrom: point.customdata[1], pos: Number(point.customdata[2])}), 550);
      };
      el.on("plotly_relayout", handlers.range);
      el.on("plotly_click", handlers.click);
      el.addEventListener("mousemove", handlers.move);
      el.addEventListener("mouseleave", handlers.leave);
      el.addEventListener("dblclick", handlers.double, true);
      el.addEventListener("mousedown", handlers.down, true);
      el._bplotHandlers = handlers;
    }
  };

  document.addEventListener("click", function (event) {
    const button = event.target.closest(".chr-btn");
    if (button && !blocked) send("chr", button.dataset.chr);
  });
  function register() {
    Shiny.addCustomMessageHandler("bplot_busy", function (value) {
      blocked = value;
      document.getElementById("igv").classList.toggle("busy", value);
      document.querySelectorAll(".chr-btn").forEach(button => { button.disabled = value; });
      if (value) status("正在准备统一坐标的 GWAS 轨道…");
    });
    Shiny.addCustomMessageHandler("bplot_view", function (value) {
      view = value; blocked = false;
      document.querySelectorAll(".chr-btn").forEach(button => {
        button.disabled = false;
        button.classList.toggle("active", button.dataset.chr === value.chr);
        button.setAttribute("aria-pressed", button.dataset.chr === value.chr ? "true" : "false");
      });
    });
    Shiny.addCustomMessageHandler("bplot_igv", function (cfg) { pending = cfg; updateIgv(); });
  }
  if (window.Shiny) register();
  else document.addEventListener("DOMContentLoaded", register, {once: true});
})();
