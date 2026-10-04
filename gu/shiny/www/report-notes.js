(() => {
  'use strict';
  let tip, active, timer;
  function hide() {
    if (active) active.removeAttribute('aria-describedby');
    active = null;
    if (tip) tip.hidden = true;
  }
  function show(mark) {
    clearTimeout(timer);
    if (!tip) {
      tip = document.createElement('div');
      tip.id = 'gu-column-note'; tip.className = 'gu-note-tooltip';
      tip.setAttribute('role', 'tooltip'); document.body.appendChild(tip);
      tip.addEventListener('mouseenter', () => clearTimeout(timer));
      tip.addEventListener('mouseleave', () => { timer = setTimeout(hide, 160); });
    }
    hide(); active = mark; tip.textContent = mark.dataset.guNote;
    mark.setAttribute('aria-describedby', tip.id); tip.hidden = false;
    const r = mark.getBoundingClientRect(), box = tip.getBoundingClientRect();
    tip.style.left = Math.max(12, Math.min(r.left, innerWidth - box.width - 12)) + 'px';
    tip.style.top = Math.max(12, r.bottom + box.height + 12 < innerHeight ? r.bottom + 8 : r.top - box.height - 8) + 'px';
  }
  // Delegation also covers DataTables' cloned scrolling headers and redraws.
  document.addEventListener('mouseover', e => {
    const mark = e.target.closest('[data-gu-note]'); if (mark) show(mark);
  });
  document.addEventListener('mouseout', e => {
    if (e.target.closest('[data-gu-note]')) timer = setTimeout(hide, 160);
  });
  document.addEventListener('focusin', e => {
    const mark = e.target.closest('[data-gu-note]'); if (mark) show(mark);
  });
  document.addEventListener('focusout', e => { if (e.target.closest('[data-gu-note]')) hide(); });
  document.addEventListener('click', e => {
    const mark = e.target.closest('[data-gu-note]');
    if (mark) { e.preventDefault(); e.stopPropagation(); show(mark); } else hide();
  }, true);
  document.addEventListener('keydown', e => {
    if (e.key === 'Escape') hide();
    const mark = e.target.closest('[data-gu-note]');
    if (mark && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); e.stopPropagation(); show(mark); }
  }, true);
  document.addEventListener('scroll', hide, true);
  window.addEventListener('resize', hide);
  // Existing inline metric notes use the same readable treatment.
  document.addEventListener('DOMContentLoaded', () => {
    document.querySelectorAll('.gu-metric-help[title]').forEach(el => {
      el.dataset.guNote = el.title; el.removeAttribute('title'); el.classList.add('gu-note-inline');
    });
  });
})();
