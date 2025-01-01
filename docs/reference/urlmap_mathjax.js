// A map between namespaces and base URLs for their online documentation
baseURLs = [
    ['GLib', 'https://docs.gtk.org/glib/'],
    ['GObject', 'https://docs.gtk.org/gobject/'],
    ['Gio', 'https://docs.gtk.org/gio/'],
    ['NumCosmoMath', '../numcosmo-math/'],
]

window.MathJax = {
    tex: { inlineMath: [['$', '$'], ['\\(', '\\)']] },
    svg: { fontCache: 'global' }
};

// Dynamically load MathJax
(function () {
    var script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js';
    script.async = true;
    document.head.appendChild(script);
})();
