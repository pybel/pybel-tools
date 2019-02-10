require.config({
    paths: {
        Ideogram: '//unpkg.com/ideogram@1.5.0/dist/js/ideogram.min'
    }
});

element.append("<div id='{{ chart }}'></div>");

require(['Ideogram'], function (Ideogram) {
    window.Ideogram = Ideogram.default;
    var ideogram = new window.Ideogram({
        container: '#{{ chart }}',
        organism: 'human',
        annotationsLayout: 'histogram',
        barWidth: 3,
        annotations: {{ annotations|safe }}
    });
});
